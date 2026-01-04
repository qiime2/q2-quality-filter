# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt

from copy import copy
import gzip
import os
from pathlib import Path
import tempfile
import unittest

import qiime2
from qiime2.sdk import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    FastqGzFormat,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    PairedEndSequencesWithQuality,
)

from q2_quality_filter._filter import (
    FastqRecord,
    _read_fastq_records,
    _find_low_quality_window,
    _truncate,
    RecordStatus,
    _process_record,
    _is_retained,
    _write_record,
)
from q2_quality_filter._format import QualityFilterStatsFmt


class HelperMethodTests(TestPluginBase):
    package = 'q2_quality_filter.tests'

    def test_read_fastq_records(self):
        exp = [
            FastqRecord(b'@foo', b'ATGC', b'+', b'IIII'),
            FastqRecord(b'@bar', b'TGCA', b'+', b'ABCD')
        ]

        obs = list(
            _read_fastq_records(self.get_data_path('simple.fastq.gz'))
        )

        self.assertEqual(len(obs), 2)

        attrs = [
            'sequence_header', 'sequence', 'quality_header', 'quality_scores'
        ]
        for exp_record, obs_record in zip(exp, obs):
            for attr in attrs:
                self.assertEqual(
                    exp_record.__getattribute__(attr),
                    obs_record.__getattribute__(attr)
                )

    def test_find_low_quality_window(self):
        # test no low quality window returns none
        # 'M' has quality score of 44 with PHRED offset of 33
        quality_scores = b'M' * 10
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=20, window_length=2
        )
        self.assertEqual(obs, None)

        # test that `min_quality` bases are not considered part of a window
        # (only scores that are lower)
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=44, window_length=2
        )
        self.assertEqual(obs, None)

        # test windows detected correctly
        # quality scores: M => 44; + => 10
        quality_scores = b'MMM++MM'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=15, window_length=2
        )
        self.assertEqual(obs, 3)

        quality_scores = b'M++MM+++MM'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=15, window_length=3
        )
        self.assertEqual(obs, 5)

        quality_scores = b'++MMMM'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=11, window_length=2
        )
        self.assertEqual(obs, 0)

        quality_scores = b'M++MMMM+++'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=11, window_length=3
        )
        self.assertEqual(obs, 7)

        # low quality window is detected when it's longer than window length
        quality_scores = b'MMMMMMM+++'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=11, window_length=2
        )
        self.assertEqual(obs, 7)

        # test when multiple windows exist, first window is returned
        quality_scores = b'ML++MMM+++'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=20, window_length=2
        )
        self.assertEqual(obs, 2)

        quality_scores = b'++ML+++M+++MM++'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=20, window_length=3
        )
        self.assertEqual(obs, 4)
        # test that when all windows are too short, None is returned
        quality_scores = b'++ML+++M+++MM++'
        obs = _find_low_quality_window(
            quality_scores, phred_offset=33, min_quality=20, window_length=4
        )
        self.assertEqual(obs, None)

    def test_truncate(self):
        fastq_record = FastqRecord(
            b'@header', b'ATTCTGTA', b'+', b'MMLMLL++'
        )

        truncated = _truncate(copy(fastq_record), position=4)
        exp = FastqRecord(
            b'@header', b'ATTC', b'+', b'MMLM'
        )
        self.assertEqual(truncated, exp)

        truncated = _truncate(copy(fastq_record), position=7)
        exp = FastqRecord(
            b'@header', b'ATTCTGT', b'+', b'MMLMLL+'
        )
        self.assertEqual(truncated, exp)

        truncated = _truncate(copy(fastq_record), position=1)
        exp = FastqRecord(
            b'@header', b'A', b'+', b'M'
        )
        self.assertEqual(truncated, exp)

        truncated = _truncate(copy(fastq_record), position=0)
        exp = FastqRecord(
            b'@header', b'', b'+', b''
        )
        self.assertEqual(truncated, exp)

    def test_process_record(self):
        # truncation
        fastq_record = FastqRecord(
            b'@header', b'ATTCTGTA', b'+', b'MMLMLL++'
        )
        processed_record, status = _process_record(
            copy(fastq_record),
            phred_offset=33,
            min_quality=15,
            window_length=2,
            min_length_fraction=0.5,
            max_ambiguous=0
        )
        exp_record = FastqRecord(
            b'@header', b'ATTCTG', b'+', b'MMLMLL'
        )
        exp_status = RecordStatus.TRUNCATED
        self.assertEqual(processed_record, exp_record)
        self.assertEqual(status, exp_status)

        # no truncation
        processed_record, status = _process_record(
            copy(fastq_record),
            phred_offset=33,
            min_quality=5,
            window_length=2,
            min_length_fraction=0.5,
            max_ambiguous=0
        )
        exp_record = fastq_record
        exp_status = RecordStatus.UNTRUNCATED
        self.assertEqual(processed_record, exp_record)
        self.assertEqual(status, exp_status)

        # ambiguous
        fastq_record = FastqRecord(
            b'@header', b'ATTCTNTN', b'+', b'MMLMLL++'
        )
        processed_record, status = _process_record(
            copy(fastq_record),
            phred_offset=33,
            min_quality=5,
            window_length=2,
            min_length_fraction=0.5,
            max_ambiguous=1
        )
        exp_record = FastqRecord(
            b'@header', b'ATTCTNTN', b'+', b'MMLMLL++'
        )
        exp_status = RecordStatus.AMBIGUOUS
        self.assertEqual(processed_record, exp_record)
        self.assertEqual(status, exp_status)

        # truncation and ambiguous
        fastq_record = FastqRecord(
            b'@header', b'ATTCTNTA', b'+', b'MMLMLL++'
        )
        processed_record, status = _process_record(
            copy(fastq_record),
            phred_offset=33,
            min_quality=15,
            window_length=2,
            min_length_fraction=0.5,
            max_ambiguous=0
        )
        exp_record = FastqRecord(
            b'@header', b'ATTCTN', b'+', b'MMLMLL'
        )
        exp_status = RecordStatus.TRUNCATED_AMBIGUOUS
        self.assertEqual(processed_record, exp_record)
        self.assertEqual(status, exp_status)

        # truncation and too short
        fastq_record = FastqRecord(
            b'@header', b'ATTCTGTA', b'+', b'MMLMLL++'
        )
        processed_record, status = _process_record(
            copy(fastq_record),
            phred_offset=33,
            min_quality=15,
            window_length=2,
            min_length_fraction=0.9,
            max_ambiguous=0
        )
        exp_record = FastqRecord(
            b'@header', b'ATTCTG', b'+', b'MMLMLL'
        )
        exp_status = RecordStatus.SHORT
        self.assertEqual(processed_record, exp_record)
        self.assertEqual(status, exp_status)

    def test_is_retained(self):
        filtering_stats = pd.Series(
            data=0,
            name='sample-a',
            index=[
                'total-input-reads',
                'total-retained-reads',
                'reads-truncated',
                'reads-too-short-after-truncation',
                'reads-exceeding-maximum-ambiguous-bases'
            ]
        )

        # retained and truncated
        retained = _is_retained(
            forward_status=RecordStatus.TRUNCATED,
            reverse_status=RecordStatus.UNTRUNCATED,
            filtering_stats=filtering_stats,
        )
        self.assertTrue(retained)
        self.assertEqual(filtering_stats['total-retained-reads'], 1)
        self.assertEqual(filtering_stats['reads-truncated'], 1)
        filtering_stats[:] = 0

        # forward read only, retained
        retained = _is_retained(
            forward_status=RecordStatus.TRUNCATED,
            reverse_status=None,
            filtering_stats=filtering_stats,
        )
        self.assertTrue(retained)
        self.assertEqual(filtering_stats['total-retained-reads'], 1)
        self.assertEqual(filtering_stats['reads-truncated'], 1)
        self.assertEqual(
            filtering_stats['reads-too-short-after-truncation'], 0
        )
        filtering_stats[:] = 0

        # forward read only, short
        retained = _is_retained(
            forward_status=RecordStatus.SHORT,
            reverse_status=None,
            filtering_stats=filtering_stats,
        )
        self.assertFalse(retained)
        self.assertEqual(filtering_stats['total-retained-reads'], 0)
        self.assertEqual(filtering_stats['reads-truncated'], 1)
        self.assertEqual(
            filtering_stats['reads-too-short-after-truncation'], 1
        )
        filtering_stats[:] = 0

        # one read untruncated, one read truncated and ambiguous
        retained = _is_retained(
            forward_status=RecordStatus.UNTRUNCATED,
            reverse_status=RecordStatus.TRUNCATED_AMBIGUOUS,
            filtering_stats=filtering_stats,
        )
        self.assertFalse(retained)
        self.assertEqual(filtering_stats['total-retained-reads'], 0)
        self.assertEqual(
            filtering_stats['reads-exceeding-maximum-ambiguous-bases'], 1
        )
        self.assertEqual(filtering_stats['reads-truncated'], 1)
        filtering_stats[:] = 0

    def test_write_record(self):
        fastq_record = FastqRecord(
            b'@header', b'ATTCTGTA', b'+', b'MMLMLL++'
        )

        with tempfile.TemporaryDirectory() as tempdir:
            fp = Path(tempdir) / 'file.fastq.gz'

            with gzip.open(fp, 'wb') as fh:
                _write_record(fastq_record, fh)

            with gzip.open(fp, 'rb') as fh:
                contents = fh.read()
                exp = b'@header\nATTCTGTA\n+\nMMLMLL++\n'
                self.assertEqual(contents, exp)

            with gzip.open(fp, 'ab') as fh:
                _write_record(fastq_record, fh)
                _write_record(fastq_record, fh)

            with gzip.open(fp, 'rb') as fh:
                contents = fh.read()
                exp = b'@header\nATTCTGTA\n+\nMMLMLL++\n' * 3
                self.assertEqual(contents, exp)


class QScoreSingleEndTests(TestPluginBase):
    package = 'q2_quality_filter.tests'

    def test_q_score_all_dropped(self):
        ar = Artifact.load(self.get_data_path('simple.qza'))

        with self.assertRaisesRegex(
            ValueError, 'All sequences from all samples were filtered'
        ):
            with redirected_stdio(stdout=os.devnull):
                self.plugin.methods['q_score'](ar, min_quality=50)

    def test_q_score_numeric_ids(self):
        ar = Artifact.load(self.get_data_path('numeric_ids.qza'))
        exp_sids = {'00123', '0.4560'}

        with redirected_stdio(stdout=os.devnull):
            obs_ar, stats_ar = self.plugin.methods['q_score'](
                ar, min_quality=2)
        obs = obs_ar.view(SingleLanePerSampleSingleEndFastqDirFmt)
        stats = stats_ar.view(pd.DataFrame)
        obs_manifest = obs.manifest.view(obs.manifest.format)
        obs_manifest = pd.read_csv(obs_manifest.open(), dtype=str, comment='#')
        obs_manifest.set_index('sample-id', inplace=True)

        obs_sids = set(obs_manifest.index)
        self.assertEqual(obs_sids, exp_sids)
        self.assertEqual(set(stats.index), exp_sids)

    def test_q_score_different_num_processes(self):
        ar = Artifact.load(self.get_data_path('simple.qza'))

        for num_processes in (1, 2):
            with redirected_stdio(stdout=os.devnull):
                obs_drop_ambig_ar, stats_ar = self.plugin.methods['q_score'](
                    ar,
                    quality_window=2,
                    min_quality=20,
                    min_length_fraction=0.25,
                    num_processes=num_processes,
                )
            obs_drop_ambig = obs_drop_ambig_ar.view(
                SingleLanePerSampleSingleEndFastqDirFmt)
            stats = stats_ar.view(pd.DataFrame)

            exp_drop_ambig = ["@foo_1", "ATGCATGC", "+", "DDDDBBDD"]
            columns = [
                'sample-id',
                'total-input-reads',
                'total-retained-reads',
                'reads-truncated',
                'reads-too-short-after-truncation',
                'reads-exceeding-maximum-ambiguous-bases'
            ]
            exp_drop_ambig_stats = pd.DataFrame(
                [('foo', 2, 1, 0, 0, 1), ('bar', 1, 0, 0, 0, 1)],
                columns=columns
            )
            exp_drop_ambig_stats = exp_drop_ambig_stats.set_index('sample-id')
            obs = []
            iterator = obs_drop_ambig.sequences.iter_views(FastqGzFormat)
            for sample_id, fp in iterator:
                obs.extend([x.strip() for x in gzip.open(str(fp), 'rt')])
            self.assertEqual(obs, exp_drop_ambig)
            pdt.assert_frame_equal(
                stats, exp_drop_ambig_stats.loc[stats.index]
            )

            with redirected_stdio(stdout=os.devnull):
                obs_trunc_ar, stats_ar = self.plugin.methods['q_score'](
                    ar,
                    quality_window=1,
                    min_quality=33,
                    min_length_fraction=0.25,
                    num_processes=num_processes,
                )
            obs_trunc = obs_trunc_ar.view(
                SingleLanePerSampleSingleEndFastqDirFmt
            )
            stats = stats_ar.view(pd.DataFrame)

            exp_trunc = [
                "@foo_1", "ATGCATGC", "+", "DDDDBBDD",
                "@bar_1", "ATA", "+", "DDD"
            ]
            exp_trunc_stats = pd.DataFrame(
                [('foo', 2, 1, 0, 0, 1), ('bar', 1, 1, 1, 0, 0)],
                columns=columns
            )
            exp_trunc_stats = exp_trunc_stats.set_index('sample-id')

            obs = []
            for sample_id, fp in obs_trunc.sequences.iter_views(FastqGzFormat):
                obs.extend([x.strip() for x in gzip.open(str(fp), 'rt')])
            self.assertEqual(sorted(obs), sorted(exp_trunc))
            pdt.assert_frame_equal(stats, exp_trunc_stats.loc[stats.index])

    def test_q_score_real(self):
        self.maxDiff = None

        ar = Artifact.load(self.get_data_path('real_data.qza'))
        with redirected_stdio(stdout=os.devnull):
            obs_ar, stats_ar = self.plugin.methods['q_score'](
                ar, min_quality=40, min_length_fraction=0.24)
        obs_result = obs_ar.view(SingleLanePerSampleSingleEndFastqDirFmt)
        stats = stats_ar.view(pd.DataFrame)

        # All input reads are represented here in their post-quality filtered
        # form. Reads that are commented out were manually identified as being
        # filtered by the q_score method. For the commented reads, the comments
        # denote why the read is not retained.

        # The first read, @HWI-EAS440_0386:1:32:15467:1432#0/1, is 25% of
        # total read length and is indicative of a sequence at the
        # min_length_fraction boundary.
        exp_result = [
                      "@HWI-EAS440_0386:1:32:15467:1432#0/1",
                      "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTT",
                      "+",
                      "hhhhhhhhhhhhfghhhghghghhhchhhahhhhhfhh",

                      # too short
                      # "@HWI-EAS440_0386:1:36:9986:17043#0/1",
                      # "TACGTAGGTGGCAAGCGTTATCCGGATTTATTG",
                      # "+",
                      # "hhhhhhhhhhhhhhhhhhhhhhhhhffhhghhh",

                      "@HWI-EAS440_0386:1:37:13343:14820#0/1",
                      ("TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGAT"
                       "GGATGTTTAAGTCAGTTGTG"),
                      "+",
                      ("hhhhhhhhhhhhhfhhhhhfhhhhghhhhghhhhhhhhhgghhhgghhhgghh"
                       "hgdhhhhghghhhdhhhhgh"),

                      "@HWI-EAS440_0386:1:41:18215:15404#0/1",
                      "TACGTAGGTGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCATGTA",
                      "+",
                      "hhhhhhhhhhhhghhhhhhhhhhhhffhhghhhhghhghgghghhhhhgh",

                      # too short
                      # "@HWI-EAS440_0386:1:42:5423:19606#0/1",
                      # "TACGTAGGGAGCAAGCGTT",
                      # "+",
                      # "hhhhghhhhhhhhhghhfh",

                      "@HWI-EAS440_0386:1:52:7507:5841#0/1",
                      "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTT",
                      "+",
                      "hhhhhhhhhghhfghhhhhhhhhhgfhhhghhhghdhh",

                      "@HWI-EAS440_0386:1:53:18599:4074#0/1",
                      "TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTG",
                      "+",
                      "hhhhfhhhhhfhhhhhhfhffhghhfgghggghdcbh",

                      # too short
                      # "@HWI-EAS440_0386:1:55:16425:9514#0/1",
                      # "TACGGAGGATCCGAGCGTTATCCGGATT",
                      # "+",
                      # "hhhhhhhhhhhhfghhhghghhhhhbgh",

                      "@HWI-EAS440_0386:1:65:12049:5619#0/1",
                      "TACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTG",
                      "+",
                      "hhhhhhhhhhhhhhhhhhhhhhhhhfhhhhhhhghdhghhhhhghcfh",

                      # @HWI-EAS440_0386:1:95:4837:16388#0/1
                      # starts off < Q40
                      ]

        columns = ['sample-id', 'total-input-reads', 'total-retained-reads',
                   'reads-truncated',
                   'reads-too-short-after-truncation',
                   'reads-exceeding-maximum-ambiguous-bases']
        exp_stats = pd.DataFrame([('foo', 10, 6, 10, 4, 0)],
                                 columns=columns)
        exp_stats = exp_stats.set_index('sample-id')
        obs = []
        iterator = obs_result.sequences.iter_views(FastqGzFormat)
        for sample_id, fp in iterator:
            obs.extend([x.strip() for x in gzip.open(str(fp), 'rt')])
        self.assertEqual(obs, exp_result)
        pdt.assert_frame_equal(stats, exp_stats.loc[stats.index])

    def test_q_score_real_joined(self):
        ar = Artifact.load(self.get_data_path('real_data_joined.qza'))
        with redirected_stdio(stdout=os.devnull):
            obs_ar, stats_ar = self.plugin.methods['q_score'](
                ar, min_quality=40, min_length_fraction=0.24)
        obs_result = obs_ar.view(SingleLanePerSampleSingleEndFastqDirFmt)
        stats = stats_ar.view(pd.DataFrame)

        # All input reads are represented here in their post-quality filtered
        # form. Reads that are commented out were manually identified as being
        # filtered by the q_score method. For the commented reads, the comments
        # denote why the read is not retained.

        # The first read, @HWI-EAS440_0386:1:32:15467:1432#0/1, is 25% of
        # total read length and is indicative of a sequence at the
        # min_length_fraction boundary.
        exp_result = [
                      "@HWI-EAS440_0386:1:32:15467:1432#0/1",
                      "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTT",
                      "+",
                      "hhhhhhhhhhhhfghhhghghghhhchhhahhhhhfhh",

                      # too short
                      # "@HWI-EAS440_0386:1:36:9986:17043#0/1",
                      # "TACGTAGGTGGCAAGCGTTATCCGGATTTATTG",
                      # "+",
                      # "hhhhhhhhhhhhhhhhhhhhhhhhhffhhghhh",

                      "@HWI-EAS440_0386:1:37:13343:14820#0/1",
                      "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGAT"
                      "GGATGTTTAAGTCAGTTGTG",
                      "+",
                      "hhhhhhhhhhhhhfhhhhhfhhhhghhhhghhhhhhhhhgghhhgghhhgghh"
                      "hgdhhhhghghhhdhhhhgh",

                      "@HWI-EAS440_0386:1:41:18215:15404#0/1",
                      "TACGTAGGTGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCATGTA",
                      "+",
                      "hhhhhhhhhhhhghhhhhhhhhhhhffhhghhhhghhghgghghhhhhgh",

                      # too short
                      # "@HWI-EAS440_0386:1:42:5423:19606#0/1",
                      # "TACGTAGGGAGCAAGCGTT",
                      # "+",
                      # "hhhhghhhhhhhhhghhfh",

                      "@HWI-EAS440_0386:1:52:7507:5841#0/1",
                      "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTT",
                      "+",
                      "hhhhhhhhhghhfghhhhhhhhhhgfhhhghhhghdhh",

                      "@HWI-EAS440_0386:1:53:18599:4074#0/1",
                      "TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTG",
                      "+",
                      "hhhhfhhhhhfhhhhhhfhffhghhfgghggghdcbh",

                      # too short
                      # "@HWI-EAS440_0386:1:55:16425:9514#0/1",
                      # "TACGGAGGATCCGAGCGTTATCCGGATT",
                      # "+",
                      # "hhhhhhhhhhhhfghhhghghhhhhbgh",

                      "@HWI-EAS440_0386:1:65:12049:5619#0/1",
                      "TACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTG",
                      "+",
                      "hhhhhhhhhhhhhhhhhhhhhhhhhfhhhhhhhghdhghhhhhghcfh",

                      # @HWI-EAS440_0386:1:95:4837:16388#0/1
                      # starts off < Q40
                      ]

        columns = ['sample-id', 'total-input-reads', 'total-retained-reads',
                   'reads-truncated',
                   'reads-too-short-after-truncation',
                   'reads-exceeding-maximum-ambiguous-bases']
        exp_stats = pd.DataFrame([('foo', 10, 6, 10, 4, 0)],
                                 columns=columns)
        exp_stats = exp_stats.set_index('sample-id')
        obs = []
        iterator = obs_result.sequences.iter_views(FastqGzFormat)
        for sample_id, fp in iterator:
            obs.extend([x.strip() for x in gzip.open(str(fp), 'rt')])
        self.assertEqual(obs, exp_result)
        pdt.assert_frame_equal(stats, exp_stats.loc[stats.index])


class QScorePairedEndTests(TestPluginBase):
    package = 'q2_quality_filter.tests'

    def _get_header_diff(
        self, forward_record: FastqRecord, reverse_record: FastqRecord
    ) -> int:
        zipped_headers = zip(
            forward_record.sequence_header, reverse_record.sequence_header
        )

        diff = 0
        for forward_byte, reverse_byte in zipped_headers:
            if forward_byte != reverse_byte:
                diff += 1

        return diff

    def _assert_records_match(self, manifest_df: pd.DataFrame):
        for forward_fp, reverse_fp in zip(
            manifest_df['forward'], manifest_df['reverse']
        ):
            forward_iterator = _read_fastq_records(forward_fp)
            reverse_iterator = _read_fastq_records(reverse_fp)
            iterator = zip(forward_iterator, reverse_iterator)

            for forward_record, reverse_record in iterator:
                # headers differ in one position to indicate read direction
                self.assertEqual(
                    self._get_header_diff(forward_record, reverse_record), 1
                )

    def test_paired_end_sequences(self):
        demux_artifact = Artifact.import_data(
            SampleData[PairedEndSequencesWithQuality],
            self.get_data_path('paired-end-data'),
        )

        output_seqs, stats = self.plugin.methods['q_score'](
            demux_artifact,
            min_quality=15,
            quality_window=2,
            min_length_fraction=0.8,
            max_ambiguous=2
        )
        output_demux_format = output_seqs.view(
            SingleLanePerSamplePairedEndFastqDirFmt
        )
        demux_manifest_df = output_demux_format.manifest.view(pd.DataFrame)

        # corresponding records should have matching headers
        self._assert_records_match(demux_manifest_df)

        # "Human-Kneecap2_S2" is dropped because the R1 reads have low q scores
        exp_sample_ids = ['Human-Kneecap', 'Human-Kneecap3']
        self.assertEqual(
            set(demux_manifest_df.index), set(exp_sample_ids)
        )
        self.assertEqual(len(demux_manifest_df), 2)

        # assert truncation positions are correct
        sample1_forward_exp = [
            # first record dropped because of R2 scores
            b'@M00899:113:000000000-A5K20:1:1101:25454:3578 1:N:0:2',
            b'CCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTA',
            b'+',
            b'8ACCCGD@AA=18=======;CEFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGF',
            b'@M00899:113:000000000-A5K20:1:1101:25177:3605 1:N:0:2',
            b'CCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGACGGAAGTCTGAACCAGCCAAGTAGCGTGCAGGATGAC', # noqa
            b'+',
            b'88BCCEDAD9018======;;CCFGGGGFGGGFGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGG', # noqa
        ]

        with gzip.open(
            demux_manifest_df.loc['Human-Kneecap', 'forward']
        ) as fh:
            sample1_forward_obs = [line.strip() for line in fh.readlines()]

        self.assertEqual(sample1_forward_exp, sample1_forward_obs)

        sample1_reverse_exp = [
            # first record dropped because of R2 scores
            b'@M00899:113:000000000-A5K20:1:1101:25454:3578 2:N:0:2',
            b'GACTACCGGGGTATCTAATCCTGTTCGATACCCGCACCTTCGAGCTTCAGCGTCAGTTGCGCTCCCGTCAGCTGC', # noqa
            b'+',
            b'CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG', # noqa
            b'@M00899:113:000000000-A5K20:1:1101:25177:3605 2:N:0:2',
            b'GACTACTGGGGTATCTAATCCTGTTTGATACCCGCACCTTCGAGCTTAAGCGTCAGTTGCGCTCCCGTCAGCTGC', # noqa
            b'+',
            b'CCCCCGG9FFGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG', # noqa
        ]

        with gzip.open(
            demux_manifest_df.loc['Human-Kneecap', 'reverse']
        ) as fh:
            sample1_reverse_obs = [line.strip() for line in fh.readlines()]

        self.assertEqual(sample1_reverse_exp, sample1_reverse_obs)


class TransformerTests(TestPluginBase):
    package = 'q2_quality_filter.tests'

    def test_stats_to_metadata(self):
        filepath = self.get_data_path('stats-1.txt')
        format = QualityFilterStatsFmt(filepath, mode='r')
        transformer = self.get_transformer(QualityFilterStatsFmt,
                                           qiime2.Metadata)
        obs = transformer(format)
        self.assertEqual(obs.id_count, 34)
        self.assertEqual(obs.column_count, 5)
        self.assertEqual(obs.id_header, 'sample-id')

    def test_numeric_ids(self):
        filepath = self.get_data_path('stats-numeric.txt')
        format = QualityFilterStatsFmt(filepath, mode='r')
        transformer = self.get_transformer(QualityFilterStatsFmt,
                                           qiime2.Metadata)
        obs = transformer(format)
        self.assertEqual(obs.id_count, 34)
        self.assertEqual(obs.column_count, 5)
        self.assertEqual(obs.id_header, 'sample-id')


class TestUsageExamples(TestPluginBase):
    package = 'q2_quality_filter.tests'

    def test_examples(self):
        self.execute_examples()


if __name__ == '__main__':
    unittest.main()
