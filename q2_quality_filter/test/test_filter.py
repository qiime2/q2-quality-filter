# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import gzip

import pandas as pd
import pandas.util.testing as pdt
from qiime2.sdk import Artifact
import numpy as np
import numpy.testing as npt
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
        FastqGzFormat,
        SingleLanePerSampleSingleEndFastqDirFmt)

from q2_quality_filter._filter import (_read_fastq_seqs, _runs_of_ones,
                                       _truncate, q_score)


class FilterTests(TestPluginBase):
    package = 'q2_quality_filter.test'

    def test_read_fastq_seqs(self):
        exp = [('@foo', 'ATGC', '+', 'IIII', np.array([40, 40, 40, 40])),
               ('@bar', 'TGCA', '+', 'ABCD', np.array([32, 33, 34, 35]))]
        obs = list(_read_fastq_seqs(self.get_data_path('simple.fastq.gz'), 33))
        self.assertEqual(len(obs), 2)

        for o, e in zip(obs, exp):
            self.assertEqual(o[:4], e[:4])
            npt.assert_equal(o[4], e[4])

    def test_runs_of_ones(self):
        data = [np.array([0, 0, 0, 0, 0, 0], dtype=bool),
                np.array([1, 0, 1, 0, 1, 0], dtype=bool),
                np.array([1, 1, 1, 1, 1, 1], dtype=bool),
                np.array([0, 1, 1, 1, 0, 0], dtype=bool),
                np.array([0, 0, 0, 0, 0, 1], dtype=bool)]

        exp_starts = [np.array([]), np.array([0, 2, 4]), np.array([0]),
                      np.array([1]), np.array([5])]
        exp_lengths = [np.array([]), np.array([1, 1, 1]), np.array([6]),
                       np.array([3]), np.array([1])]

        for i, d in enumerate(data):
            o_starts, o_lengths = _runs_of_ones(d)
            npt.assert_equal(o_starts, exp_starts[i])
            npt.assert_equal(o_lengths, exp_lengths[i])

    def test_truncate(self):
        data = [('@x', 'ATGCG', '+', 'IIIIA', np.array([40, 40, 40, 40, 32])),
                ('@y', 'TGCAC', '+', 'ABCDA', np.array([32, 33, 34, 35, 32]))]

        exp1 = [('@x', 'A', '+', 'I', np.array([40])),
                ('@y', 'T', '+', 'A', np.array([32]))]

        exp2 = [('@x', 'AT', '+', 'II', np.array([40, 40])),
                ('@y', 'TG', '+', 'AB', np.array([32, 33]))]

        for i, d in enumerate(data):
            o1 = _truncate(d, 1)
            o2 = _truncate(d, 2)
            self.assertEqual(o1[:4], exp1[i][:4])
            npt.assert_equal(o1[4], exp1[i][4])
            self.assertEqual(o2[:4], exp2[i][:4])
            npt.assert_equal(o2[4], exp2[i][4])

    def test_q_score(self):
        ar = Artifact.load(self.get_data_path('simple.qza'))
        view = ar.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_drop_ambig, stats = q_score(view, quality_window=2,
                                        min_length_fraction=0.25)

        exp_drop_ambig = ["@foo_1",
                          "ATGCATGC",
                          "+",
                          "DDDDBBDD"]
        columns = ['sample-id', 'total-input-reads', 'total-retained-reads',
                   'reads-truncated',
                   'reads-too-short-after-truncation',
                   'reads-exceeding-maximum-ambiguous-bases']
        exp_drop_ambig_stats = pd.DataFrame([('foo', 2, 1, 0, 0, 1),
                                             ('bar', 1, 0, 0, 0, 1)],
                                            columns=columns)
        exp_drop_ambig_stats = exp_drop_ambig_stats.set_index('sample-id')
        obs = []
        iterator = obs_drop_ambig.sequences.iter_views(FastqGzFormat)
        for sample_id, fp in iterator:
            obs.extend([l.strip() for l in gzip.open(str(fp), 'rt')])
        self.assertEqual(obs, exp_drop_ambig)
        pdt.assert_frame_equal(stats, exp_drop_ambig_stats.loc[stats.index])

        obs_trunc, stats = q_score(view, quality_window=1, min_quality=33,
                                   min_length_fraction=0.25)
        exp_trunc = ["@foo_1",
                     "ATGCATGC",
                     "+",
                     "DDDDBBDD",
                     "@bar_1",
                     "ATA",
                     "+",
                     "DDD"]
        exp_trunc_stats = pd.DataFrame([('foo', 2, 1, 0, 0, 1),
                                        ('bar', 1, 1, 1, 0, 0)],
                                       columns=columns)
        exp_trunc_stats = exp_trunc_stats.set_index('sample-id')

        obs = []
        for sample_id, fp in obs_trunc.sequences.iter_views(FastqGzFormat):
            obs.extend([l.strip() for l in gzip.open(str(fp), 'rt')])
        self.assertEqual(sorted(obs), sorted(exp_trunc))
        pdt.assert_frame_equal(stats, exp_trunc_stats.loc[stats.index])

    def test_q_score_real(self):
        ar = Artifact.load(self.get_data_path('real_data.qza'))
        view = ar.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_result, stats = q_score(view, min_quality=40,
                                    min_length_fraction=0.24)

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
            obs.extend([l.strip() for l in gzip.open(str(fp), 'rt')])
        self.assertEqual(obs, exp_result)
        pdt.assert_frame_equal(stats, exp_stats.loc[stats.index])


if __name__ == '__main__':
    unittest.main()
