# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import gzip
import yaml
import pandas as pd

import numpy as np
from q2_types.per_sample_sequences import (
            SingleLanePerSampleSingleEndFastqDirFmt,
            FastqManifestFormat, YamlFormat, FastqGzFormat)


def _read_fastq_seqs(filepath, phred_offset):
    # This function is adapted from @jairideout's SO post:
    # http://stackoverflow.com/a/39302117/3424666
    fh = gzip.open(filepath, 'rt')
    for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
        qual = qual.strip()
        qual_parsed = np.fromstring(qual, dtype=np.uint8) - phred_offset
        yield (seq_header.strip(), seq.strip(), qual_header.strip(),
               qual, qual_parsed)


def _runs_of_ones(arr):
    """Find the location and length of runs

    This method assumes the input array is boolean
    """
    # http://stackoverflow.com/a/1066838
    # make sure all runs of ones are well-bounded
    bounded = np.hstack(([0], arr, [0]))
    # get 1 at run starts and -1 at run ends
    difs = np.diff(bounded)
    run_starts, = np.where(difs > 0)
    run_ends, = np.where(difs < 0)
    return run_starts, run_ends - run_starts


def _truncate(sequence_record, position):
    """Truncate a record up to a specified position"""
    seq = sequence_record[1][:position]
    qual = sequence_record[3][:position]
    qual_parsed = sequence_record[4][:position]
    return (sequence_record[0], seq, sequence_record[2], qual, qual_parsed)


# defaults as used Bokulich et al, Nature Methods 2013,
# same as QIIME 1.9.1
_default_params = {
    'min_quality': 4,
    'quality_window': 3,
    'min_length_fraction': 0.75,
    'max_ambiguous': 0
}


# TODO: fix up demux fmt writing a la q2-cutadapt
def q_score(demux: SingleLanePerSampleSingleEndFastqDirFmt,
            min_quality: int = _default_params['min_quality'],
            quality_window: int = _default_params['quality_window'],
            min_length_fraction:
            float = _default_params['min_length_fraction'],
            max_ambiguous: int = _default_params['max_ambiguous']) \
                  -> (SingleLanePerSampleSingleEndFastqDirFmt,
                      pd.DataFrame):
    result = SingleLanePerSampleSingleEndFastqDirFmt()

    manifest = FastqManifestFormat()
    manifest_fh = manifest.open()
    manifest_fh.write('sample-id,filename,direction\n')
    manifest_fh.write('# direction is not meaningful in this file as these\n')
    manifest_fh.write('# data may be derived from forward, reverse, or \n')
    manifest_fh.write('# joined reads\n')

    log_records_truncated_counts = {}
    log_records_max_ambig_counts = {}
    log_records_tooshort_counts = {}
    log_records_totalread_counts = {}
    log_records_totalkept_counts = {}

    metadata_view = demux.metadata.view(YamlFormat).open()
    phred_offset = yaml.load(metadata_view)['phred-offset']
    demux_manifest = demux.manifest.view(demux.manifest.format)
    demux_manifest = pd.read_csv(demux_manifest.open(), dtype=str)
    demux_manifest.set_index('filename', inplace=True)

    iterator = demux.sequences.iter_views(FastqGzFormat)
    for bc_id, (fname, fp) in enumerate(iterator):
        sample_id = demux_manifest.loc[str(fname)]['sample-id']

        log_records_truncated_counts[sample_id] = 0
        log_records_max_ambig_counts[sample_id] = 0
        log_records_tooshort_counts[sample_id] = 0
        log_records_totalread_counts[sample_id] = 0
        log_records_totalkept_counts[sample_id] = 0

        # per q2-demux, barcode ID, lane number and read number are not
        # relevant here
        path = result.sequences.path_maker(sample_id=sample_id,
                                           barcode_id=bc_id,
                                           lane_number=1,
                                           read_number=1)

        # we do not open a writer by default in the event that all sequences
        # for a sample are filtered out; an empty fastq file is not a valid
        # fastq file.
        writer = None
        for sequence_record in _read_fastq_seqs(str(fp), phred_offset):
            log_records_totalread_counts[sample_id] += 1

            # determine the length of the runs below quality threshold
            # NOTE: QIIME 1.x used <= the quality threshold and the parameter
            #   -q was interpreted as the maximum unacceptable PHRED score. In
            #   QIIME 2.x, we're now interpreting this as the minimum
            #   acceptable score.
            qual_below_threshold = sequence_record[4] < min_quality
            run_starts, run_lengths = _runs_of_ones(qual_below_threshold)
            bad_windows = np.argwhere(run_lengths > quality_window)

            # if there is a run of sufficient size, truncate it
            if bad_windows.size > 0:
                log_records_truncated_counts[sample_id] += 1

                full_length = len(sequence_record[1])
                sequence_record = _truncate(sequence_record,
                                            run_starts[bad_windows[0]][0])
                trunc_length = len(sequence_record[1])

                # do not keep the read if it is too short following truncation
                if round(trunc_length / full_length, 3) <= min_length_fraction:
                    log_records_tooshort_counts[sample_id] += 1
                    continue

            # do not keep the read if there are too many ambiguous bases
            if sequence_record[1].count('N') > max_ambiguous:
                log_records_max_ambig_counts[sample_id] += 1
                continue

            fastq_lines = '\n'.join(sequence_record[:4]) + '\n'
            fastq_lines = fastq_lines.encode('utf-8')

            if writer is None:
                writer = gzip.open(str(path), mode='w')
            writer.write(fastq_lines)

            log_records_totalkept_counts[sample_id] += 1

        if writer is not None:
            manifest_fh.write('%s,%s,%s\n' % (sample_id, path.name, 'forward'))
            writer.close()

    if set(log_records_totalkept_counts.values()) == {0, }:
        raise ValueError("All sequences from all samples were filtered out. "
                         "The parameter choices may be too stringent for the "
                         "data.")

    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': phred_offset}))
    result.metadata.write_data(metadata, YamlFormat)

    columns = ['sample-id', 'total-input-reads', 'total-retained-reads',
               'reads-truncated',
               'reads-too-short-after-truncation',
               'reads-exceeding-maximum-ambiguous-bases']
    stats = []
    for id_, _ in sorted(log_records_truncated_counts.items()):
        stats.append([id_, log_records_totalread_counts[id_],
                      log_records_totalkept_counts[id_],
                      log_records_truncated_counts[id_],
                      log_records_tooshort_counts[id_],
                      log_records_max_ambig_counts[id_]])

    stats = pd.DataFrame(stats, columns=columns).set_index('sample-id')

    return result, stats
