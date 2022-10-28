# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import qiime2

epoch = qiime2.__release__
m_p_base = f'https://docs.qiime2.org/{epoch}/data/tutorials/moving-pictures/'
demuxed_seqs_url = m_p_base + 'demux.qza'


def q_score_example(use):
    demuxed = use.init_artifact_from_url('demuxed_seqs', demuxed_seqs_url)

    filtered_seqs, stats = use.action(
        use.UsageAction('quality_filter', 'q_score'),
        use.UsageInputs(
            demux=demuxed
        ),
        use.UsageOutputNames(
            filtered_sequences='dumux-filtered',
            filter_stats='demux_filter_stats'
        )
    )

    filtered_seqs.assert_output_type('SampleData[SequencesWithQuality]')
    stats.assert_output_type('QualityFilterStats')
