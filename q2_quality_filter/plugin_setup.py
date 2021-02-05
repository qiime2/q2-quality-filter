# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

import qiime2.plugin
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality)

import q2_quality_filter
from q2_quality_filter._type import QualityFilterStats
from q2_quality_filter._format import (QualityFilterStatsFmt,
                                       QualityFilterStatsDirFmt)

citations = qiime2.plugin.Citations.load(
    'citations.bib', package='q2_quality_filter')
plugin = qiime2.plugin.Plugin(
    name='quality-filter',
    version=q2_quality_filter.__version__,
    website='https://github.com/qiime2/q2-quality-filter',
    package='q2_quality_filter',
    user_support_text=None,
    description=('This QIIME 2 plugin supports filtering and trimming of '
                 'sequence reads based on PHRED scores and ambiguous '
                 'nucleotide characters.'),
    short_description='Plugin for PHRED-based filtering and trimming.',
    citations=citations
)

plugin.register_formats(QualityFilterStatsFmt, QualityFilterStatsDirFmt)

plugin.register_semantic_types(QualityFilterStats)
plugin.register_semantic_type_to_format(
    QualityFilterStats,
    artifact_format=QualityFilterStatsDirFmt)

InputMap, OutputMap = qiime2.plugin.TypeMap({
    SampleData[SequencesWithQuality | PairedEndSequencesWithQuality]:
        SampleData[SequencesWithQuality],

    SampleData[JoinedSequencesWithQuality]:
        SampleData[JoinedSequencesWithQuality],
})

_q_score_parameters = {
    'min_quality': qiime2.plugin.Int,
    'quality_window': qiime2.plugin.Int,
    'min_length_fraction': qiime2.plugin.Float,
    'max_ambiguous': qiime2.plugin.Int
}

_q_score_input_descriptions = {
    'demux': 'The demultiplexed sequence data to be quality filtered.'
}

_q_score_parameter_descriptions = {
    'min_quality': 'The minimum acceptable PHRED score. All PHRED scores '
                   'less that this value are considered to be low PHRED '
                   'scores.',
    'quality_window': 'The maximum number of low PHRED scores that '
                      'can be observed in direct succession before '
                      'truncating a sequence read.',
    'min_length_fraction': 'The minimum length that a sequence read can '
                           'be following truncation and still be '
                           'retained. This length should be provided '
                           'as a fraction of the input sequence length.',
    'max_ambiguous': 'The maximum number of ambiguous (i.e., N) base '
                     'calls. This is applied after trimming sequences '
                     'based on `min_length_fraction`.'
}

_q_score_output_descriptions = {
    'filtered_sequences': 'The resulting quality-filtered sequences.',
    'filter_stats': 'Summary statistics of the filtering process.'
}


plugin.methods.register_function(
    function=q2_quality_filter.q_score,
    inputs={'demux': InputMap},
    parameters=_q_score_parameters,
    outputs=[
        ('filtered_sequences', OutputMap),
        ('filter_stats', QualityFilterStats)
    ],
    input_descriptions=_q_score_input_descriptions,
    parameter_descriptions=_q_score_parameter_descriptions,
    output_descriptions=_q_score_output_descriptions,
    name='Quality filter based on sequence quality scores.',
    description=('This method filters sequence based on quality scores and '
                 'the presence of ambiguous base calls.')
)

importlib.import_module('q2_quality_filter._transformer')
