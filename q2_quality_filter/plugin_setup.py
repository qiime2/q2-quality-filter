# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
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
import q2_quality_filter._examples as ex

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
    SampleData[SequencesWithQuality]:
        SampleData[SequencesWithQuality],
    SampleData[PairedEndSequencesWithQuality]:
        SampleData[PairedEndSequencesWithQuality],
    SampleData[JoinedSequencesWithQuality]:
        SampleData[JoinedSequencesWithQuality],
})

_q_score_parameters = {
    'min_quality': qiime2.plugin.Int,
    'quality_window': qiime2.plugin.Int,
    'min_length_fraction': qiime2.plugin.Float,
    'max_ambiguous': qiime2.plugin.Int,
    'num_processes': qiime2.plugin.Threads,
}

_q_score_input_descriptions = {
    'demux': 'The demultiplexed sequence data to be quality filtered.'
}

_q_score_parameter_descriptions = {
    'min_quality': (
        'The minimum acceptable PHRED score. All PHRED scores less that this '
        'value are considered to be low PHRED scores.'
    ),
    'quality_window': (
        'The maximum number of low PHRED scores that can be observed in '
        'direct succession before truncating a sequence read. Note that '
        'truncation is performed such that the entirety of the low quality '
        'window is truncated.'
    ),
    'min_length_fraction': (
        'Filter truncated reads whose length fraction (truncated length '
        'divided by original length) is less than or equal to this value.'
    ),
    'max_ambiguous': (
        'The maximum number of ambiguous (i.e., N) base calls. This is '
        'applied after trimming sequences based on `min_length_fraction`.'
    ),
    'num_processes': (
        'The number of processes to use. A higher number will improve '
        'response time if the hardware resources are available. Request no '
        'more than one process per available core.'
    )
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
                 'the presence of ambiguous base calls.'),
    examples={
        'q_score': ex.q_score_example,
    },
)

importlib.import_module('q2_quality_filter._transformer')
