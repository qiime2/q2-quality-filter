import qiime2.plugin

from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import SequencesWithQuality

import q2_quality_filter
from q2_quality_filter._filter import basic


plugin = qiime2.plugin.Plugin(
    name='quality-filter',
    version=q2_quality_filter.__version__,
    website='Website for q2-quality-filter',
    package='q2_quality_filter',
    user_support_text=None,
    citation_text=None
)


plugin.methods.register_function(
    function=basic,
    inputs={'demux': SampleData[SequencesWithQuality]},
    parameters={
        'minimum_quality': qiime2.plugin.Int,
        'quality_window': qiime2.plugin.Int,
        'min_length_fraction': qiime2.plugin.Float,
        'maximum_ambiguous': qiime2.plugin.Int
    },
    outputs=[
        ('filtered_sequences', SampleData[SequencesWithQuality])
    ],
    name='Quality filterer',
    description=('This method filters sequence based on quality scores and '
                 'the presence of ambiguous base calls.')
)
