# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from ._filter import q_score
from ._viz_stats import visualize_stats

__version__ = pkg_resources.get_distribution('q2-quality-filter').version

__all__ = ['q_score', 'visualize_stats']
