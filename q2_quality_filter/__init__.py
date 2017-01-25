import pkg_resources

from ._filter import q_score
from ._viz_stats import visualize_stats

__version__ = pkg_resources.get_distribution('q2-quality-filter').version

__all__ = ['q_score', 'visualize_stats']
