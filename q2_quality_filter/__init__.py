import pkg_resources

from ._filter import basic


__version__ = pkg_resources.get_distribution('q2-quality-filter').version

__all__ = ['basic']
