__version__ = '0.1'
__all__ = []

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
