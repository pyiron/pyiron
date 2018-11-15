from pyiron.vasp.vasp import Vasp
__all__ = ['Vasp']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
