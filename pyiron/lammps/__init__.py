from lammps.lammps import Lammps
__all__ = ['Lammps']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
