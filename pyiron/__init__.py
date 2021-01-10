__version__ = "0.1"
__all__ = []

from pyiron_atomistics.project import Project
from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, Atoms
from pyiron_base import Notebook, install_dialog, JOB_CLASS_DICT

# To maintain backwards compatibility until we deprecate the old structure creation functions:
from pyiron_atomistics.atomistics.structure.factory import StructureFactory as _StructureFactory
create_surface = _StructureFactory.surface
create_ase_bulk = _StructureFactory.ase_bulk
create_structure = _StructureFactory.crystal

# Fix modules for backwards compatibility
from pyiron_atomistics import \
    atomistics, dft, gaussian, gpaw, interactive, lammps, quickff, \
    sphinx, table, testing, thermodynamics, vasp, yaff, project

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def install():
    install_dialog()
