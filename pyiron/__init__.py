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
import sys
import pkgutil
import importlib
from pyiron_atomistics import \
    atomistics, dft, gaussian, gpaw, interactive, lammps, quickff, \
    sphinx, table, testing, thermodynamics, vasp, yaff, project
sys.modules["pyiron.atomistics"] = atomistics
sys.modules["pyiron.dft"] = dft
sys.modules["pyiron.gaussian"] = gaussian
sys.modules["pyiron.gpaw"] = gpaw
sys.modules["pyiron.interactive"] = interactive
sys.modules["pyiron.lammps"] = lammps
sys.modules["pyiron.quickff"] = quickff
sys.modules["pyiron.sphinx"] = sphinx
sys.modules["pyiron.table"] = table
sys.modules["pyiron.testing"] = testing
sys.modules["pyiron.thermodynamics"] = thermodynamics
sys.modules["pyiron.vasp"] = vasp
sys.modules["pyiron.yaff"] = yaff
sys.modules["pyiron.project"] = project

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

_ = [
    importlib.import_module(name) 
    for finder, name, ispkg in pkgutil.iter_modules() 
    if name.startswith('pyiron_') and name not in ['pyiron_base', 'pyiron_atomistics']
]

def install():
    install_dialog()
