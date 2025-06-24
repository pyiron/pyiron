__version__ = "0.1"
__all__ = []

from pyiron_atomistics.project import Project
from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, Atoms
from pyiron_base import Notebook, install_dialog, JOB_CLASS_DICT

# To maintain backwards compatibility until we deprecate the old structure creation functions:
from pyiron_atomistics.atomistics.structure.factory import StructureFactory as _StructureFactory
create_surface = _StructureFactory().surface
create_ase_bulk = _StructureFactory().ase_bulk
create_structure = _StructureFactory().crystal

# Fix modules for backwards compatibility
import sys
import pkgutil
import importlib
import warnings
from pyiron_atomistics import \
    atomistics, dft, gpaw, interactive, lammps, sphinx, \
    table, testing, thermodynamics, vasp, project
sys.modules["pyiron.atomistics"] = atomistics
sys.modules["pyiron.dft"] = dft
sys.modules["pyiron.gpaw"] = gpaw
sys.modules["pyiron.interactive"] = interactive
sys.modules["pyiron.lammps"] = lammps
sys.modules["pyiron.sphinx"] = sphinx
sys.modules["pyiron.table"] = table
sys.modules["pyiron.testing"] = testing
sys.modules["pyiron.thermodynamics"] = thermodynamics
sys.modules["pyiron.vasp"] = vasp
sys.modules["pyiron.project"] = project

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

for finder, name, ispkg in pkgutil.iter_modules():
    if name.startswith('pyiron_') and name not in ['pyiron_base', 'pyiron_atomistics']:
        try:
            _ = importlib.import_module(name) 
        except Exception as err:
            warnings.warn(f"Err: Failed to import pyiron module {name} due to {err}.")

def install():
    install_dialog()
