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

_pyiron_modules = [
    importlib.import_module(name) 
    for finder, name, ispkg in pkgutil.iter_modules() 
    if name.startswith('pyiron_') and name not in ['pyiron_base', 'pyiron_atomistics']
]

__pyiron_versions__ = {'pyiron': __version__}
for pyiron_module in _pyiron_modules:
    name = pyiron_module.__name__
    try:
        version = pyiron_module.__version__
    except AttributeError:
        version = "unknown"
    __pyiron_versions__[name] = version

del _pyiron_modules

def get_pyiron_versions():
    return __pyiron_versions__


def install():
    install_dialog()
