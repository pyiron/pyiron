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
sys.modules["pyiron.atomistics"] = __import__(pyiron_atomistics.atomistics)
sys.modules["pyiron.dft"] = __import__(pyiron_atomistics.dft)
sys.modules["pyiron.gaussian"] = __import__(pyiron_atomistics.gaussian)
sys.modules["pyiron.gpaw"] = __import__(pyiron_atomistics.gpaw)
sys.modules["pyiron.interactive"] = __import__(pyiron_atomistics.interactive)
sys.modules["pyiron.lammps"] = __import__(pyiron_atomistics.lammps)
sys.modules["pyiron.quickff"] = __import__(pyiron_atomistics.quickff)
sys.modules["pyiron.sphinx"] = __import__(pyiron_atomistics.sphinx)
sys.modules["pyiron.table"] = __import__(pyiron_atomistics.table)
sys.modules["pyiron.testing"] = __import__(pyiron_atomistics.testing)
sys.modules["pyiron.thermodynamics"] = __import__(pyiron_atomistics.thermodynamics)
sys.modules["pyiron.vasp"] = __import__(pyiron_atomistics.vasp)
sys.modules["pyiron.yaff"] = __import__(pyiron_atomistics.yaff)
sys.modules["pyiron.project"] = __import__(pyiron_atomistics.project)

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def install():
    install_dialog()
