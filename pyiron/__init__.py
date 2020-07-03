__version__ = "0.1"
__all__ = []

from pyiron.project import Project
from pyiron.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, Atoms
from pyiron.atomistics.structure.generator import create_surface, create_ase_bulk, create_structure
from pyiron.base.job.external import Notebook
from pyiron.base.settings.install import install_dialog
from pyiron.base.generic.jedi import fix_ipython_autocomplete

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# jedi fix
fix_ipython_autocomplete()


def install():
    install_dialog()
