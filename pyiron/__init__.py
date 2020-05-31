__version__ = "0.1"
__all__ = []

from pyiron.project import Project
from pyiron.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, Atoms
from pyiron.atomistics.structure.generator import create_surface, create_ase_bulk, create_structure
from pyiron.base.job.script import Notebook
from pyiron.base.settings.install import install_dialog
from pyiron.base.generic.jedi import fix_ipython_autocomplete

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# jedi fix
fix_ipython_autocomplete()

# populate namespace with default job classes, this also populates the job
# class dict in pyiron.base.job.jobtype
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.testing.randomatomistic import AtomisticExampleJob
from pyiron.dft.master.convergence_encut_parallel import ConvEncutParallel
from pyiron.dft.master.convergence_encut_serial import ConvEncutSerial
from pyiron.atomistics.master.convergence_volume import ConvergenceVolume
from pyiron.dft.master.convergence_kpoint_parallel import ConvKpointParallel
from pyiron.testing.randomatomistic import ExampleJob
from pyiron.base.master.flexible import FlexibleMaster
from pyiron.gaussian.gaussian import Gaussian
from pyiron.gpaw.gpaw import GpawJob
from pyiron.thermodynamics.hessian import HessianJob
from pyiron.lammps.lammps import Lammps
from pyiron.atomistics.master.parallel import MapMaster
from pyiron.atomistics.master.murnaghan import Murnaghan
from pyiron.dft.master.murnaghan_dft import MurnaghanDFT
from pyiron.atomistics.master.phonopy import PhonopyJob
from pyiron.quickff.quickff import QuickFF
from pyiron.interactive.scipy_minimizer import ScipyMinimizer
from pyiron.base.job.script import ScriptJob
from pyiron.atomistics.master.serial import SerialMaster
from pyiron.base.master.serial import SerialMasterBase
from pyiron.sphinx.sphinx import Sphinx
from pyiron.atomistics.job.structurecontainer import StructureContainer
from pyiron.atomistics.master.structure import StructureListMaster
from pyiron.thermodynamics.sxphonons import SxDynMat
from pyiron.interactive.sxextoptint import SxExtOptInteractive
from pyiron.thermodynamics.sxphonons import SxHarmPotTst
from pyiron.thermodynamics.sxphonons import SxPhonons
from pyiron.thermodynamics.sxphonons import SxUniqDispl
from pyiron.table.datamining import TableJob
from pyiron.vasp.vasp import Vasp
from pyiron.vasp.metadyn import VaspMetadyn
from pyiron.vasp.vaspsol import VaspSol
from pyiron.yaff.yaff import Yaff

def install():
    install_dialog()
