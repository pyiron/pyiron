__version__ = "0.1"
__all__ = []

from pyiron.project import Project
from pyiron.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, Atoms
from pyiron.atomistics.structure.generator import create_surface, create_ase_bulk, create_structure
from pyiron_base import Notebook, install_dialog, JOB_CLASS_DICT


# Make classes available for new pyiron version
JOB_CLASS_DICT["ART"] = "pyiron.interactive.activation_relaxation_technique"
JOB_CLASS_DICT["Atoms"] = "pyiron.atomistics.structure.atoms"
JOB_CLASS_DICT["AtomisticExampleJob"] = "pyiron.testing.randomatomistic"
JOB_CLASS_DICT["ConvEncutParallel"] = "pyiron.dft.master.convergence_encut_parallel"
JOB_CLASS_DICT["ConvEncutSerial"] = "pyiron.dft.master.convergence_encut_serial"
JOB_CLASS_DICT["ConvergenceVolume"] = "pyiron.atomistics.master.convergence_volume"
JOB_CLASS_DICT["ConvKpointParallel"] = "pyiron.dft.master.convergence_kpoint_parallel"
JOB_CLASS_DICT["ElasticTensor"] = "pyiron.atomistics.master.elastic"
JOB_CLASS_DICT["ExampleJob"] = "pyiron.testing.randomatomistic"
JOB_CLASS_DICT["Gaussian"] = "pyiron.gaussian.gaussian"
JOB_CLASS_DICT["Gpaw"] = "pyiron.gpaw.gpaw"
JOB_CLASS_DICT["HessianJob"] = "pyiron.thermodynamics.hessian"
JOB_CLASS_DICT["Lammps"] = "pyiron.lammps.lammps"
JOB_CLASS_DICT["MapMaster"] = "pyiron.atomistics.master.parallel"
JOB_CLASS_DICT["Murnaghan"] = "pyiron.atomistics.master.murnaghan"
JOB_CLASS_DICT["MurnaghanDFT"] = "pyiron.dft.master.murnaghan_dft"
JOB_CLASS_DICT["PhonopyJob"] = "pyiron.atomistics.master.phonopy"
JOB_CLASS_DICT["QuasiHarmonicJob"] = "pyiron.atomistics.master.quasi"
JOB_CLASS_DICT["QuickFF"] = "pyiron.quickff.quickff"
JOB_CLASS_DICT["ScipyMinimizer"] = "pyiron.interactive.scipy_minimizer"
JOB_CLASS_DICT["SerialMaster"] = "pyiron.atomistics.master.serial"
JOB_CLASS_DICT["Sphinx"] = "pyiron.sphinx.sphinx"
JOB_CLASS_DICT["StructureContainer"] = "pyiron.atomistics.job.structurecontainer"
JOB_CLASS_DICT["StructureListMaster"] = "pyiron.atomistics.master.structure"
JOB_CLASS_DICT["SQSJob"] = "pyiron.atomistics.job.sqs"
JOB_CLASS_DICT["SQSMaster"] = "pyiron.atomistics.master.sqsmaster"
JOB_CLASS_DICT["SxDynMat"] = "pyiron.thermodynamics.sxphonons"
JOB_CLASS_DICT["SxExtOptInteractive"] = "pyiron.interactive.sxextoptint"
JOB_CLASS_DICT["SxHarmPotTst"] = "pyiron.thermodynamics.sxphonons"
JOB_CLASS_DICT["SxPhonons"] = "pyiron.thermodynamics.sxphonons"
JOB_CLASS_DICT["SxUniqDispl"] = "pyiron.thermodynamics.sxphonons"
JOB_CLASS_DICT["TableJob"] = "pyiron.table.datamining"
JOB_CLASS_DICT["Vasp"] = "pyiron.vasp.vasp"
JOB_CLASS_DICT["VaspMetadyn"] = "pyiron.vasp.metadyn"
JOB_CLASS_DICT["VaspSol"] = "pyiron.vasp.vaspsol"
JOB_CLASS_DICT["Yaff"] = "pyiron.yaff.yaff"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def install():
    install_dialog()
