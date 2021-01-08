__version__ = "0.1"
__all__ = []

from pyiron_atomistic.project import Project
from pyiron_atomistic.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, Atoms
from pyiron_base import Notebook, install_dialog, JOB_CLASS_DICT

# To maintain backwards compatibility until we deprecate the old structure creation functions:
from pyiron_atomistic.atomistics.structure.factory import StructureFactory as _StructureFactory
create_surface = _StructureFactory.surface
create_ase_bulk = _StructureFactory.ase_bulk
create_structure = _StructureFactory.crystal

# Make classes available for new pyiron version
JOB_CLASS_DICT["ART"] = "pyiron_atomistic.interactive.activation_relaxation_technique"
JOB_CLASS_DICT["AtomisticExampleJob"] = "pyiron_atomistic.testing.randomatomistic"
JOB_CLASS_DICT["ConvEncutParallel"] = "pyiron_atomistic.dft.master.convergence_encut_parallel"
JOB_CLASS_DICT["ConvEncutSerial"] = "pyiron_atomistic.dft.master.convergence_encut_serial"
JOB_CLASS_DICT["ConvergenceVolume"] = "pyiron_atomistic.atomistics.master.convergence_volume"
JOB_CLASS_DICT["ConvKpointParallel"] = "pyiron_atomistic.dft.master.convergence_kpoint_parallel"
JOB_CLASS_DICT["ElasticTensor"] = "pyiron_atomistic.atomistics.master.elastic"
JOB_CLASS_DICT["ExampleJob"] = "pyiron_atomistic.testing.randomatomistic"
JOB_CLASS_DICT["Gaussian"] = "pyiron_atomistic.gaussian.gaussian"
JOB_CLASS_DICT["Gpaw"] = "pyiron_atomistic.gpaw.gpaw"
JOB_CLASS_DICT["HessianJob"] = "pyiron_atomistic.thermodynamics.hessian"
JOB_CLASS_DICT["Lammps"] = "pyiron_atomistic.lammps.lammps"
JOB_CLASS_DICT["MapMaster"] = "pyiron_atomistic.atomistics.master.parallel"
JOB_CLASS_DICT["Murnaghan"] = "pyiron_atomistic.atomistics.master.murnaghan"
JOB_CLASS_DICT["MurnaghanDFT"] = "pyiron_atomistic.dft.master.murnaghan_dft"
JOB_CLASS_DICT["PhonopyJob"] = "pyiron_atomistic.atomistics.master.phonopy"
JOB_CLASS_DICT["QuasiHarmonicJob"] = "pyiron_atomistic.atomistics.master.quasi"
JOB_CLASS_DICT["QuickFF"] = "pyiron_atomistic.quickff.quickff"
JOB_CLASS_DICT["ScipyMinimizer"] = "pyiron_atomistic.interactive.scipy_minimizer"
JOB_CLASS_DICT["SerialMaster"] = "pyiron_atomistic.atomistics.master.serial"
JOB_CLASS_DICT["Sphinx"] = "pyiron_atomistic.sphinx.sphinx"
JOB_CLASS_DICT["StructureContainer"] = "pyiron_atomistic.atomistics.job.structurecontainer"
JOB_CLASS_DICT["StructureListMaster"] = "pyiron_atomistic.atomistics.master.structure"
JOB_CLASS_DICT["SQSJob"] = "pyiron_atomistic.atomistics.job.sqs"
JOB_CLASS_DICT["SQSMaster"] = "pyiron_atomistic.atomistics.master.sqsmaster"
JOB_CLASS_DICT["SxDynMat"] = "pyiron_atomistic.thermodynamics.sxphonons"
JOB_CLASS_DICT["SxExtOptInteractive"] = "pyiron_atomistic.interactive.sxextoptint"
JOB_CLASS_DICT["SxHarmPotTst"] = "pyiron_atomistic.thermodynamics.sxphonons"
JOB_CLASS_DICT["SxPhonons"] = "pyiron_atomistic.thermodynamics.sxphonons"
JOB_CLASS_DICT["SxUniqDispl"] = "pyiron_atomistic.thermodynamics.sxphonons"
JOB_CLASS_DICT["TableJob"] = "pyiron_atomistic.table.datamining"
JOB_CLASS_DICT["Vasp"] = "pyiron_atomistic.vasp.vasp"
JOB_CLASS_DICT["VaspMetadyn"] = "pyiron_atomistic.vasp.metadyn"
JOB_CLASS_DICT["VaspSol"] = "pyiron_atomistic.vasp.vaspsol"
JOB_CLASS_DICT["Yaff"] = "pyiron_atomistic.yaff.yaff"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def install():
    install_dialog()
