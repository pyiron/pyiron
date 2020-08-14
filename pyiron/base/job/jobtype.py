# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import importlib
import inspect
import pkgutil
from six import with_metaclass
from pyiron.base.generic.util import static_isinstance

"""
Jobtype class to create GenericJob type objects
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


JOB_CLASS_DICT = {
    "Atoms": "pyiron.atomistics.structure.atoms",
    "AtomisticExampleJob": "pyiron.testing.randomatomistic",
    "ConvEncutParallel": "pyiron.dft.master.convergence_encut_parallel",
    "ConvEncutSerial": "pyiron.dft.master.convergence_encut_serial",
    "ConvergenceVolume": "pyiron.atomistics.master.convergence_volume",
    "ConvKpointParallel": "pyiron.dft.master.convergence_kpoint_parallel",
    "ExampleJob": "pyiron.testing.randomatomistic",
    "FlexibleMaster": "pyiron.base.master.flexible",
    "Gaussian": "pyiron.gaussian.gaussian",
    "GpawJob": "pyiron.gpaw.gpaw",
    "HessianJob": "pyiron.thermodynamics.hessian",
    "Lammps": "pyiron.lammps.lammps",
    "MapMaster": "pyiron.atomistics.master.parallel",
    "Murnaghan": "pyiron.atomistics.master.murnaghan",
    "MurnaghanDFT": "pyiron.dft.master.murnaghan_dft",
    "PhonopyJob": "pyiron.atomistics.master.phonopy",
    "QuickFF": "pyiron.quickff.quickff",
    "ScipyMinimizer": "pyiron.interactive.scipy_minimizer",
    "ScriptJob": "pyiron.base.job.script",
    "SerialMaster": "pyiron.atomistics.master.serial",
    "SerialMasterBase": "pyiron.base.master.serial",
    "Sphinx": "pyiron.sphinx.sphinx",
    "StructureContainer": "pyiron.atomistics.job.structurecontainer",
    "StructureListMaster": "pyiron.atomistics.master.structure",
    "SQSJob": "pyiron.atomistics.job.sqs",
    "SxDynMat": "pyiron.thermodynamics.sxphonons",
    "SxExtOptInteractive": "pyiron.interactive.sxextoptint",
    "SxHarmPotTst": "pyiron.thermodynamics.sxphonons",
    "SxPhonons": "pyiron.thermodynamics.sxphonons",
    "SxUniqDispl": "pyiron.thermodynamics.sxphonons",
    "TableJob": "pyiron.table.datamining",
    "Vasp": "pyiron.vasp.vasp",
    "VaspMetadyn": "pyiron.vasp.metadyn",
    "VaspSol": "pyiron.vasp.vaspsol",
    "Yaff": "pyiron.yaff.yaff",
}


class Singleton(type):
    """
    Implemented with suggestions from

    http://stackoverflow.com/questions/6760685/creating-a-singleton-in-python

    """

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class JobTypeChoice(with_metaclass(Singleton)):
    """
    Helper class to choose the job type directly from the project, autocompletion is enabled by overwriting the
    __dir__() function.
    """

    def __init__(self):
        self._job_class_dict = None
        self.job_class_dict = JOB_CLASS_DICT

    @property
    def job_class_dict(self):
        return self._job_class_dict

    @job_class_dict.setter
    def job_class_dict(self, job_class_dict):
        self._job_class_dict = job_class_dict
        for item in list(self._job_class_dict.keys()):
            self.__setattr__(item, item)

    def __dir__(self):
        """
        Enable autocompletion by overwriting the __dir__() function.
        """
        return list(self.job_class_dict.keys())


class JobType(object):
    """
    The JobTypeBase class creates a new object of a given class type.
    """

    def __new__(cls, class_name, project, job_name, job_class_dict, delete_existing_job=False):
        """
        The __new__() method allows to create objects from other classes - the class selected by class_name

        Args:
            class_name (str): The specific class name of the class this object belongs to.
            project (Project): Project object (defines path where job will be created and stored)
            job_name (str): name of the job (must be unique within this project path)
            job_class_dict (dict): dictionary with the jobtypes to choose from.

        Returns:
            GenericJob: object of type class_name
        """
        cls.job_class_dict = job_class_dict
        if isinstance(class_name, str):
            job_class = cls.convert_str_to_class(
                job_class_dict=cls.job_class_dict, class_name=class_name
            )
        elif inspect.isclass(class_name):
            job_class = class_name
        else:
            raise TypeError()
        job = job_class(project, job_name)
        if delete_existing_job:
            job.remove()
            job = job_class(project, job_name)
        if job.status.aborted:
            job.logger.warning(
                "Job aborted - please remove it and run again! {}".format(job.job_name)
            )
        if not job.status.initialized:
            job.from_hdf()
        if job.status.finished or job.status.collect:
            job.set_input_to_read_only()
        return job

    @staticmethod
    def convert_str_to_class(job_class_dict, class_name):
        """
        convert the name of a class to the corresponding class object - only for pyiron internal classes.

        Args:
            job_class_dict (dict):
            class_name (str):

        Returns:
            (class):
        """
        job_type_lst = class_name.split(".")
        if len(job_type_lst) > 1:
            class_name = class_name.split()[-1][1:-2]
            job_type = class_name.split(".")[-1]
        else:
            job_type = job_type_lst[-1]
        for job_class_name in list(
            job_class_dict.keys()
        ):  # for job_class in cls.JOB_CLASSES:
            if job_type == job_class_name:
                job_module = importlib.import_module(job_class_dict[job_class_name])
                job_class = getattr(job_module, job_class_name)
                return job_class
        raise ValueError(
            "Unknown job type: ",
            class_name,
            [job for job in list(job_class_dict.keys())],
        )
