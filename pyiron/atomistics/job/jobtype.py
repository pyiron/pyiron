# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import importlib
import inspect
import pkgutil
from pyiron.base.objects.job.generic import GenericJob
from pyiron.base.objects.job.jobtypebase import JobTypeChoiceBase, JobTypeBase

"""
Jobtype class to create GenericJob type objects - the choice of job types is limited to the JOB_CLASS_DICT
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

JOB_CLASS_DICT = {'GenericMaster': 'pyiron.base.objects.job.master',
                  'ListMaster': 'pyiron.base.objects.job.list',
                  'ScriptJob': 'pyiron.base.objects.job.script',
                  'SerialMaster': 'pyiron.atomistics.job.serial',
                  'ParallelMaster': 'pyiron.atomistics.job.parallel',
                  'Murnaghan': 'pyiron.atomistics.hamilton.murnaghan',
                  'MurnaghanInt': 'pyiron.atomistics.hamilton.murninterface',
                  'PhonopyJob': 'pyiron.atomistics.hamilton.phonopy',
                  'ConvergenceVolume': 'pyiron.atomistics.hamilton.convergence_volume',
                  'StructureContainer': 'pyiron.atomistics.structure.structurecontainer',
                  'ConvergenceEncutParallel': 'pyiron.dft.convergence_encut_parallel',
                  'ConvergenceEncutSerial': 'pyiron.dft.convergence_encut_serial',
                  'ConvergenceKpointParallel': 'pyiron.dft.convergence_kpoint_parallel',
                  'MurnaghanDFT': 'pyiron.dft.murnaghan_dft',
                  'Lammps': 'pyiron.lammps.lammps',
                  'LammpsInt': 'pyiron.lammps.lammpsinterface',
                  'AtomisticExampleJob': 'pyiron.testing.randomatomistic',
                  'ExampleJob': 'pyiron.testing.randomatomistic',
                  'Vasp': 'pyiron.vasp.vasp',
                  'VaspInt': 'pyiron.vasp.vaspinterface',
}

for d in [{name: obj.__module__
           for name, obj in inspect.getmembers(importlib.import_module(name))
           if inspect.isclass(obj) and issubclass(obj, GenericJob)}
          for finder, name, ispkg in pkgutil.iter_modules()
          if name.startswith('pyiron_')]:
    JOB_CLASS_DICT.update(d)


class JobTypeChoice(JobTypeChoiceBase):
    """
    Helper class to choose the job type directly from the project, autocompletion is enabled by overwriting the
    __dir__() function.
    """
    def __init__(self):
        super(JobTypeChoice, self).__init__(job_class_dict=JOB_CLASS_DICT)


class JobType(JobTypeBase):
    """
    The JobType class creates a new object of a given class type.
    """
    def __new__(cls, class_name, project, job_name):
        """
        The __new__() method allows to create objects from other classes - the class selected by class_name

        Args:
            class_name (str): The specific class name of the class this object belongs to.
            project (Project): Project object (defines path where job will be created and stored)
            job_name (str): name of the job (must be unique within this project path)

        Returns:
            GenericJob: object of type class_name
        """
        return super(JobType, cls).__new__(cls=cls, class_name=class_name, project=project, job_name=job_name,
                                           job_class_dict=JOB_CLASS_DICT)
