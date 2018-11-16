# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import importlib
import inspect
import pkgutil
from six import with_metaclass
from pyiron.base.core.settings.generic import Settings

"""
Jobtype class to create GenericJob type objects 
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
                  'SerialMasterBase': 'pyiron.base.objects.job.serial',
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

    Args:
        job_class_dict: dictionary with the jobtypes to choose from.
    """
    def __init__(self):
        for d in [{name: obj.__module__
                   for name, obj in inspect.getmembers(importlib.import_module(name))
                   if inspect.isclass(obj) and 'GenericJob' in [subcls.__name__ for subcls in obj.__mro__]}
                  for finder, name, ispkg in pkgutil.iter_modules()
                  if name.startswith('pyiron_')]:
            JOB_CLASS_DICT.update(d)
        self.job_class_dict = JOB_CLASS_DICT
        for item in list(self.job_class_dict.keys()):
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
    def __new__(cls, class_name, project, job_name, job_class_dict):
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
        job_type_lst = class_name.split(".")
        if len(job_type_lst) > 1:
            class_name = class_name.split()[-1][1:-2]
            job_type = class_name.split(".")[-1]
        else:
            job_type = job_type_lst[-1]
        for job_class_name in list(cls.job_class_dict.keys()):  # for job_class in cls.JOB_CLASSES:
            if job_type == job_class_name:
                job_module = importlib.import_module(cls.job_class_dict[job_class_name])
                job_class = getattr(job_module, job_class_name)
                job = job_class(project, job_name)
                if job.status.aborted:
                    job.logger.warn('Job aborted - please remove it and run again! {}'.format(job.job_name))
                if not job.status.initialized:
                    Settings().logger.info("job_exists -> load from hdf")
                    job.from_hdf()
                return job
        raise ValueError("Unknown job type: ", class_name, [job.__name__ for job in list(cls.job_class_dict.keys())])
