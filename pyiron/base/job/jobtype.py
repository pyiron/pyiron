# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import importlib
import inspect
import pkgutil
from six import with_metaclass
from pyiron.base.settings.generic import Settings

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


JOB_CLASS_DICT = {'GenericMaster': 'pyiron.base.master.generic',
                  'ListMaster': 'pyiron.base.master.list',
                  'ScriptJob': 'pyiron.base.job.script',
                  'SerialMasterBase': 'pyiron.base.master.serial',
                  'SerialMaster': 'pyiron.atomistics.master.serial',
                  'ParallelMaster': 'pyiron.atomistics.master.parallel',
                  'Murnaghan': 'pyiron.atomistics.master.murnaghan',
                  'MurnaghanInt': 'pyiron.atomistics.master.murninterface',
                  'PhonopyJob': 'pyiron.atomistics.master.phonopy',
                  'PhonopyMaster': 'pyiron.atomistics.master.phonopy',
                  'PhonopyMaster2': 'pyiron.atomistics.master.phonopy',
                  'ConvergenceVolume': 'pyiron.atomistics.master.convergence_volume',
                  'StructureContainer': 'pyiron.atomistics.job.structurecontainer',
                  'ConvergenceEncutParallel': 'pyiron.dft.master.convergence_encut_parallel',
                  'ConvergenceEncutSerial': 'pyiron.dft.master.convergence_encut_serial',
                  'ConvergenceKpointParallel': 'pyiron.dft.master.convergence_kpoint_parallel',
                  'MurnaghanDFT': 'pyiron.dft.master.murnaghan_dft',
                  'Lammps': 'pyiron.lammps.lammps',
                  'LammpsInt': 'pyiron.lammps.lammpsinterface',
                  'LammpsInt2': 'pyiron.lammps.lammpsinterface',
                  'AtomisticExampleJob': 'pyiron.testing.randomatomistic',
                  'ExampleJob': 'pyiron.testing.randomatomistic',
                  'Vasp': 'pyiron.vasp.vasp',
                  'VaspInt': 'pyiron.vasp.vaspinterface',
                  'VaspInt2': 'pyiron.vasp.vaspinterface',
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
        self._job_class_dict = None
        self.job_class_dict = self._extend_job_dict(JOB_CLASS_DICT)

    @property
    def job_class_dict(self):
        return self._job_class_dict

    @job_class_dict.setter
    def job_class_dict(self, job_class_dict):
        self._job_class_dict = job_class_dict
        for item in list(self._job_class_dict.keys()):
            self.__setattr__(item, item)

    @staticmethod
    def _extend_job_dict(job_dict):
        def derived_from_generic_job(obj):
            return 'pyiron.base.job.generic.GenericJob' in [subcls.__module__ + '.' + subcls.__name__
                                                            for subcls in obj.__mro__]

        for d in [{name: obj.__module__
                   for name, obj in inspect.getmembers(importlib.import_module(name))
                   if inspect.isclass(obj) and derived_from_generic_job(obj)}
                  for finder, name, ispkg in pkgutil.iter_modules()
                  if name.startswith('pyiron_')]:
            job_dict.update(d)
        return job_dict

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
        if isinstance(class_name, str):
            job_class = cls.convert_str_to_class(job_class_dict=cls.job_class_dict, class_name=class_name)
        elif inspect.isclass(class_name):
            job_class = class_name
        else:
            raise TypeError()
        job = job_class(project, job_name)
        if job.status.aborted:
            job.logger.warn('Job aborted - please remove it and run again! {}'.format(job.job_name))
        if not job.status.initialized:
            Settings().logger.info("job_exists -> load from hdf")
            job.from_hdf()
        return job

    @staticmethod
    def convert_str_to_class(job_class_dict, class_name):
        job_type_lst = class_name.split(".")
        if len(job_type_lst) > 1:
            class_name = class_name.split()[-1][1:-2]
            job_type = class_name.split(".")[-1]
        else:
            job_type = job_type_lst[-1]
        for job_class_name in list(job_class_dict.keys()):  # for job_class in cls.JOB_CLASSES:
            if job_type == job_class_name:
                job_module = importlib.import_module(job_class_dict[job_class_name])
                job_class = getattr(job_module, job_class_name)
                return job_class
        raise ValueError("Unknown job type: ", class_name, [job for job in list(job_class_dict.keys())])
