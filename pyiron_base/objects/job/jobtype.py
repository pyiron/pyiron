# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron_base.objects.job.jobtypebase import JobTypeChoiceBase, JobTypeBase

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


# Dictionary of the available job types and the corresponding classes
JOB_CLASS_DICT = {'GenericMaster': 'pyiron_base.objects.job.master',
                  'ListMaster': 'pyiron_base.objects.job.list',
                  'ParallelMaster': 'pyiron_base.objects.job.parallel',
                  'ScriptJob': 'pyiron_base.objects.job.script',
                  'SerialMasterBase': 'pyiron_base.objects.job.serial'}


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
