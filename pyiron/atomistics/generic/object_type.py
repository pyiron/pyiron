# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import importlib
from pyiron.base.settings.generic import Settings

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


s = Settings()

OBJECT_CLASS_DICT = {'ThermoBulk': 'pyiron.atomistics.thermodynamics.thermo_bulk'}


class ObjectTypeChoice(object):
    def __init__(self):
        for item in list(OBJECT_CLASS_DICT.keys()):
            self.__setattr__(item, item)

    @staticmethod
    def __dir__():
        return list(OBJECT_CLASS_DICT.keys())


class ObjectType(object):

    def __new__(cls, class_name, project=None, job_name=None):
        """
        The ObjectType class is used to create light weight pyiron object that do not represent jobs (i.e., they do
        not have a status, etc.)
        
        Args:
            cls: The JobType class which contains a list of all available JOB_TYPES these are the objects which 
                 can be restored from an HDF5 File.
            class_name: The specific class name of the class this object belongs to.
            project: Project object (defines path where job will be created and stored)
            job_name (str): name of the job (must be unique within this project path needed only if object is stored in 
                            db)

        Returns:

        """
        type_lst = class_name.split(".")
        if len(type_lst) > 1:
            class_name = class_name.split()[-1][1:-2]
            object_type = class_name.split(".")[-1]
        else:
            object_type = type_lst[-1]
        for class_name in list(OBJECT_CLASS_DICT.keys()):
            if object_type == class_name:
                object_module = importlib.import_module(OBJECT_CLASS_DICT[class_name])
                object_class = getattr(object_module, class_name)
                return object_class(project, job_name)

        raise ValueError("Unknown object type: ", class_name, [job.__name__ for job in list(OBJECT_CLASS_DICT.keys())])
