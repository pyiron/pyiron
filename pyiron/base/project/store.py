# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from datetime import datetime
from pyiron.base.job.generic import GenericJob

"""
Class for storing user aggregated information in an pyiron object 
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "testing"
__date__ = "Sep 1, 2017"


class ProjectStore(GenericJob):
    """
    The ProjectStore object, is derived from the GenericJob class and allows the user to store
    aggregated information in an HDF5 file associated with the corresponding project. To the user
    the ProjectStore object behaves like a dictionary.

    Args:
        project: Project object (defines path where job will be created and stored)
        job_name: name of the job (must be unique within this project path)

    Attributes:

        .. attribute:: key

            keys of the ProjectStore object (like a dictionary)

        .. attribute:: items

            items of the ProjectStore object (like a dictionary)

        .. attribute:: values

            values of the ProjectStore object (like a dictionary)

        .. attribute:: time_created

            date when the ProjecStore object was created

        .. attribute:: time_modified

            date when the ProjecStore object was modified
    """
    def __init__(self, project, job_name):
        super(ProjectStore, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "ProjectStore"
        self._lib = {'available': True, 'enabled': True}
        self._store = {}

    @property
    def keys(self):
        """
        a set-like object providing a view on ProjectStore's dictionary keys
        """
        return self._store.keys()

    @property
    def items(self):
        """
        a set-like object providing a view on ProjectStore's dictionary items
        """
        return self._store.items()

    @property
    def values(self):
        """
        an object providing a view on ProjectStore's dictionary values
        """
        return self._store.values()

    @property
    def time_created(self):
        """
        Return the date when the ProjectStore object was created

        Returns:
            DateTime: the date when the ProjectStore object was created
        """
        if self.job_id:
            return self.project.db.get_item_by_id(self._job_id)["timestart"]
        return None

    @property
    def time_modified(self):
        """
        Return the date when the ProjectStore object was modified

        Returns:
            DateTime: the date when the ProjectStore object was modified
        """
        if self.job_id:
            return self.project.db.get_item_by_id(self._job_id)["timestop"]
        return None

    def run_if_lib(self):
        """
        Internal function to handle jobs with Python based executables

        Returns:
            int: Database ID of the ProjectStore object
        """
        job_id = self.save()
        self._write_to_database()
        return job_id

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore object from hdf5 format

        Args:
            hdf: Optional hdf5 file, otherwise self is used.
            group_name (str): Optional hdf5 group in the hdf5 file.

        """
        super(ProjectStore, self).from_hdf(hdf=hdf, group_name=group_name)
        for node in self.project_hdf5.list_nodes():
            value = self.project_hdf5[node]
            self._store[node] = value
            self.__setattr__(node, value)

    def __setitem__(self, key, value):
        """
        Store values in the ProjectStore

        Args:
            key (str): Key for the dictionary
            value: corresponding values

        """
        self._hdf5[key] = value
        self._store[key] = value
        self.__setattr__(key, value)
        self.run()

    def _run_if_finished(self, run_again=False):
        """
        Internal function overwriting the default behaviour when the job is finished,
        to update the database entry when the job was modified

        Args:
            run_again (bool): not used for this job type

        """
        self._write_to_database()

    def _write_to_database(self):
        """
        If a job_id exists update the timestop entry in the database, to validate when this object was updated.
        """
        if self.job_id:
            self.project.db.item_update({'timestop': datetime.now()}, self._job_id)

    def write_input(self):
        """
        Implement required function template - even though it is not required for this job type.
        """
        pass

    def collect_output(self):
        """
        Implement required function template - even though it is not required for this job type.
        """
        pass

    def collect_logfiles(self):
        """
        Implement required function template - even though it is not required for this job type.
        """
        pass
