# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
Template class to list the required properties and functions for every pyiron object.
"""

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class PyironObject(object):
    """
    Template class to list the required properties and functions for every pyiron object.
    """

    @property
    def id(self):
        """
        Every pyiron object should have the ability to be stored in the database

        Returns:
            int: object id
        """
        raise NotImplementedError("id should be implemented in the derived class")

    @id.setter
    def id(self, new_id):
        """
        Every pyiron object should have the ability to be stored in the database

        Args:
            new_id (int): object id
        """
        raise NotImplementedError("id should be implemented in the derived class")

    @property
    def master_id(self):
        """
        If the pyiron object belongs to a series of objects, series object is linked by the master id

        Returns:
            int: master id
        """
        raise NotImplementedError(
            "master_id should be implemented in the derived class"
        )

    @master_id.setter
    def master_id(self, master_id):
        """
        If the pyiron object belongs to a series of objects, series object is linked by the master id

        Args:
            master_id (int): master id
        """
        raise NotImplementedError(
            "master_id should be implemented in the derived class"
        )

    @property
    def parent_id(self):
        """
        If the pyiron object belongs to a serial series of objects, the predecessor is linked by the parent id

        Returns:
            int: parent id
        """
        raise NotImplementedError(
            "master_id should be implemented in the derived class"
        )

    @parent_id.setter
    def parent_id(self, parent_id):
        """
        If the pyiron object belongs to a serial series of objects, the predecessor is linked by the parent id

        Args:
            parent_id (int): parent id
        """
        raise NotImplementedError(
            "master_id should be implemented in the derived class"
        )

    @property
    def child_ids(self):
        """
        If the pyiron object is a meta object which includes a series of objects these objects ids are listed as child ids.

        Returns:
            list: list of child ids
        """
        raise NotImplementedError(
            "child_ids should be implemented in the derived class"
        )

    def save(self):
        """
        Store the pyiron object in the HDF5 file an create a corresponding database entry.
        """
        raise NotImplementedError("save() should be implemented in the derived class")

    def remove(self):
        """
        Remove the pyiron obect from the database and delete the HDF5 file
        """
        raise NotImplementedError("remove() should be implemented in the derived class")

    def load(self, job_specifier, convert_to_object=True):
        """
        Load a pyiron object from the database

        Args:
            job_specifier (str, int): identifier of the pyiron object - this needs to be unique in the project.
            convert_to_object (bool): [True/False] - it is faster to only load the HDF5 access but the pyiron object
                                      offers a whole more functionality.

        Returns:
            PyironObject: the pyiron object
        """
        raise NotImplementedError("load() should be implemented in the derived class")

    def inspect(self, job_specifier):
        """
        Inspect a pyiron object from the database - the same like load(self, job_specifier, convert_to_object=False).
        The inspect mode provides quick access to the HDF5 file without loading the full object, which is especially
        helpful if only a specific section of the output is required.

        Args:
            job_specifier (str, int): identifier of the pyiron object - this needs to be unique in the project.

        Returns:
            PyironObject: a reduces pyiron object like JobPath
        """
        raise NotImplementedError("insect() should be implemented in the derived class")

    def copy(self):
        """
        Copy the pyiron object - this copies only the link to the HDF5 file - not the content of the HDF5 file.

        Returns:
            PyironObject: a copy of the pyiron object
        """
        raise NotImplementedError("copy() should be implemented in the derived class")

    def copy_to(self, new_project=None):
        """
        Copy the pyiron object including the HDF5 file to an new location

        Args:
            new_project (ProjectHDFio): new location

        Returns:
            PyironObject: a copy of the pyiron object
        """
        raise NotImplementedError("copy() should be implemented in the derived class")

    def move_to(self, new_project):
        """
        Move the pyiron object including the HDF5 file to an new location

        Args:
            new_project (ProjectHDFio): new location
        """
        raise NotImplementedError("move() should be implemented in the derived class")

    def to_hdf(self, hdf, group_name="group"):
        """
        Store the PyironObject in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object
            group_name (str): HDF5 subgroup name - optional
        """
        raise NotImplementedError("to_hdf() should be implemented in the derived class")

    def from_hdf(self, hdf, group_name="group"):
        """
        Restore the PyironObject from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object
            group_name (str): HDF5 subgroup name - optional
        """
        raise NotImplementedError(
            "from_hdf() should be implemented in the derived class"
        )
