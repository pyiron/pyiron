# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import posixpath
from pyiron.base.project.path import GenericPath
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.base.job.core import JobCore
from pyiron.base.project.generic import Project

"""
The JobPath class enables quick access to the HDF5 data file without loading the full object 
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class JobPath(JobCore):
    """
    The JobPath class is derived from the JobCore and is used as a lean version of the GenericJob class. Instead of
    loading the full pyiron object the JobPath class only provides access to the HDF5 file, which should be enough
    for most analysis.

    Args:
        db (DatabaseAccess): database object
        job_id (int): Job ID - optional, but either a job ID or a database entry db_entry has to be provided.
        db_entry (dict): database entry {"job":, "subjob":, "projectpath":, "project":, "hamilton":, "hamversion":,
                                         "status":} and optional entries are {"id":, "masterid":, "parentid":}
        user (str): current unix/linux/windows user who is running pyiron

    Attributes:

        .. attribute:: job_name

            name of the job, which has to be unique within the project

        .. attribute:: status

            execution status of the job, can be one of the following [initialized, appended, created, submitted, running,
                                                                      aborted, collect, suspended, refresh, busy, finished]

        .. attribute:: job_id

            unique id to identify the job in the pyiron database

        .. attribute:: parent_id

            job id of the predecessor job - the job which was executed before the current one in the current job series

        .. attribute:: master_id

            job id of the master job - a meta job which groups a series of jobs, which are executed either in parallel or in
            serial.

        .. attribute:: child_ids

            list of child job ids - only meta jobs have child jobs - jobs which list the meta job as their master

        .. attribute:: project

            Project instance the jobs is located in

        .. attribute:: project_hdf5

            ProjectHDFio instance which points to the HDF5 file the job is stored in

        .. attribute:: job_info_str

            short string to describe the job by it is job_name and job ID - mainly used for logging

        .. attribute:: working_directory

            working directory of the job is executed in - outside the HDF5 file

        .. attribute:: path

            path to the job as a combination of absolute file system path and path within the HDF5 file.

        .. attribute:: is_root

            boolean if the HDF5 object is located at the root level of the HDF5 file

        .. attribute:: is_open

            boolean if the HDF5 file is currently opened - if an active file handler exists

        .. attribute:: is_empty

            boolean if the HDF5 file is empty

        .. attribute:: base_name

            name of the HDF5 file but without any file extension

        .. attribute:: file_path

            directory where the HDF5 file is located

        .. attribute:: h5_path

            path inside the HDF5 file - also stored as absolute path
    """
    def __init__(self, db, job_id=None, db_entry=None, user=None):
        if not db_entry:
            db_entry = db.get_item_by_id(job_id)
        if db_entry is None:
            raise ValueError("job ID {0} does not exist!".format(job_id))

        job_name = db_entry["job"]

        h5_path = None
        sub_job = db_entry["subjob"]

        if sub_job is not None:
            if len(sub_job.strip()) > 0:
                h5_path = '/'.join(sub_job.split('/')[:-1])
        hdf5_file = sub_job.split('/')[1] + '.h5'

        gp = GenericPath(root_path=db_entry["projectpath"], project_path=db_entry["project"])
        hdf_project = ProjectHDFio(project=Project(path=gp, user=user), file_name=hdf5_file, h5_path=h5_path, mode="r")
        super(JobPath, self).__init__(hdf_project, job_name)

        self.__name__ = db_entry["hamilton"]
        self.__version__ = db_entry["hamversion"]

        if 'id' in db_entry:
            self._job_id = db_entry["id"]
        self._status = db_entry["status"]
        self._master_id = db_entry["masterid"]
        self._parent_id = db_entry["parentid"]

    @property
    def is_root(self):
        """
        Check if the current h5_path is pointing to the HDF5 root group.

        Returns:
            bool: [True/False]
        """
        return self.project_hdf5.is_root

    # @property
    # def is_open(self):
    #     """
    #     Check if the HDF5 file is currently opened in h5py
    #
    #     Returns:
    #         bool: [True/False]
    #     """
    #     return self.project_hdf5.is_open

    @property
    def is_empty(self):
        """
        Check if the HDF5 file is empty

        Returns:
            bool: [True/False]
        """
        return self.project_hdf5.is_empty

    @property
    def base_name(self):
        """
        Name of the HDF5 file - but without the file extension .h5

        Returns:
            str: file name without the file extension
        """
        return self.project_hdf5.base_name

    @property
    def file_path(self):
        """
        Path where the HDF5 file is located - posixpath.dirname()

        Returns:
            str: HDF5 file location
        """
        return self.project_hdf5.file_path

    @property
    def h5_path(self):
        """
        Get the path in the HDF5 file starting from the root group - meaning this path starts with '/'

        Returns:
            str: HDF5 path
        """
        return self.project_hdf5.h5_path

    @h5_path.setter
    def h5_path(self, path):
        """
        Set the path in the HDF5 file starting from the root group

        Args:
            path (str): HDF5 path
        """
        self.project_hdf5.h5_path = path

    def create_group(self, name):
        """
        Create an HDF5 group - similar to a folder in the filesystem - the HDF5 groups allow the users to structure their
        data.

        Args:
            name (str): name of the HDF5 group

        Returns:
            FileHDFio: FileHDFio object pointing to the new group
        """
        return self.project_hdf5.create_group(name)

    def open(self, h5_rel_path):
        """
        Create an HDF5 group and enter this specific group. If the group exists in the HDF5 path only the h5_path is
        set correspondingly otherwise the group is created first.

        Args:
            h5_rel_path (str): relative path from the current HDF5 path - h5_path - to the new group

        Returns:
            FileHDFio: FileHDFio object pointing to the new group
        """
        return self.project_hdf5.open(h5_rel_path)

    def close(self):
        """
        Close the current HDF5 path and return to the path before the last open
        """
        self.project_hdf5.close()

    def remove_file(self):
        """
        Remove the HDF5 file with all the related content
        """
        self.project_hdf5.remove_file()

    def put(self, key, value):
        """
        Store data inside the HDF5 file

        Args:
            key (str): key to store the data
            value (pandas.DataFrame, pandas.Series, dict, list, float, int): basically any kind of data is supported
        """
        self.project_hdf5.__setitem__(key, value)

    def listdirs(self):
        """
        equivalent to os.listdirs (consider groups as equivalent to dirs)

        Returns:
            (list): list of groups in pytables for the path self.h5_path

        """
        return self.project_hdf5.list_groups()

    def list_dirs(self):
        """
        equivalent to os.listdirs (consider groups as equivalent to dirs)

        Returns:
            (list): list of groups in pytables for the path self.h5_path
        """
        return self.project_hdf5.list_groups()

    def keys(self):
        """
        List all groups and nodes of the HDF5 file - where groups are equivalent to directories and nodes to files.

        Returns:
            list: all groups and nodes
        """
        return self.project_hdf5.keys()

    def values(self):
        """
        List all values for all groups and nodes of the HDF5 file

        Returns:
            list: list of all values
        """
        return self.project_hdf5.values()

    def items(self):
        """
        List all keys and values as items of all groups and nodes of the HDF5 file

        Returns:
            list: list of sets (key, value)
        """
        return self.project_hdf5.items()

    def groups(self):
        """
        Filter HDF5 file by groups

        Returns:
            FileHDFio: an HDF5 file which is filtered by groups
        """
        return self.project_hdf5.groups()

    def nodes(self):
        """
        Filter HDF5 file by nodes

        Returns:
            FileHDFio: an HDF5 file which is filtered by nodes
        """
        return self.project_hdf5.nodes()

    def __enter__(self):
        """
        Compatibility function for the with statement
        """
        return self.project_hdf5.__enter__()

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Compatibility function for the with statement
        """
        self.project_hdf5.__exit__(exc_type=exc_type, exc_val=exc_val, exc_tb=exc_tb)

    def __setitem__(self, key, value):
        """
        Store data inside the HDF5 file

        Args:
            key (str): key to store the data
            value (pandas.DataFrame, pandas.Series, dict, list, float, int): basically any kind of data is supported
        """
        self.project_hdf5.__setitem__(key, value)

    def __delitem__(self, key):
        """
        Delete item from the HDF5 file

        Args:
            key (str): key of the item to delete
        """
        self.project_hdf5.__delitem__(key)

    def __str__(self):
        """
        Machine readable string representation

        Returns:
            str: list all nodes and groups as string
        """
        return self.project_hdf5.__str__()

    def __repr__(self):
        """
        Human readable string representation

        Returns:
            str: list all nodes and groups as string
        """
        return self.project_hdf5.__repr__()

    def __del__(self):
        """
        When the object is deleted the HDF5 file has to be closed
        """
        try:
            self.project_hdf5._store.close()
        except AttributeError:
            pass

    def __getitem__(self, item):
        """
        Get/ read data from the HDF5 file

        Args:
            item (str, slice): path to the data or key of the data object

        Returns:
            dict, list, float, int: data or data object
        """
        if item in self.list_files():
            file_name = posixpath.join(self.working_directory, "{}".format(item))
            with open(file_name) as f:
                return f.readlines()
        return self.project_hdf5.__getitem__(item)
