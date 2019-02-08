# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import inspect
from pyiron.base.job.generic import GenericJob

"""
The GenericMaster is the template class for all meta jobs
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


try:
    FileExistsError = FileExistsError
except NameError:
    class FileExistsError(OSError):
        pass


class GenericMaster(GenericJob):
    """
    The GenericMaster is the template class for all meta jobs - meaning all jobs which contain multiple other jobs. It
    defines the shared functionality of the different kind of job series.

    Args:
        project (ProjectHDFio): ProjectHDFio instance which points to the HDF5 file the job is stored in
        job_name (str): name of the job, which has to be unique within the project

    Attributes:

        .. attribute:: job_name

            name of the job, which has to be unique within the project

        .. attribute:: status

            execution status of the job, can be one of the following [initialized, appended, created, submitted,
                                                                      running, aborted, collect, suspended, refresh,
                                                                      busy, finished]

        .. attribute:: job_id

            unique id to identify the job in the pyiron database

        .. attribute:: parent_id

            job id of the predecessor job - the job which was executed before the current one in the current job series

        .. attribute:: master_id

            job id of the master job - a meta job which groups a series of jobs, which are executed either in parallel
            or in serial.

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

        .. attribute:: version

            Version of the hamiltonian, which is also the version of the executable unless a custom executable is used.

        .. attribute:: executable

            Executable used to run the job - usually the path to an external executable.

        .. attribute:: library_activated

            For job types which offer a Python library pyiron can use the python library instead of an external
            executable.

        .. attribute:: server

            Server object to handle the execution environment for the job.

        .. attribute:: queue_id

            the ID returned from the queuing system - it is most likely not the same as the job ID.

        .. attribute:: logger

            logger object to monitor the external execution and internal pyiron warnings.

        .. attribute:: restart_file_list

            list of files which are used to restart the calculation from these files.

        .. attribute:: job_type

            Job type object with all the available job types: ['ExampleJob', 'SerialMaster', 'ParallelMaster',
                                                               'ScriptJob', 'ListMaster']

        .. attribute:: child_names

            Dictionary matching the child ID to the child job name.
    """
    def __init__(self, project, job_name):
        super(GenericMaster, self).__init__(project, job_name=job_name)
        self._job_name_lst = []
        self._job_object_lst = []
        self._child_id_func = None
        self._child_id_func_str = None

    @property
    def child_names(self):
        """
        Dictionary matching the child ID to the child job name

        Returns:
            dict: {child_id: child job name }
        """
        child_dict = {}
        for child_id in self.child_ids:
            child_dict[child_id] = self.project.db.get_item_by_id(child_id)["job"]
        return child_dict

    @property
    def child_ids(self):
        """
        list of child job ids - only meta jobs have child jobs - jobs which list the meta job as their master

        Returns:
            list: list of child job ids
        """
        if self._child_id_func:
            return self._child_id_func(self)
        else:
            return super(GenericMaster, self).child_ids

    def first_child_name(self):
        """
        Get the name of the first child job

        Returns:
            str: name of the first child job
        """
        return self.project.db.get_item_by_id(self.child_ids[0])['job']

    def validate_ready_to_run(self):
        """
        Validate that the calculation is ready to be executed. By default no generic checks are performed, but one could
        check that the input information is complete or validate the consistency of the input at this point.
        """
        pass

    def append(self, job):
        """
        Append a job to the GenericMaster - just like you would append an element to a list.

        Args:
            job (GenericJob): job to append
        """
        if job.job_name not in self._job_name_lst:
            self._job_name_lst.append(job.job_name)
            self._child_job_update_hdf(parent_job=self, child_job=job)
            setattr(self, job.job_name, job)
            self._job_object_lst.append(job)

    def pop(self, i):
        """
        Pop a job from the GenericMaster - just like you would pop an element from a list

        Args:
            i (int): position of the job

        Returns:
            GenericJob: job
        """
        job_to_return = self._job_object_lst[i]
        del self._job_name_lst[i]
        del self._job_object_lst[i]
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input["job_list"] = self._job_name_lst
        job_to_return.project_hdf5.remove_group()
        job_to_return.project_hdf5 = self.project_hdf5.__class__(self.project, job_to_return.job_name,
                                                                 h5_path='/' + job_to_return.job_name)
        if isinstance(job_to_return, GenericMaster):
            for sub_job in job_to_return._job_object_lst:
                self._child_job_update_hdf(parent_job=job_to_return, child_job=sub_job)
        job_to_return.status.initialized = True
        return job_to_return

    def move_to(self, project):
        """
        Move the content of the job including the HDF5 file to a new location

        Args:
            project (ProjectHDFio): project to move the job to

        Returns:
            JobCore: JobCore object pointing to the new location.
        """
        if self._job_id:
            for child_id in self.child_ids:
                child = self.project.load(child_id)
                child.move_to(project.open(self.job_name + '_hdf5'))
        super(GenericMaster, self).move_to(project)

    def copy_to(self, project=None, new_job_name=None, input_only=False, new_database_entry=True):
        """
        Copy the content of the job including the HDF5 file to a new location

        Args:
            project (ProjectHDFio): project to copy the job to
            new_job_name (str): to duplicate the job within the same porject it is necessary to modify the job name
                                - optional
            input_only (bool): [True/False] to copy only the input - default False
            new_database_entry (bool): [True/False] to create a new database entry - default True

        Returns:
            GenericJob: GenericJob object pointing to the new location.
        """
        new_generic_job = super(GenericMaster, self).copy_to(project=project, new_job_name=new_job_name,
                                                             input_only=input_only,
                                                             new_database_entry=new_database_entry)
        if new_generic_job.job_id and new_database_entry and self._job_id:
            for child_id in self.child_ids:
                child = self.project.load(child_id)
                new_child = child.copy_to(project.open(self.job_name + '_hdf5'),
                                          new_database_entry=new_database_entry)
                if new_database_entry and child.parent_id:
                    new_child.parent_id = new_generic_job.job_id
                if new_database_entry and child.master_id:
                    new_child.master_id = new_generic_job.job_id
        return new_generic_job

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the GenericMaster in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(GenericMaster, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input["job_list"] = self._job_name_lst
            if self._child_id_func is not None:
                try:
                    hdf5_input["child_id_func"] = inspect.getsource(self._child_id_func)
                except IOError:
                    hdf5_input["child_id_func"] = self._child_id_func_str
            else:
                hdf5_input["child_id_func"] = "None"
        for job in self._job_object_lst:
            job.to_hdf()

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the GenericMaster from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(GenericMaster, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            job_list_tmp = hdf5_input["job_list"]
            try:
                child_id_func_str = hdf5_input["child_id_func"]
            except ValueError:
                child_id_func_str = "None"
            if child_id_func_str == "None":
                self._child_id_func = None
            else:
                self._child_id_func_str = child_id_func_str
                self._child_id_func = self.get_function_from_string(child_id_func_str)
        for ham in job_list_tmp:
            # try:
            ham_obj = self.project_hdf5.create_object(class_name=self._hdf5[ham + '/TYPE'], project=self._hdf5,
                                                      job_name=ham)
            ham_obj.from_hdf()
            setattr(self, ham_obj.job_name, ham_obj)
            self._job_object_lst.append(ham_obj)
            self._job_name_lst.append(ham_obj.job_name)
            # self.append(ham_obj)
            # except ValueError:
            #     pass

    def set_child_id_func(self, child_id_func):
        """
        Add an external function to derive a list of child IDs - experimental feature

        Args:
            child_id_func (Function): Python function which returns the list of child IDs
        """
        self._child_id_func = child_id_func
        self.save()
        self.status.finished = True

    def get_child_cores(self):
        """
        Calculate the currently active number of cores, by summarizing all childs which are neither finished nor
        aborted.

        Returns:
            (int): number of cores used
        """
        return sum([int(db_entry['computer'].split('#')[1]) for db_entry in
                    self.project.db.get_items_dict({'masterid': self.job_id})
                    if db_entry['status'] not in ['finished', 'aborted']])

    def _child_job_update_hdf(self, parent_job, child_job):
        """

        Args:
            parent_job:
            child_job:
        """
        child_job.project_hdf5.file_name = parent_job.project_hdf5.file_name
        child_job.project_hdf5.h5_path = parent_job.project_hdf5.h5_path + '/' + child_job.job_name
        if isinstance(child_job, GenericMaster):
            for sub_job in child_job._job_object_lst:
                self._child_job_update_hdf(parent_job=child_job, child_job=sub_job)

    def _executable_activate_mpi(self):
        """
        Internal helper function to switch the executable to MPI mode
        """
        pass

    def __len__(self):
        """
        Length of the GenericMaster equal the number of childs appended.

        Returns:
            int: length of the GenericMaster
        """
        return len(self._job_object_lst)

    def __getitem__(self, item):
        """
        Get/ read data from the GenericMaster

        Args:
            item (str, slice): path to the data or key of the data object

        Returns:
            dict, list, float, int: data or data object
        """
        if isinstance(item, str):
            name_lst = item.split("/")
            if name_lst[0] in self._job_name_lst:
                child = getattr(self, name_lst[0])
                if len(name_lst) == 1:
                    return child
                else:
                    return child['/'.join(name_lst[1:])]
            return super(GenericMaster, self).__getitem__(item)
        elif isinstance(item, int):
            return self._job_object_lst[item]

    @staticmethod
    def get_function_from_string(function_str):
        """
        Convert a string of source code to a function

        Args:
            function_str: function source code

        Returns:
            function:
        """
        exec(function_str)
        return eval(function_str.split("(")[0][4:])  # get function name
