# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import inspect
import textwrap
from pyiron.base.job.generic import GenericJob

"""
The GenericMaster is the template class for all meta jobs
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
        self._job_object_dict = {}
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

    @property
    def job_object_dict(self):
        """
        internal cache of currently loaded jobs

        Returns:
            dict: Dictionary of currently loaded jobs
        """
        return self._job_object_dict

    def first_child_name(self):
        """
        Get the name of the first child job

        Returns:
            str: name of the first child job
        """
        return self.project.db.get_item_by_id(self.child_ids[0])["job"]

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
        if self.status.initialized and not job.status.initialized:
            raise ValueError(
                "GenericMaster requires reference jobs to have status initialized, rather than ",
                job.status.string
            )
        if job.server.cores >= self.server.cores:
            self.server.cores = job.server.cores
        if job.job_name not in self._job_name_lst:
            self._job_name_lst.append(job.job_name)
            self._child_job_update_hdf(parent_job=self, child_job=job)

    def pop(self, i=-1):
        """
        Pop a job from the GenericMaster - just like you would pop an element from a list

        Args:
            i (int): position of the job. (Default is last element, -1.)

        Returns:
            GenericJob: job
        """
        job_name_to_return = self._job_name_lst[i]
        job_to_return = self._load_all_child_jobs(
            self._load_job_from_cache(job_name_to_return)
        )
        del self._job_name_lst[i]
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input["job_list"] = self._job_name_lst
        job_to_return.project_hdf5.remove_group()
        job_to_return.project_hdf5 = self.project_hdf5.__class__(
            self.project, job_to_return.job_name, h5_path="/" + job_to_return.job_name
        )
        if isinstance(job_to_return, GenericMaster):
            for sub_job in job_to_return._job_object_dict.values():
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
                child.move_to(project.open(self.job_name + "_hdf5"))
        super(GenericMaster, self).move_to(project)

    def copy_to(
        self, project=None, new_job_name=None, input_only=False, new_database_entry=True
    ):
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
        new_generic_job = super(GenericMaster, self).copy_to(
            project=project,
            new_job_name=new_job_name,
            input_only=input_only,
            new_database_entry=new_database_entry,
        )
        if new_generic_job.job_id and new_database_entry and self._job_id:
            for child_id in self.child_ids:
                child = self.project.load(child_id)
                new_child = child.copy_to(
                    project=project.open(self.job_name + "_hdf5"),
                    new_database_entry=new_database_entry,
                )
                if new_database_entry and child.parent_id:
                    new_child.parent_id = new_generic_job.job_id
                if new_database_entry and child.master_id:
                    new_child.master_id = new_generic_job.job_id
        return new_generic_job

    def update_master(self):
        """
        After a job is finished it checks whether it is linked to any metajob - meaning the master ID is pointing to
        this jobs job ID. If this is the case and the master job is in status suspended - the child wakes up the master
        job, sets the status to refresh and execute run on the master job. During the execution the master job is set to
        status refresh. If another child calls update_master, while the master is in refresh the status of the master is
        set to busy and if the master is in status busy at the end of the update_master process another update is
        triggered.
        """
        master_id = self.master_id
        project = self.project
        self._logger.info("update master: {} {} {}".format(master_id, self.get_job_id(), self.server.run_mode))
        if master_id is not None:
            self._reload_update_master(
                project=project,
                master_id=master_id
            )

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
            self._to_hdf_child_function(hdf=hdf5_input)
        for job in self._job_object_dict.values():
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
            self._from_hdf_child_function(hdf=hdf5_input)
            self._job_name_lst = job_list_tmp

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
        return sum(
            [
                int(db_entry["computer"].split("#")[1])
                for db_entry in self.project.db.get_items_dict(
                    {"masterid": self.job_id}
                )
                if db_entry["status"] not in ["finished", "aborted"]
            ]
        )

    def write_input(self):
        """
        Write the input files for the external executable. This method has to be implemented in the individual
        hamiltonians.
        """
        raise NotImplementedError(
            "write procedure must be defined for derived Hamilton!"
        )

    def collect_output(self):
        """
        Collect the output files of the external executable and store the information in the HDF5 file. This method has
        to be implemented in the individual hamiltonians.
        """
        raise NotImplementedError(
            "read procedure must be defined for derived Hamilton!"
        )

    def run_if_interactive(self):
        """
        For jobs which executables are available as Python library, those can also be executed with a library call
        instead of calling an external executable. This is usually faster than a single core python job.
        """
        raise NotImplementedError(
            "This function needs to be implemented in the specific class."
        )

    def interactive_close(self):
        """
        interactive close is not implemtned for MetaJobs
        """
        pass

    def interactive_fetch(self):
        """
        interactive fetch is not implemtned for MetaJobs
        """
        pass

    def interactive_flush(self, path="generic", include_last_step=True):
        """
        interactive flush is not implemtned for MetaJobs
        """
        pass

    def run_if_interactive_non_modal(self):
        """
        Run if interactive non modal is not implemented for MetaJobs
        """
        pass

    def __len__(self):
        """
        Length of the GenericMaster equal the number of childs appended.

        Returns:
            int: length of the GenericMaster
        """
        return len(self._job_name_lst)

    def __getitem__(self, item):
        """
        Get/ read data from the GenericMaster

        Args:
            item (str, slice): path to the data or key of the data object

        Returns:
            dict, list, float, int: data or data object
        """
        child_id_lst = self.child_ids
        child_name_lst = [
            self.project.db.get_item_by_id(child_id)["job"]
            for child_id in self.child_ids
        ]
        if isinstance(item, int):
            item = self._job_name_lst[item]
        return self._get_item_when_str(
            item=item, child_id_lst=child_id_lst, child_name_lst=child_name_lst
        )

    def __getattr__(self, item):
        """
        CHeck if a job with the specific name exists

        Args:
            item (str): name of the job

        Returns:

        """
        item_from_get_item = self.__getitem__(item=item)
        if item_from_get_item is not None:
            return item_from_get_item
        else:
            raise AttributeError

    def _load_all_child_jobs(self, job_to_load):
        """
        Helper function to load all child jobs to memory - like it was done in the previous implementation

        Args:
            job_to_load (GenericJob): job to be reloaded

        Returns:
            GenericJob: job to be reloaded - including all the child jobs and their child jobs
        """
        if isinstance(job_to_load, GenericMaster):
            for sub_job_name in job_to_load._job_name_lst:
                job_to_load._job_object_dict[sub_job_name] = self._load_all_child_jobs(
                    job_to_load._load_job_from_cache(sub_job_name)
                )
        return job_to_load

    def _load_job_from_cache(self, job_name):
        """
        Helper funcction to load a job either from the _job_object_dict or from the HDF5 file

        Args:
            job_name (str): name of the job

        Returns:
            GenericJob: the reloaded job
        """
        if job_name in self._job_object_dict.keys():
            return self._job_object_dict[job_name]
        else:
            ham_obj = self.project_hdf5.create_object(
                class_name=self._hdf5[job_name + "/TYPE"],
                project=self._hdf5,
                job_name=job_name,
            )
            ham_obj.from_hdf()
            return ham_obj

    def _to_hdf_child_function(self, hdf):
        """
        Helper function to store the child function in HDF5

        Args:
            hdf: HDF5 file object
        """
        hdf["job_list"] = self._job_name_lst
        if self._child_id_func is not None:
            try:
                hdf["child_id_func"] = inspect.getsource(self._child_id_func)
            except IOError:
                hdf["child_id_func"] = self._child_id_func_str
        else:
            hdf["child_id_func"] = "None"

    def _from_hdf_child_function(self, hdf):
        """
        Helper function to load the child function from HDF5

        Args:
            hdf: HDF5 file object
        """
        try:
            child_id_func_str = hdf["child_id_func"]
        except ValueError:
            child_id_func_str = "None"
        if child_id_func_str == "None":
            self._child_id_func = None
        else:
            self._child_id_func_str = child_id_func_str
            self._child_id_func = get_function_from_string(child_id_func_str)

    def _get_item_when_str(self, item, child_id_lst, child_name_lst):
        """
        Helper function for __get_item__ when item is type string

        Args:
            item (str):
            child_id_lst (list): a list containing all child job ids
            child_name_lst (list): a list containing the names of all child jobs

        Returns:
            anything
        """
        name_lst = item.split("/")
        item_obj = name_lst[0]
        if item_obj in child_name_lst:
            child_id = child_id_lst[child_name_lst.index(item_obj)]
            if len(name_lst) > 1:
                return self.project.inspect(child_id)["/".join(name_lst[1:])]
            else:
                return self.project.load(child_id, convert_to_object=True)
        elif item_obj in self._job_name_lst:
            child = self._load_job_from_cache(job_name=item_obj)
            if len(name_lst) == 1:
                return child
            else:
                return child["/".join(name_lst[1:])]
        else:
            return super(GenericMaster, self).__getitem__(item)

    def _child_job_update_hdf(self, parent_job, child_job):
        """

        Args:
            parent_job:
            child_job:
        """
        child_job.project_hdf5.file_name = parent_job.project_hdf5.file_name
        child_job.project_hdf5.h5_path = (
            parent_job.project_hdf5.h5_path + "/" + child_job.job_name
        )
        if isinstance(child_job, GenericMaster):
            for sub_job_name in child_job._job_name_lst:
                self._child_job_update_hdf(
                    parent_job=child_job,
                    child_job=child_job._load_job_from_cache(sub_job_name),
                )
        parent_job.job_object_dict[child_job.job_name] = child_job

    def _executable_activate_mpi(self):
        """
        Internal helper function to switch the executable to MPI mode
        """
        pass

    def run_if_refresh(self):
        """
        Internal helper function the run if refresh function is called when the job status is 'refresh'. If the job was
        suspended previously, the job is going to be started again, to be continued.
        """
        raise NotImplementedError(
            "Refresh is not supported for this job type for job  " + str(self.job_id)
        )

    def _run_if_busy(self):
        """
        Run if busy is not implemented for MetaJobs
        """
        pass


def get_function_from_string(function_str):
    """
    Convert a string of source code to a function

    Args:
        function_str: function source code

    Returns:
        function:
    """
    function_dedent_str = textwrap.dedent(function_str)
    exec(function_dedent_str)
    return eval(function_dedent_str.split("(")[0][4:])
