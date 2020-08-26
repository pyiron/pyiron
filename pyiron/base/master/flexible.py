# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import inspect
from pyiron.base.master.generic import GenericMaster
from pyiron.base.job.jobstatus import job_status_finished_lst

"""
The Flexible master uses a list of functions to connect multiple jobs in a series.
"""

__author__ = "Jan Janssen, Liam Huber"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Mar 24, 2019"


class FlexibleMaster(GenericMaster):
    """
    The FlexibleMaster uses a list of functions to connect multiple jobs in a series.

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
        super(FlexibleMaster, self).__init__(project, job_name=job_name)
        self.__name__ = "FlexibleMaster"
        self.__version__ = "0.1"
        self._step_function_lst = []

    @property
    def function_lst(self):
        return self._step_function_lst

    def validate_ready_to_run(self):
        """
        Validate that the calculation is ready to be executed. By default no generic checks are performed, but one could
        check that the input information is complete or validate the consistency of the input at this point.
        """
        super(FlexibleMaster, self).validate_ready_to_run()
        if len(self._job_name_lst) != len(self._step_function_lst) + 1:
            raise ValueError

    def is_finished(self):
        """
        Check if the ParallelMaster job is finished - by checking the job status and the submission status.

        Returns:
            bool: [True/False]
        """
        if self.status.finished:
            return True
        if len(self._job_name_lst) > 0:
            return False
        return self.check_all_childs_finished()

    def check_all_childs_finished(self):
        return set(
            [
                self.project.db.get_job_status(job_id=child_id)
                for child_id in self.child_ids
            ]
        ) < set(job_status_finished_lst)

    def run_static(self):
        """
        The FlexibleMaster uses functions to connect multiple Jobs.
        """
        self.status.running = True
        max_steps = len(self.child_ids + self._job_name_lst)
        ind = max_steps - 1
        if self.check_all_childs_finished():
            for ind in range(len(self.child_ids), max_steps):
                job = self.pop(0)
                job._master_id = self.job_id
                if ind != 0:
                    prev_job = self[ind - 1]
                    if ind < max_steps:
                        mod_funct = self._step_function_lst[ind - 1]
                        mod_funct(prev_job, job)
                    job._parent_id = prev_job.job_id
                job.run()
                if job.server.run_mode.interactive and not isinstance(job, GenericMaster):
                    job.interactive_close()
                if self.server.run_mode.non_modal and job.server.run_mode.non_modal:
                    break
                if job.server.run_mode.queue:
                    break
        if ind == max_steps - 1 and self.is_finished():
            self.status.finished = True
            self.project.db.item_update(self._runtime(), self.job_id)
        else:
            self.status.suspended = True

    def run_if_refresh(self):
        """
        Internal helper function the run if refresh function is called when the job status is 'refresh'. If the job was
        suspended previously, the job is going to be started again, to be continued.
        """
        if self.is_finished():
            self.status.collect = True
            self.run()  # self.run_if_collect()
        elif self.server.run_mode.non_modal or self.server.run_mode.queue or self.server.run_mode.modal:
            self.run_static()
        else:
            self.refresh_job_status()
            if self.status.refresh:
                self.status.suspended = True
            if self.status.busy:
                self.status.refresh = True
                self.run_if_refresh()

    def write_input(self):
        """
        write_input is not implemented for FlexibleMaster jobs
        """
        pass

    def collect_output(self):
        """
        Collect output is not implemented for FlexibleMaster jobs
        """
        pass

    def run_if_interactive(self):
        """
        run_if_interactive() is not implemented for FlexibleMaster jobs
        """
        pass

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the FlexibleMaster in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(FlexibleMaster, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if self._step_function_lst is not []:
                try:
                    hdf5_input["funct_lst"] = [
                        inspect.getsource(funct) for funct in self._step_function_lst
                    ]
                except IOError:
                    pass

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the FlexibleMaster from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(FlexibleMaster, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if "funct_lst" in hdf5_input.list_nodes() and self._step_function_lst == []:
                funct_str_lst = hdf5_input["funct_lst"]
                for funct_str in funct_str_lst:
                    exec(funct_str)
                    self._step_function_lst.append(eval(funct_str.split("(")[0][4:]))

    def __getitem__(self, item):
        child_id_lst = self.child_ids
        child_name_lst = [
            self.project.db.get_item_by_id(child_id)["job"]
            for child_id in self.child_ids
        ]
        if isinstance(item, int):
            total_lst = child_name_lst + self._job_name_lst
            item = total_lst[item]
        return self._get_item_when_str(
            item=item, child_id_lst=child_id_lst, child_name_lst=child_name_lst
        )
