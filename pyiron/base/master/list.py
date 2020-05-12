# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron.base.generic.parameters import GenericParameters
from pyiron.base.job.core import JobCore
from pyiron.base.job.generic import GenericJob
from pyiron.base.master.generic import GenericMaster
from pyiron.base.master.submissionstatus import SubmissionStatus

"""
The ListMaster behaves like a list, just for job objects.
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


class ListMaster(GenericMaster):
    """
    The ListMaster is the most simple MetaJob derived from the GenericMaster. It behaves like a Python list object. Jobs
    can be append to the ListMaster just like elements are added to a list and then all jobs can be executed together.
    This also works for already executed jobs, unless they are already linked to a different MetaJob - meaning they
    already have a master ID assigned to them.

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

        .. attribute:: submission_status

            Monitors how many jobs have been submitted and how many have to be submitted in future.
    """

    def __init__(self, project, job_name):
        self._input = GenericParameters("parameters")
        super(ListMaster, self).__init__(project, job_name=job_name)
        self.__name__ = "ListMaster"
        self.__version__ = "0.1"
        self._input["mode"] = "parallel"
        self.submission_status = SubmissionStatus(db=project.db, job_id=self.job_id)
        self.refresh_submission_status()

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        self._input.read_only = True

    def reset_job_id(self, job_id=None):
        """
        Reset the job id sets the job_id to None as well as all connected modules like JobStatus and SubmissionStatus.
        """
        super(ListMaster, self).reset_job_id(job_id=job_id)
        self.submission_status = SubmissionStatus(db=self.project.db, job_id=job_id)

    def refresh_submission_status(self):
        """
        Refresh the submission status - if a job ID job_id is set then the submission status is loaded from the
        database.
        """
        if self.job_id:
            self.submission_status = SubmissionStatus(
                db=self._hdf5.db, job_id=self.job_id
            )
            self.submission_status.refresh()

    def save(self):
        """
        Save the object, by writing the content to the HDF5 file and storing an entry in the database.

        Returns:
            (int): Job ID stored in the database
        """
        job_id = super(ListMaster, self).save()
        self.refresh_submission_status()
        return job_id

    def append(self, job):
        """
        Append a job to the ListMaster - just like you would append an element to a list.

        Args:
            job (JobCore, GenericJob, int): job to append
        """
        if isinstance(job, JobCore) and job.job_id:
            job = job.job_id
        if isinstance(job, int):
            job = self._hdf5.project.load(job)
        if isinstance(job, GenericJob):
            if job.status.created or job.status.initialized:
                super(ListMaster, self).append(job=job)
            elif job.job_id and job.status.finished:
                if not job.master_id:
                    if not self.job_id:
                        self._job_id = self.save()
                    child_db_entry = self.project.db.get_item_by_id(job.job_id)
                    self.project.db.delete_item(job.job_id)
                    job._job_id = None
                    job.master_id = self._job_id
                    job.save()
                    del child_db_entry["id"]
                    del child_db_entry["masterid"]
                    self.project.db.item_update(child_db_entry, job.job_id)
                    self.submission_status.submit_next()
                    if len(self._job_name_lst) == 0:
                        self.status.finished = True
                        self.project.db.item_update(self._runtime(), self.job_id)
                else:
                    raise ValueError(
                        "This job ",
                        job.job_name,
                        " is already connected to a master ",
                        job.master_id,
                        " and can not be appended here.",
                    )
        else:
            raise TypeError(
                "job has to be either GenericJob, JobCore or int, but it - ",
                job,
                " is ",
                type(job),
            )

    def is_finished(self):
        """
        Check if the ListMaster job is finished - by checking the job status and the submission status.

        Returns:
            bool: [True/False]
        """
        if self.status.finished:
            return True
        # self.status.busy = True
        self.submission_status.refresh()
        if not self.submission_status.finished:
            return False
        else:
            status_set = set(
                [
                    self._hdf5.db.get_item_by_id(child_id)["status"]
                    for child_id in self.child_ids
                ]
            )
            # status_set = set([job.get_status() for job in self.iter_jobs(convert_to_object=False)])
            if "finished" in status_set:
                return len(status_set) == 1
            else:
                return False

    def run_static(self):
        """
        The run static function is called by run to execute the simulation. For the
        ListMaster this means executing all the childs appened in parallel.
        """
        self._input["num_points"] = len(self)
        self._logger.info("{} run parallel master (modal)".format(self.job_info_str))
        self.status.running = True
        if len(self._job_name_lst) > 0:
            job_lst = []
            for i in range(len(self._job_name_lst)):
                ham = self.pop(i=0)
                if (
                    ham.server.run_mode.non_modal
                    and self.get_child_cores() + ham.server.cores > self.server.cores
                ):
                    break
                self.submission_status.submit_next()
                if not ham.status.finished:
                    ham.run()
                self._logger.info("ListMaster: finished job {}".format(ham.job_name))
                if ham.server.run_mode.thread:
                    job_lst.append(ham._process)
                else:
                    self.refresh_job_status()
            _ = [process.communicate() for process in job_lst if process]
            self.status.suspended = True
        if self.server.run_mode.modal or (
            (self.server.run_mode.non_modal or self.server.run_mode.queue)
            and self.is_finished()
        ):
            self.status.finished = True

    def write_input(self):
        """
        Write the input files - for the ListMaster this only contains the execution mode, which is 'parallel' by
        default.
        """
        self._input.write_file(file_name="input.inp", cwd=self.working_directory)

    def copy(self):
        """
        Copy the ListMaster object which links to the job and its HDF5 file

        Returns:
            ListMaster: New ListMaster object pointing to the same job
        """
        new_job = super(ListMaster, self).copy()
        new_job._child_ids = self.child_ids[:]
        return new_job

    def iter_jobs(self, convert_to_object=True):
        """
        Iterate over the jobs within the ListMaster

        Args:
            convert_to_object (bool): load the full GenericJob object (default) or just the HDF5 / JobCore object

        Returns:
            yield: Yield of GenericJob or JobCore
        """
        for job_id in self.child_ids:
            yield self._hdf5.load(job_id, convert_to_object=convert_to_object)

    def run_if_refresh(self):
        """
        Internal helper function the run if refresh function is called when the job status is 'refresh'. If the job was
        suspended previously, the job is going to be started again, to be continued.
        """
        self._logger.info(
            "{}, status: {}, finished: {} parallel master "
            "refresh".format(self.job_info_str, self.status, self.is_finished())
        )
        if self.is_finished() and not self.server.run_mode.modal:
            self.status.finished = True
        elif (
            self.server.run_mode.non_modal or self.server.run_mode.queue
        ) and not self.submission_status.finished:
            self.run_static()
        else:
            self.refresh_job_status()
            if self.status.refresh:
                self.status.suspended = True

    def collect_output(self):
        """
        Collect output is not implemented for ListMaster jobs
        """
        pass

    def run_if_interactive(self):
        """
        run_if_interactive() is not implemented for ListMaster jobs
        """
        pass

    def __getitem__(self, item):
        """
        Get/ read data from the HDF5 file

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
            total_lst = child_name_lst + self._job_name_lst
            item = total_lst[item]
        return self._get_item_when_str(
            item=item, child_id_lst=child_id_lst, child_name_lst=child_name_lst
        )

    def __len__(self):
        """
        Length of the ListMaster equal the number of childs appended.

        Returns:
            int: length of the ListMaster
        """
        return len(self.child_ids + self._job_name_lst)
