# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
import inspect
import time
import numpy as np
from pyiron.base.master.generic import GenericMaster, get_function_from_string
from pyiron.base.generic.parameters import GenericParameters

"""
The serial master class is a metajob consisting of a dynamic list of jobs which are executed in serial mode.
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


class SerialMasterBase(GenericMaster):
    """
    The serial master class is a metajob consisting of a dynamic list of jobs which are executed in serial mode. The job
    is derived from the GenericMaster.

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

        .. attribute:: start_job

            The first job of the series.

        .. attribute:: input

            The input of the start job - the first job of the series.
    """

    def __init__(self, project, job_name):
        self._input = GenericParameters("parameters")  # e.g. convergence goal

        super(SerialMasterBase, self).__init__(project, job_name=job_name)
        self.__name__ = "SerialMaster"
        self.__version__ = "0.3"

        self._output = GenericOutput()
        self._max_iterations = 100
        self._start_job = None
        self._run_fast = False
        self._logger.debug("run_fast: {}".format(self._run_fast))
        self._convergence_goal = None
        self._convergence_goal_qwargs = {}
        self._convergence_goal_str = None

    @property
    def start_job(self):
        """
        Get the first job of the series.

        Returns:
            GenericJob: start job
        """
        if self._start_job:
            return self._start_job
        elif len(self) > 0:
            self._start_job = self[-1]
            return self._start_job
        else:
            return None

    @start_job.setter
    def start_job(self, job):
        """
        Set the first job of the series - that is the same like appending the job.

        Args:
            job (GenericJob): start job
        """
        self.append(job)

    @property
    def ref_job(self):
        return self.start_job

    @ref_job.setter
    def ref_job(self, job):
        self.append(job)

    @property
    def input(self):
        """
        Get the input of the start job - the first job of the series.

        Returns:
            GenericParameters: input of the start job
        """
        if self.start_job:
            return self._start_job.input
        else:
            return None

    @input.setter
    def input(self, value):
        """
        Set the input of the start job - the first job of the series.

        Args:
            value (GenericParameters): input of the start job
        """
        if self.start_job:
            self._start_job.input = value
        else:
            raise ValueError(
                "Input can only be set after a start job has been assinged."
            )

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        self._input.read_only = True

    def get_initial_child_name(self):
        """
        Get name of the initial child.

        Returns:
            str: name of the initial child
        """
        return self.project.db.get_item_by_id(self.child_ids[0])["job"]

    def create_next(self, job_name=None):
        """
        Create the next job in the series by duplicating the previous job.

        Args:
            job_name (str): name of the new job - optional - default='job_<index>'

        Returns:
            GenericJob: next job
        """
        if len(self) == 0:
            raise ValueError("No job available in job list, please append a job first.")
        if len(self._job_name_lst) > len(self.child_ids):
            return self.pop(-1)
        ham_old = self.project.load(self.child_ids[-1], convert_to_object=True)

        if ham_old.status.aborted:
            ham_old.status.created = True
            return ham_old
        elif not ham_old.status.finished:
            return None
        if job_name is None:
            job_name = "_".join(
                ham_old.job_name.split("_")[:-1] + [str(len(self.child_ids))]
            )
        new_job = ham_old.restart(job_name=job_name)
        new_job.server.cores = self.server.cores
        return new_job

    def collect_output(self):
        """
        Collect the output files of the individual jobs and set the output of the last job to be the output of the
        SerialMaster - so the SerialMaster contains the same output as its last child.
        """
        ham_lst = [self.project_hdf5.inspect(child_id) for child_id in self.child_ids]
        if (
            "output" in ham_lst[0].list_groups()
            and "generic" in ham_lst[0]["output"].list_groups()
        ):
            nodes = ham_lst[0]["output/generic"].list_nodes()
            with self.project_hdf5.open("output/generic") as hh:
                for node in nodes:
                    hh[node] = np.concatenate(
                        [ham["output/generic/{}".format(node)] for ham in ham_lst],
                        axis=0,
                    )

    def collect_logfiles(self):
        """
        The collect logfiles function is required by the GenericJob class, therefore we use an empty template here.
        """
        pass

    def copy(self):
        """
        Copy the GenericJob object which links to the job and its HDF5 file

        Returns:
            GenericJob: New GenericJob object pointing to the same job
        """
        new_job = super(SerialMasterBase, self).copy()
        new_job.start_job = self.start_job
        return new_job

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the SerialMaster from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(SerialMasterBase, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self._input.from_hdf(hdf5_input)
            convergence_goal_str = hdf5_input["convergence_goal"]
            if convergence_goal_str == "None":
                self._convergence_goal = None
            else:
                self._convergence_goal_str = convergence_goal_str
                self._convergence_goal = get_function_from_string(convergence_goal_str)
                self._convergence_goal_qwargs = hdf5_input["convergence_goal_qwargs"]

    def get_from_childs(self, path):
        """
        Extract the output from all child jobs and appending it to a list

        Args:
            path (str): path inside the HDF5 files of the individual jobs like 'output/generic/volume'

        Returns:
            list: list of output from the child jobs
        """
        var_lst = []
        for child_id in self.child_ids:
            ham = self.project.load(child_id, convert_to_object=False)
            var = ham.__getitem__(path)
            var_lst.append(var)
        return np.array(var_lst)

    def iter_jobs(self, convert_to_object=True):
        """
        Iterate over the jobs within the SerialMaster

        Args:
            convert_to_object (bool): load the full GenericJob object (default) or just the HDF5 / JobCore object

        Returns:
            yield: Yield of GenericJob or JobCore
        """
        for job_id in self.child_ids:
            yield self.project.load(job_id, convert_to_object=convert_to_object)

    def run_if_interactive(self):
        pass

    def _get_job_template(self):
        self._logger.info("run serial master {}".format(self.job_info_str))
        job = self.pop(-1)
        job._master_id = self.job_id
        if self.server.new_hdf:
            job._hdf5 = self.project_hdf5.create_hdf(
                path=self.project.open(self.job_name + "_hdf5").path,
                job_name=job.job_name,
            )
        else:
            job._hdf5 = self.project_hdf5.open(job.job_name)
        self._logger.info("SerialMaster: run job {}".format(job.job_name))
        return job

    @staticmethod
    def _run_child_job(job):
        job.run()

    def _run_if_master_queue(self, job):
        job.server.run_mode.modal = True
        job.run()
        self.run_if_refresh()

    def _run_if_master_non_modal_child_non_modal(self, job):
        job.run()
        if self.master_id:
            del self

    def _run_if_master_modal_child_modal(self, job):
        job.run()
        self.run_if_refresh()

    def _run_if_master_modal_child_non_modal(self, job):
        job.run()
        while not job.status.finished and not job.status.aborted:
            job.refresh_job_status()
            time.sleep(5)
        self.run_if_refresh()

    def run_static(self, **qwargs):
        self.status.running = True
        if len(self) > len(self.child_ids):
            job = self._get_job_template()
            self.status.suspended = True
            if self.server.run_mode.queue:
                self._run_if_master_queue(job)
            elif self.server.run_mode.non_modal and job.server.run_mode.non_modal:
                self._run_if_master_non_modal_child_non_modal(job)
            elif self.server.run_mode.modal and job.server.run_mode.modal:
                self._run_if_master_modal_child_modal(job)
            elif self.server.run_mode.modal and job.server.run_mode.non_modal:
                self._run_if_master_modal_child_non_modal(job)
            else:
                raise TypeError()
        else:
            self.status.collect = True
            self.run()

    def set_goal(self, convergence_goal, **qwargs):
        """
        Set a convergence goal for the SerialMaster - this is necessary to stop the series.

        Args:
            convergence_goal (Function): the convergence goal can be any Python function, but if external packages are
                                         used like numpy they have to be imported within the function.
            **qwargs: arguments of the convergence goal function.
        """
        self._convergence_goal = convergence_goal
        self._convergence_goal_qwargs = qwargs
        self._convergence_goal_str = inspect.getsource(convergence_goal)
        if self.project_hdf5.file_exists:
            self.to_hdf()

    def show(self):
        """
        list all jobs in the SerialMaster

        Returns:
            list: list of jobs ['job', <index>, <GenericJob>]
        """
        return [["job", str(i), str(job)] for i, job in enumerate(self)]

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the SerialMaster in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(SerialMasterBase, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self._input.to_hdf(hdf5_input)
            if self._convergence_goal is not None:
                try:
                    hdf5_input["convergence_goal"] = inspect.getsource(
                        self._convergence_goal
                    )
                except IOError:
                    hdf5_input["convergence_goal"] = self._convergence_goal_str

                hdf5_input["convergence_goal_qwargs"] = self._convergence_goal_qwargs
            else:
                hdf5_input["convergence_goal"] = "None"

    def write_input(self):
        """
        Write the input files - for the SerialMaster this only contains convergence goal.
        """
        self._input.write_file(file_name="input.inp", cwd=self.working_directory)

    def __len__(self):
        """
        Length of the SerialMaster equal the number of childs appended.

        Returns:
            int: length of the SerialMaster
        """
        return len(self.child_ids + self._job_name_lst)

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
            total_lst = child_name_lst + self._job_name_lst
            item = total_lst[item]
        return self._get_item_when_str(
            item=item, child_id_lst=child_id_lst, child_name_lst=child_name_lst
        )

    def run_if_refresh(self):
        """
        Internal helper function the run if refresh function is called when the job status is 'refresh'. If the job was
        suspended previously, the job is going to be started again, to be continued.
        """
        conv_goal_exists = bool(self._convergence_goal)
        self._logger.info("Does the convergence goal exit: {}".format(conv_goal_exists))
        if not conv_goal_exists:
            self.status.collect = True
            self.run()
        else:
            subjobs_statuses = set(
                [
                    self.project.db.get_job_status(job_id=child_id)
                    for child_id in self.child_ids
                ]
            )
            if len(subjobs_statuses) == 0 or subjobs_statuses == {"finished"}:
                ham = self._convergence_goal(self, **self._convergence_goal_qwargs)
                if ham is not True:
                    self.append(ham)
                    self.to_hdf()
                    self.run_static()
                else:
                    self.status.collect = True
                    self.run()


class GenericOutput(OrderedDict):
    """
    Generic Output just a place holder to store the output of the last child directly in the SerialMaster.
    """

    def __init__(self):
        super(GenericOutput, self).__init__()
