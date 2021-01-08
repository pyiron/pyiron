# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
from pyiron_base import SerialMasterBase
from pyiron_atomistic.atomistics.job.atomistic import AtomisticGenericJob

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


class GenericOutput(OrderedDict):
    def __init__(self):
        super(GenericOutput, self).__init__()


class SerialMaster(SerialMasterBase, AtomisticGenericJob):
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

        .. attribute:: version

            Version of the hamiltonian, which is also the version of the executable unless a custom executable is used.

        .. attribute:: executable

            Executable used to run the job - usually the path to an external executable.

        .. attribute:: library_activated

            For job types which offer a Python library pyiron can use the python library instead of an external executable.

        .. attribute:: server

            Server object to handle the execution environment for the job.

        .. attribute:: queue_id

            the ID returned from the queuing system - it is most likely not the same as the job ID.

        .. attribute:: logger

            logger object to monitor the external execution and internal pyiron warnings.

        .. attribute:: restart_file_list

            list of files which are used to restart the calculation from these files.

        .. attribute:: job_type

            Job type object with all the available job types: ['ExampleJob', 'SerialMaster', 'ParallelMaster', 'ScriptJob',
                                                               'ListMaster']

        .. attribute:: child_names

            Dictionary matching the child ID to the child job name.

        .. attribute:: start_job

            The first job of the series.

        .. attribute:: input

            The input of the start job - the first job of the series.
    """

    def __init__(self, project, job_name):
        super(SerialMaster, self).__init__(project, job_name=job_name)

    @property
    def structure(self):
        if self.start_job is not None:
            return self._start_job.structure
        else:
            return None

    @structure.setter
    def structure(self, basis):
        if self.start_job is not None:
            self._start_job.structure = basis
        else:
            raise ValueError(
                "A structure can only be set after a start job has been assinged."
            )

    def get_structure(self, iteration_step=-1):
        """

        Returns:

        """
        if len(self.child_ids) > 0:
            return self.project.load(self.child_ids[-1]).get_structure(
                iteration_step=iteration_step
            )
        else:
            return None
