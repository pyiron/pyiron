# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.master.parallel import ParallelMaster
from pyiron_base import JobGenerator

"""
The StructureListMaster class is a parallel master consisting of a list of structures which are executed in parallel.
"""

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


class StructureJobGenerator(JobGenerator):
    """
    JobGenerator for the StructureListMaster - this class implements the functions to generate the parameter list,
    modify the individual jobs according to the parameter list and generate the new job names according to the
    parameter list.
    """

    @property
    def parameter_list(self):
        """

        Returns:
            (list): [[index(int), pyiron.atomistics.structure.atoms.Atoms], ...]
        """
        return list(enumerate(self._master.structure_lst))

    @staticmethod
    def job_name(parameter):
        """
        Generate job name for a give set of parameters

        Args:
            parameter: For the StructureListMaster the structures are simply numbered - struct_0, struct_1, ...

        Returns:
            str: job name for the next job
        """
        return "struct_" + str(parameter[0])

    def modify_job(self, job, parameter):
        """
        Modify the next job by setting the structure for the specific parameter

        Args:
            job (GenericJob): next job object to be executed
            parameter: includes the atomistic structure

        Returns:
            GenericJob:
        """
        job.structure = parameter[1]
        return job


class StructureListMaster(ParallelMaster):
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
        super(StructureListMaster, self).__init__(project, job_name)
        self.__name__ = "StructureListMaster"
        self.__version__ = "0.0.1"
        self._job_generator = StructureJobGenerator(self)
        self._structure_lst = []

    @property
    def structure_lst(self):
        return self._structure_lst

    @structure_lst.setter
    def structure_lst(self, structure_lst):
        self._structure_lst = structure_lst

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the StructureListMaster object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(StructureListMaster, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input/structures") as hdf5_input:
            for ind, struct in enumerate(self.structure_lst):
                struct.to_hdf(hdf=hdf5_input, group_name="s_" + str(ind))

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the StructureListMaster object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(StructureListMaster, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input/structures") as hdf5_input:
            self._structure_lst = [
                Atoms().from_hdf(hdf5_input, group_name)
                for group_name in sorted(hdf5_input.list_groups())
            ]

    def collect_output(self):
        """
        Implemented for compatibilty
        """
        pass

    def run_if_interactive_non_modal(self):
        """
        Implemented for compatibilty
        """
        raise TypeError

    def _run_if_busy(self):
        """
        Implemented for compatibilty
        """
        pass
