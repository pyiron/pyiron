# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import os
from pyiron.base.generic.parameters import GenericParameters
from pyiron.base.job.generic import GenericJob
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.testing.interface import InteractiveRandomInterface

"""
Example Job class for testing the pyiron classes  
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class ExampleJob(GenericJob):
    """
    ExampleJob generating a list of random numbers to simulate energy fluctuations.

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
    """

    def __init__(self, project, job_name):
        super(ExampleJob, self).__init__(project, job_name)
        self.__version__ = "0.3"
        self.__name__ = "ExampleJob"
        self.input = ExampleInput()
        self.executable = "python " + str(os.path.dirname(os.path.realpath(__file__))) + \
                          "/executable.py"
        self._interface = InteractiveRandomInterface(job=self)

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(ExampleJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(ExampleJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.from_hdf(hdf5_input)


class ExampleInput(GenericParameters):
    """
    Input class for the ExampleJob based on the GenericParameters class.

    Args:
        input_file_name (str): Name of the input file - optional
    """

    def __init__(self, input_file_name=None):
        super(ExampleInput, self).__init__(input_file_name=input_file_name,
                                           table_name="input_inp",
                                           comment_char="#")

    def load_default(self):
        """
        Loading the default settings for the input file.
        """
        input_str = '''\
alat  3.2     # lattice constant (would be in a more realistic example in the structure file)
alpha 0.1     # noise amplitude
a_0   3       # equilibrium lattice constant
a_1   0
a_2   1.      # 2nd order in energy (corresponds to bulk modulus)
a_3   0.      # 3rd order
a_4   0.      # 4th order
count 10      # number of calls (dummy)
write_restart True
read_restart False
'''
        self.load_string(input_str)


class AtomisticExampleJob(AtomisticGenericJob, ExampleJob):
    """
    ExampleJob generating a list of random numbers to simulate energy fluctuations.

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
    """
    def __init__(self, project, job_name):
        super(AtomisticExampleJob, self).__init__(project, job_name)
        self.__version__ = "0.3"
        self.__name__ = "AtomisticExampleJob"
        self.input = ExampleInput()
        self.executable = "python " + str(os.path.dirname(os.path.realpath(__file__))) + \
                          "/executable.py"

    @property
    def structure(self):
        """
        
        Returns:

        """
        return self._structure

    @structure.setter
    def structure(self, structure):
        """
        
        Args:
            structure: 

        Returns:

        """
        self._structure = structure
        if self._structure.cell.any():
            self.input["alat"] = self._structure.cell[0, 0]
            # print("set alat: {}".format(self.input["alat"]))

    def get_structure(self, iteration_step=-1):
        return self.structure

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(AtomisticExampleJob, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(AtomisticExampleJob, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()
