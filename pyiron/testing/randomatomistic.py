# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import numpy as np
import os
import posixpath
from pyiron.base.generic.parameters import GenericParameters
from pyiron.base.job.generic import GenericJob
from pyiron.base.pyio.parser import Logstatus
from pyiron.atomistics.job.interactive import GenericInteractive

"""
Example Job class for testing the pyiron classes
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
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
        self.executable = "python -m pyiron.testing.executable"
        self._interactive_cache = {"alat": [], "count": [], "energy": []}

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        self.input.read_only = True

    # define routines that create all necessary input files
    def write_input(self):
        """
        Call routines that generate the codespecifc input files
        """
        self.input.write_file(file_name="input.inp", cwd=self.working_directory)

    # define routines that collect all output files
    def collect_output(self):
        """
        Parse the output files of the example job and store the results in the HDF5 File.
        """
        self.collect_output_log()
        self.collect_warnings()
        self.collect_logfiles()

    def collect_output_log(self, file_name="output.log"):
        """
        general purpose routine to extract output from logfile

        Args:
            file_name (str): output.log - optional
        """
        tag_dict = {
            "alat": {"arg": "0", "rows": 0},
            "count": {"arg": "0", "rows": 0},
            "energy": {"arg": "0", "rows": 0},
        }
        lf = Logstatus()
        file_name = posixpath.join(self.working_directory, file_name)
        lf.extract_file(file_name=file_name, tag_dict=tag_dict)
        with self.project_hdf5.open("output/generic") as h5:
            lf.to_hdf(h5)
            h5["energy_tot"] = np.array(h5["energy"])
            h5["volume"] = np.array(h5["alat"])

    def collect_warnings(self):
        """
        Collect the warnings if any were written to the info.log file and store them in the HDF5 file
        """
        warnings_lst = []
        with open(posixpath.join(self.working_directory, "info.log"), "r") as f:
            lines = f.readlines()
        for line in lines:
            if "WARNING" in line:
                warnings_lst.append(line.split("WARNING"))
                warnings_lst[-1][-1] = warnings_lst[-1][-1].rstrip()
        if len(warnings_lst) > 0:
            warnings_dict = {
                "Module": [warnings_lst[i][0] for i in range(len(warnings_lst))],
                "Message": [warnings_lst[i][1] for i in range(len(warnings_lst))],
            }
            print("module: ", warnings_lst[:][:])
            with self.project_hdf5.open("output"):
                self._hdf5["WARNINGS"] = warnings_dict

    def collect_logfiles(self):
        """
        Collect the errors from the info.log file and store them in the HDF5 file
        """
        errors_lst = []
        with open(posixpath.join(self.working_directory, "info.log"), "r") as f:
            lines = f.readlines()
        for line in lines:
            if "ERROR" in line:
                errors_lst.append(line)
        if len(errors_lst) > 0:
            with self.project_hdf5.open("output") as hdf_output:
                hdf_output["ERRORS"] = errors_lst

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

    def run_if_interactive(self):
        """
        Run the job as Python library and store the result in the HDF5 File.

        Returns:
            int: job ID
        """
        from pyiron.testing.executable import ExampleExecutable

        self._interactive_library = True
        self.status.running = True
        alat, count, energy = ExampleExecutable().run_lib(self.input)
        self._interactive_cache["alat"].append(alat)
        self._interactive_cache["count"].append(count)
        self._interactive_cache["energy"].append(energy)

    def interactive_close(self):
        self._interactive_library = False
        self.to_hdf()
        with self.project_hdf5.open("output") as h5:
            h5["generic/energy"] = np.array(self._interactive_cache["energy"])
            h5["generic/volume"] = np.array(self._interactive_cache["alat"])
            h5["generic/alat"] = np.array(self._interactive_cache["alat"])
            h5["generic/count"] = np.array(self._interactive_cache["count"])
            h5["generic/energy_tot"] = np.array(self._interactive_cache["energy"])
        self.project.db.item_update(self._runtime(), self._job_id)
        self.status.finished = True


class ExampleInput(GenericParameters):
    """
    Input class for the ExampleJob based on the GenericParameters class.

    Args:
        input_file_name (str): Name of the input file - optional
    """

    def __init__(self, input_file_name=None):
        super(ExampleInput, self).__init__(
            input_file_name=input_file_name, table_name="input_inp", comment_char="#"
        )

    def load_default(self):
        """
        Loading the default settings for the input file.
        """
        input_str = """\
alat  3.2     # lattice constant (would be in a more realistic example in the structure file)
alpha 0.1     # noise amplitude
a_0   3       # equilibrium lattice constant
a_1   0
a_2   1.0     # 2nd order in energy (corresponds to bulk modulus)
a_3   0.0     # 3rd order
a_4   0.0     # 4th order
count 10      # number of calls (dummy)
write_restart True
read_restart False
"""
        self.load_string(input_str)


class AtomisticExampleJob(ExampleJob, GenericInteractive):
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
        self.executable = "python -m pyiron.testing.executable"
        self.interactive_cache = {
            "cells": [],
            "energy_pot": [],
            "energy_tot": [],
            "forces": [],
            "positions": [],
            "pressures": [],
            "stress": [],
            "steps": [],
            "temperature": [],
            "indices": [],
            "computation_time": [],
            "unwrapped_positions": [],
            "atom_spin_constraints": [],
            "atom_spins": [],
            "magnetic_forces": [],
            "volume": [],
        }

    @property
    def structure(self):
        """

        Returns:

        """
        return self._structure

    def get_structure(self, iteration_step=-1, wrap_atoms=True):
        structure = super(AtomisticExampleJob, self).get_structure(
            iteration_step=iteration_step, wrap_atoms=wrap_atoms
        )
        if structure is None:
            return self.structure
        return structure

    @structure.setter
    def structure(self, structure):
        """

        Args:
            structure:

        Returns:

        """
        self._structure = structure
        if structure is not None:
            self.input["alat"] = self._structure.cell[0, 0]
            # print("set alat: {}".format(self.input["alat"]))

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(AtomisticExampleJob, self).set_input_to_read_only()
        self.input.read_only = True

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

    def run_if_interactive(self):
        """
        Run the job as Python library and store the result in the HDF5 File.

        Returns:
            int: job ID
        """
        super(AtomisticExampleJob, self).run_if_interactive()
        self.interactive_cache["cells"].append(self._structure.cell)
        self.interactive_cache["energy_pot"].append(
            self._interactive_cache["energy"][-1][-1]
        )
        self.interactive_cache["energy_tot"].append(
            self._interactive_cache["energy"][-1][-1]
        )
        self.interactive_cache["forces"].append(
            np.random.random((len(self._structure), 3))
        )
        self.interactive_cache["positions"].append(self._structure.positions)
        self.interactive_cache["pressures"].append(np.random.random((3, 3)))
        self.interactive_cache["stress"].append(
            np.random.random((len(self._structure), 3, 3))
        )
        self.interactive_cache["steps"].append(len(self.interactive_cache["steps"]))
        self.interactive_cache["temperature"].append(np.random.random())
        self.interactive_cache["indices"].append(self._structure.indices)
        self.interactive_cache["computation_time"].append(np.random.random())
        self.interactive_cache["unwrapped_positions"].append(self._structure.positions)
        self.interactive_cache["volume"].append(self._structure.get_volume())
