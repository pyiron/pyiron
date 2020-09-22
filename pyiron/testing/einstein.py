# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import numpy as np
import posixpath
from pyiron_base import GenericParameters, GenericJob, Logstatus
from pyiron.atomistics.job.interactive import GenericInteractive
from pyiron.testing.randomatomistic import ExampleJob

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Osamu Waseda"
__email__ = "waseda@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2020"


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
spring_constant  1.0     # spring constant
volume_per_atom  10.0    # ideal volume per atom
bulk_modulus     100.0   # bulk modulus
displacement_mag 1       # magnitude of displacements for the ideal structure
"""
        self.load_string(input_str)

class EinsteinCrystal(object):
    def __init__(self, structure, spring_constant=0, bulk_modulus=0, volume_per_atom=0):
        self._ideal_structure = None
        self.structure = structure
        self.spring_constant = spring_constant
        self.bulk_modulus = bulk_modulus
        self.volume_per_atom = volume_per_atom
        self.energy = None
        self.forces = None
        self.pressures = None

    def create_ideal_structure(self, displacement_magnitude):
        self._ideal_structure = self.structure.copy()
        self._ideal_structure.positions += displacement_magnitude*(np.random.random((len(self.structure), 3))-0.5)
        self._ideal_structure.center_coordinates_in_unit_cell()

    def run(self):
        dr = self._ideal_structure.get_scaled_positions()
        dr -= self.structure.get_scaled_positions()
        dr = np.einsum('ij,nj->ni', self.structure.cell, dr)
        volume = self.structure.get_volume(per_atom=True)
        self.pressures = -np.eye(3)*(volume-self.volume_per_atom)/volume*self.bulk_modulus
        self.energy = 0.5*self.spring_constant*np.sum(dr**2)-self.pressures[0,0]*(volume-self.volume_per_atom)*len(dr)
        self.forces = self.spring_constant*dr


class EinsteinExampleJob(ExampleJob, GenericInteractive):
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
        super(EinsteinExampleJob, self).__init__(project, job_name)
        self.__version__ = "0.3"
        self.__name__ = "EinsteinExampleJob"
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

    def interactive_initialize_interface(self):
        self._interactive_library = EinsteinCrystal(self._structure,
                                                    self.input['spring_constant'],
                                                    self.input['volume_per_atom'],
                                                    self.input['bulk_modulus'])
        self._interactive_library.create_ideal_structure(self.input['displacement_mag'])

    @property
    def structure(self):
        """

        Returns:

        """
        return self._structure

    def get_structure(self, iteration_step=-1, wrap_atoms=True):
        structure = super(EinsteinExampleJob, self).get_structure(
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

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(EinsteinExampleJob, self).set_input_to_read_only()
        self.input.read_only = True

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(EinsteinExampleJob, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(EinsteinExampleJob, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()

    def run_if_interactive(self):
        """
        Run the job as Python library and store the result in the HDF5 File.

        Returns:
            int: job ID
        """
        super(EinsteinExampleJob, self).run_if_interactive()
        self._interactive_library.run()
        self.interactive_collect()
