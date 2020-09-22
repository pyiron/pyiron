# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import numpy as np
from pyiron_base import GenericParameters
from pyiron.atomistics.job.interactive import GenericInteractive

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
    def __init__(self, structure, spring_constant=0, volume_per_atom=0, bulk_modulus=0):
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
        dr -= np.rint(dr)
        dr = np.einsum('ji,nj->ni', self.structure.cell, dr)
        volume = self.structure.get_volume(per_atom=True)
        self.pressures = (volume-self.volume_per_atom)/volume
        self.pressures = -np.eye(3)*self.pressures*self.bulk_modulus
        self.energy = 0.5*self.spring_constant*np.sum(dr**2)-self.pressures[0,0]*(volume-self.volume_per_atom)*len(dr)
        self.forces = self.spring_constant*dr


class EinsteinExampleJob(GenericInteractive):
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
        self.__version__ = "0.0.1"
        self.__name__ = "EinsteinExampleJob"
        self.input = ExampleInput()
        self.server.run_mode.interactive = True

    def interactive_initialize_interface(self):
        self._interactive_library = EinsteinCrystal(self.structure,
                                                    self.input['spring_constant'],
                                                    self.input['volume_per_atom'],
                                                    self.input['bulk_modulus'])
        self._interactive_library.create_ideal_structure(self.input['displacement_mag'])

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

    def write_input(self):
        """
        Call routines that generate the codespecifc input files
        """
        self.input.write_file(file_name="input.inp", cwd=self.working_directory)

    def interactive_pressures_getter(self):
        return self._interactive_library.pressures

    def interactive_cells_setter(self, cell):
        self._interactive_library.structure.cell = cell

    def interactive_energy_pot_getter(self):
        return self._interactive_library.energy

    def interactive_energy_tot_getter(self):
        return self._interactive_library.energy

    def interactive_forces_getter(self):
        return self._interactive_library.forces

    def interactive_positions_setter(self, positions):
        self._interactive_library.structure.positions = positions

    def interactive_structure_setter(self, structure):
        self._interactive_library.structure = structure

    def interactive_indices_setter(self, indices):
        self._interactive_library.structure.indices = indices

    def interactive_volume_getter(self):
        return self.structure.get_volume()

    def interactive_close(self):
        if self.interactive_is_activated():
            super(EinsteinExampleJob, self).interactive_close()
            with self.project_hdf5.open("output") as h5:
                if "interactive" in h5.list_groups():
                    for key in h5["interactive"].list_nodes():
                        h5["generic/" + key] = h5["interactive/" + key]
