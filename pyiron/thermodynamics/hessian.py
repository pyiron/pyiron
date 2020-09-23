# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.job.interactive import GenericInteractive

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Feb 20, 2020"


class HessianJob(GenericInteractive):
    def __init__(self, project, job_name):
        super(HessianJob, self).__init__(project, job_name)
        self.__version__ = "0.0.1"
        self.__name__ = "HessianJob"
        self.server.run_mode.interactive = True
        self._force_constants = None
        self._reference_structure = None
        self._next_positions = None
        self._energy = None
        self._forces = None

    def set_force_constants(self, force_constants):
        if self.structure is None:
            raise ValueError('Set reference structure via set_reference_structure() first')
        n_atom = len(self.structure.positions)
        if len(np.array([force_constants]).flatten())==1:
            self._force_constants = force_constants*np.eye(3*n_atom)
        elif np.array(force_constants).shape==(3*n_atom, 3*n_atom):
            self._force_constants = force_constants
        elif np.array(force_constants).shape==(n_atom, n_atom):
            na = np.newaxis
            self._force_constants = (np.array(force_constants)[:,na,:,na]*np.eye(3)[na,:,na,:]).flatten()
        elif len(np.shape(force_constants))==4:
            force_shape = np.shape(force_constants)
            if force_shape[2]==3 and force_shape[3]==3:
                force_reshape = force_shape[0] * force_shape[2]
                self._force_constants = np.transpose(force_constants, (0, 2, 1, 3)).reshape((force_reshape, force_reshape))
            elif force_shape[1]==3 and force_shape[3]==3:
                self._force_constants = np.array(force_constants).reshape(3*n_atom, 3*n_atom)
            else:
                raise AssertionError('force constant shape not recognized')
        else:
            raise AssertionError('force constant shape not recognized')

    def set_reference_structure(self, structure):
        self._reference_structure = structure.copy()
        if self.structure is None:
            self.structure = structure.copy()

    def interactive_position_setter(self, positions):
        self.structure.positions = positions.copy()
        self._next_positions = self.structure.get_scaled_positions()
        self._next_positions -= self._reference_structure.get_scaled_positions()
        self._next_positions -= np.rint(self._next_positions)
        self._next_positions = np.einsum('ji,ni->nj', self.structure.cell, self._next_positions)
        self._next_positions = positions

    def validate_ready_to_run(self):
        super(HessianJob, self).validate_ready_to_run()
        if self._force_constants is None:
            raise AssertionError('set force constants by set_force_constants before run')
        if self._reference_structure is None:
            raise AssertionError('set reference structure by set_reference_structure before run')

    def interactive_forces_getter(self):
        return self._forces

    def interactive_energy_pot_getter(self):
        return self._energy

    def interactive_positions_getter(self):
        return self.structure.positions

    def interactive_cells_getter(self):
        return self.structure.cell

    def interactive_volume_getter(self):
        return self.structure.get_volume()

    def calculate_forces(self):
        self.interactive_cache["displacements"].append(self._next_positions)
        position_transformed = self._next_positions.reshape(
            self._next_positions.shape[0] * self._next_positions.shape[1])
        forces_transformed = -np.dot(self._force_constants, position_transformed)
        self._forces = forces_transformed.reshape(self._next_positions.shape)
        self._energy_pot = -1 / 2 * np.dot(position_transformed, forces_transformed)
        self._next_positions = None

    def run_if_interactive(self):
        """
        Run the job as Python library and store the result in the HDF5 File.

        Returns:
            int: job ID
        """
        self._interactive_library = True
        self.status.running = True
        self.interactive_position_setter(positions=self.structure.positions)
        self.calculate_forces()
        self.interactive_collect()

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(HessianJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if self._force_constants is not None:
                hdf5_input["force_constants"] = self._force_constants
            if self._reference_structure is not None:
                self._reference_structure.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(HessianJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if "structure" in hdf5_input.list_groups():
                self._reference_structure = Atoms().from_hdf(hdf5_input)
            if "force_constants" in hdf5_input.list_nodes():
                self._force_constants = hdf5_input["force_constants"]

    def interactive_close(self):
        if self.interactive_is_activated():
            super(HessianJob, self).interactive_close()
            with self.project_hdf5.open("output") as h5:
                if "interactive" in h5.list_groups():
                    for key in h5["interactive"].list_nodes():
                        h5["generic/" + key] = h5["interactive/" + key]

