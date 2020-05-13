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
        super(GenericInteractive, self).__init__(project, job_name)
        self.__version__ = "0.0.1"
        self.__name__ = "HessianJob"
        self.interactive_cache = {
            "forces": [],
            "positions": [],
            "energy_pot": []
        }
        self._force_constants = None
        self._reference_structure = None
        self._next_positions = None

    def set_force_constants(self, force_constants):
        force_shape = np.shape(force_constants)
        force_reshape = force_shape[0] * force_shape[2]
        self._force_constants = np.transpose(force_constants, (0, 2, 1, 3)).reshape((force_reshape, force_reshape))

    def set_reference_structure(self, structure):
        self._reference_structure = structure
        self.structure = structure

    def interactive_position_setter(self, positions):
        positions -= self._reference_structure.positions
        positions[positions > self._reference_structure.cell[0, 0] / 2] -= self._reference_structure.cell[0, 0]
        self._next_positions = positions

    def calculate_forces(self):
        self.interactive_cache["positions"].append(self._next_positions)
        position_transformed = self._next_positions.reshape(
            self._next_positions.shape[0] * self._next_positions.shape[1])
        forces_transformed = -np.dot(self._force_constants, position_transformed)
        self.interactive_cache["forces"].append(forces_transformed.reshape(self._next_positions.shape))
        self.interactive_cache["energy_pot"].append(-1 / 2 * np.dot(position_transformed, forces_transformed))
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
        self._interactive_library = False
        self.to_hdf()
        with self.project_hdf5.open("output") as h5:
            h5["generic/forces"] = np.array(self.interactive_cache["forces"])
            h5["generic/energy_pot"] = np.array(self.interactive_cache["energy_pot"])
            h5["generic/positions"] = np.array(self.interactive_cache["positions"])
        self.project.db.item_update(self._runtime(), self._job_id)
        self.status.finished = True
