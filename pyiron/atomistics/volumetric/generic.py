# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.vasp.structure import write_poscar

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH "
    "- Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class VolumetricData(object):
    """
    A new class to handle 3-dimensional volumetric data elegantly (charge densities, electrostatic potentials etc) based
    on the numpy.ndarray instance. This module is adapted from the pymatgen vasp VolumtricData class

    http://pymatgen.org/_modules/pymatgen/io/vasp/outputs.html#VolumetricData

    Attributes:

        total_data (numpy.ndarray): A 3D array containing the data

    """

    def __init__(self):
        self._total_data = None
        self._atoms = None

    @property
    def atoms(self):
        """
        The structure related to the volumeric data

        Returns:
            pyiron.atomistics.structure.Atoms: The structure associated with the data

        """
        return self._atoms

    @atoms.setter
    def atoms(self, val):
        self._atoms = val

    @property
    def total_data(self):
        """
        numpy.ndarray: The Nx x Ny x Nz sized array for the total data
        """
        return self._total_data

    @total_data.setter
    def total_data(self, val):
        if not (isinstance(val, (np.ndarray, list))):
            raise TypeError(
                "Attribute total_data should be a numpy.ndarray instance or a list and "
                "not {}".format(type(val))
            )
        val = np.array(val)
        shape = np.array(np.shape(val))
        if not (len(shape) == 3):
            raise ValueError("Attribute total_data should be a 3D array")
        self._total_data = val

    def tmp(self, ind=2):
        print("hello world")

    def get_average_along_axis(self, ind=2):
        """
        Get the lateral average along a certain axis direction. This function is adapted from the pymatgen vasp
        VolumtricData class

        http://pymatgen.org/_modules/pymatgen/io/vasp/outputs.html#VolumetricData.get_average_along_axis

        Args:
            ind (int): Index of axis (0, 1 and 2 for the x, y, and z axis respectively)

        Returns:
            numpy.ndarray: A 1D vector with the laterally averaged values of the volumetric data
        """
        if ind == 0:
            return np.average(np.average(self._total_data, axis=1), 1)
        elif ind == 1:
            return np.average(np.average(self._total_data, axis=0), 1)
        else:
            return np.average(np.average(self._total_data, axis=0), 0)

    def to_hdf(self, hdf5, group_name="volumetric_data"):
        """
        Writes the data as a group to a HDF5 file

        Args:
            hdf5 (pyiron.base.generic.hdfio.ProjectHDFio): The HDF file/path to write the data to
            group_name (str): The name of the group under which the data must be stored as

        """
        with hdf5.open(group_name) as hdf_vd:
            hdf_vd["TYPE"] = str(type(self))
            hdf_vd["total"] = self.total_data

    def from_hdf(self, hdf5, group_name="volumetric_data"):
        """
        Recreating the VolumetricData instance by reading data from the HDF5 files

        Args:
            hdf5 (pyiron.base.generic.hdfio.ProjectHDFio): The HDF file/path to write the data to
            group_name (str): The name of the group under which the data must be stored as

        Returns:
            pyiron.atomistics.volumetric.generic.VolumetricData: The VolumetricData instance

        """
        with hdf5.open(group_name) as hdf_vd:
            self._total_data = hdf_vd["total"]

    def write_cube_file(self, filename="cube_file.cube", cell_scaling=1.0):
        """
        Write the volumetric data into the CUBE file format

        Args:
            filename (str): Filename
            cell_scaling (float): Scale the cell by this fraction

        """
        if self._atoms is None:
            raise ValueError(
                "The volumetric data object must have a valid structure assigned to it before writing "
                "to the cube format"
            )
        data = self.total_data
        n_x, n_y, _ = data.shape
        origin = np.zeros(3)
        flattened_data = np.hstack(
            [data[i, j, :] for i in range(n_x) for j in range(n_y)]
        )
        n_atoms = len(self.atoms)
        total_lines = int(len(flattened_data) / 6) * 6
        reshaped_data = np.reshape(flattened_data[0:total_lines], (-1, 6))
        last_line = [flattened_data[total_lines:]]
        head_array = np.zeros((4, 4))
        head_array[0] = np.append([n_atoms], origin)
        head_array[1:, 0] = data.shape
        head_array[1:, 1:] = self.atoms.cell / data.shape * cell_scaling
        position_array = np.zeros((len(self.atoms.positions), 5))
        position_array[:, 0] = self.atoms.get_atomic_numbers()
        position_array[:, 2:] = self.atoms.positions
        with open(filename, "w") as f:
            f.write("Cube file generated by pyiron (http://pyiron.org) \n")
            f.write("z is the fastest index \n")
            np.savetxt(f, head_array, fmt="%4d %.6f %.6f %.6f")
            np.savetxt(f, position_array, fmt="%4d %.6f %.6f %.6f %.6f")
            np.savetxt(f, reshaped_data, fmt="%.5e")
            np.savetxt(f, last_line, fmt="%.5e")

    def read_cube_file(self, filename="cube_file.cube"):
        """
        Generate data from a CUBE file

        Args:
            filename (str): Filename to parse

        """
        with open(filename, "r") as f:
            lines = f.readlines()
            n_atoms = int(lines[2].strip().split()[0])
            cell_data = np.genfromtxt(lines[3:6])
            cell_grid = cell_data[:, 1:]
            grid_shape = np.array(cell_data[:, 0], dtype=int)
            # total_data = np.zeros(grid_shape)
            cell = np.array([val * grid_shape[i] for i, val in enumerate(cell_grid)])
            pos_data = np.genfromtxt(lines[6 : n_atoms + 6])
            if n_atoms == 1:
                pos_data = np.array([pos_data])
            atomic_numbers = np.array(pos_data[:, 0], dtype=int)
            positions = pos_data[:, 2:]
            self._atoms = Atoms(numbers=atomic_numbers, positions=positions, cell=cell)
            end_int = n_atoms + 6 + int(np.prod(grid_shape) / 6)
            data = np.genfromtxt(lines[n_atoms + 6 : end_int])
            data_flatten = np.hstack(data)
            if np.prod(grid_shape) % 6 > 0:
                data_flatten = np.append(
                    data_flatten, [float(val) for val in lines[end_int].split()]
                )
            n_x, n_y, n_z = grid_shape
            self._total_data = data_flatten.reshape((n_x, n_y, n_z))

    def write_vasp_volumetric(self, filename="CHGCAR", normalize=False):
        """
        Writes volumetric data into a VASP CHGCAR format

        Args:
            filename (str): Filename of the new file
            normalize (bool): True if the data is to be normalized by the volume

        """
        write_poscar(structure=self.atoms, filename=filename)
        with open(filename, "a") as f:
            f.write("\n")
            f.write(" ".join(list(np.array(self.total_data.shape, dtype=str))))
            f.write("\n")
            _, n_y, n_z = self.total_data.shape
            flattened_data = np.hstack(
                [self.total_data[:, i, j] for j in range(n_z) for i in range(n_y)]
            )
            if normalize:
                flattened_data /= self.atoms.get_volume()
            num_lines = int(len(flattened_data) / 5) * 5
            reshaped_data = np.reshape(flattened_data[0:num_lines], (-1, 5))
            np.savetxt(f, reshaped_data, fmt="%.12f")
            if len(flattened_data) % 5 > 0:
                np.savetxt(f, [flattened_data[num_lines:]], fmt="%.12f")
