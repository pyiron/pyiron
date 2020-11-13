# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.vasp.structure import write_poscar

__author__ = "Sudarsan Surendralal, Su-Hyun Yoo"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH "
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

    @staticmethod
    def gauss_f(d, fwhm=0.529177):
        """
        Generates a Gaussian distribution for a given distance and full width half maximum value

        Args:
            d (float): distance between target point and reference point
            fwhm (float): Full width half maximum in angstrom

        Returns:
            float: Gaussian reduction constant

        """
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        d2 = d * d
        return np.exp(-1 / (2 * sigma ** 2) * d2)

    @staticmethod
    def dist_between_two_grid_points(target_grid_point, n_grid_at_center, lattice, grid_shape):
        """
        Calculates the distance between a target grid point and another grid point

        Args:
            target_grid_point (numpy.ndarray/list): Target grid point
            n_grid_at_center (numpy.ndarray/list): coordinate of center of sphere
            lattice (numpy.ndarray/list): lattice vector
            grid_shape (tuple/list/numpy.ndarray): size of grid

        Returns:

            float: Distance between target grid and center of sphere in angstrom

        """
        unit_dist_in_grid = [np.sqrt(np.dot(lattice[0], lattice[0])) / grid_shape[0],
                             np.sqrt(np.dot(lattice[1], lattice[1])) / grid_shape[1],
                             np.sqrt(np.dot(lattice[2], lattice[2])) / grid_shape[2]]
        dn = np.multiply(np.subtract(target_grid_point, n_grid_at_center), unit_dist_in_grid)
        dist = np.linalg.norm(dn)
        return dist

    def spherical_average_potential(self, structure, spherical_center, rad=2, fwhm=0.529177):
        """
        Calculates the spherical average about a given point in space

        Args:
            structure (pyiron.atomistics.structure.Atoms): Input structure
            spherical_center (list/numpy.ndarray): position of spherical_center in direct coordinate
            rad (float): radius of sphere to be considered in Angstrom (recommended value: 2)
            fwhm (float): Full width half maximum of gaussian function in Angstrom (recommended value: 0.529177)

        Returns:
            float: Spherical average at the target center

        """
        grid_shape = self._total_data.shape

        # Position of center of sphere at grid coordinates
        n_grid_at_center = [int(np.ceil(spherical_center[0] * grid_shape[0])),
                            int(np.ceil(spherical_center[1] * grid_shape[1])),
                            int(np.ceil(spherical_center[2] * grid_shape[2]))]

        # Unit distance between grids
        dist_in_grid = [np.linalg.norm(structure.cell[0]) / grid_shape[0],
                        np.linalg.norm(structure.cell[1]) / grid_shape[1],
                        np.linalg.norm(structure.cell[2]) / grid_shape[2]]

        # Range of grids to be considered within the provided radius w.r.t. center of sphere
        num_grid_in_sph = [[], []]
        for i, dist in enumerate(dist_in_grid):
            num_grid_in_sph[0].append(n_grid_at_center[i] - int(np.ceil(rad / dist)))
            num_grid_in_sph[1].append(n_grid_at_center[i] + int(np.ceil(rad / dist)))

        sph_avg_tmp = []
        weight = 0
        for k in range(num_grid_in_sph[0][0], num_grid_in_sph[1][0]):
            for l in range(num_grid_in_sph[0][1], num_grid_in_sph[1][1]):
                for m in range(num_grid_in_sph[0][2], num_grid_in_sph[1][2]):
                    target_grid_point = [k, l, m]
                    dist = self.dist_between_two_grid_points(target_grid_point,
                                                             n_grid_at_center, structure.cell, grid_shape)
                    if dist <= rad:
                        sph_avg_tmp.append(
                            self._total_data[k % grid_shape[0], l % grid_shape[1], m % grid_shape[2]]
                            * self.gauss_f(dist, fwhm))
                        weight += self.gauss_f(dist, fwhm)
                    else:
                        pass
        sum_list = np.sum(sph_avg_tmp)
        sph_avg = sum_list / weight
        return sph_avg

    @staticmethod
    def dist_between_two_grid_points_cyl(target_grid_point, n_grid_at_center, lattice, grid_shape, direction_of_cyl):
        """
        Distance between a target grid point and the center of a cylinder

        Args:
            target_grid_point (numpy.ndarray/list): Target grid point
            n_grid_at_center (numpy.ndarray/list): coordinate of center of sphere
            lattice (numpy.ndarray/list): lattice vector
            grid_shape (tuple/list/numpy.ndarray): size of grid
            direction_of_cyl (int): Axis of cylinder (0 (x) or 1 (y) or 2 (z))

        Returns:
            float: Distance between target grid and in-plane center of cylinder

        """
        unit_dist_in_grid = [np.sqrt(np.dot(lattice[0], lattice[0])) / grid_shape[0],
                             np.sqrt(np.dot(lattice[1], lattice[1])) / grid_shape[1],
                             np.sqrt(np.dot(lattice[2], lattice[2])) / grid_shape[2]]
        dn = np.multiply(np.subtract(target_grid_point, n_grid_at_center), unit_dist_in_grid)
        if direction_of_cyl == 0:
            dn[0] = 0
        elif direction_of_cyl == 1:
            dn[1] = 0
        elif direction_of_cyl == 2:
            dn[2] = 0
        else:
            print("check the direction of cylindrical axis")
        dist = np.linalg.norm(dn)
        return dist

    def cylindrical_average_potential(self, structure, spherical_center, axis_of_cyl, rad=2, fwhm=0.529177):
        """
        Calculates the cylindrical average about a given point in space

        Args:
            structure (pyiron.atomistics.structure.Atoms): Input structure
            spherical_center (list/numpy.ndarray): position of spherical_center in direct coordinate
            rad (float): radius of sphere to be considered in Angstrom (recommended value: 2)
            fwhm (float): Full width half maximum of gaussian function in Angstrom (recommended value: 0.529177)
            axis_of_cyl (int): Axis of cylinder (0 (x) or 1 (y) or 2 (z))

        Returns:
            float: Cylindrical average at the target center

        """
        grid_shape = self._total_data.shape

        # Position of center of sphere at grid coordinates
        n_grid_at_center = [int(np.ceil(spherical_center[0] * grid_shape[0])),
                            int(np.ceil(spherical_center[1] * grid_shape[1])),
                            int(np.ceil(spherical_center[2] * grid_shape[2]))]

        # Unit distance between grids
        dist_in_grid = [np.linalg.norm(structure.cell[0]) / grid_shape[0],
                        np.linalg.norm(structure.cell[1]) / grid_shape[1],
                        np.linalg.norm(structure.cell[2]) / grid_shape[2]]

        # Range of grids to be considered within the provided radius w.r.t. center of sphere
        num_grid_in_cyl = [[], []]

        for i, dist in enumerate(dist_in_grid):
            if i == axis_of_cyl:
                num_grid_in_cyl[0].append(0)
                num_grid_in_cyl[1].append(grid_shape[i])
            else:
                num_grid_in_cyl[0].append(n_grid_at_center[i] - int(np.ceil(rad / dist)))
                num_grid_in_cyl[1].append(n_grid_at_center[i] + int(np.ceil(rad / dist)))

        cyl_avg_tmp = []
        weight = 0
        for k in range(num_grid_in_cyl[0][0], num_grid_in_cyl[1][0]):
            for l in range(num_grid_in_cyl[0][1], num_grid_in_cyl[1][1]):
                for m in range(num_grid_in_cyl[0][2], num_grid_in_cyl[1][2]):
                    target_grid_point = [k, l, m]
                    dist = self.dist_between_two_grid_points_cyl(target_grid_point, n_grid_at_center, structure.cell,
                                                                 grid_shape, axis_of_cyl)
                    if dist <= rad:
                        cyl_avg_tmp.append(
                            self._total_data[k % grid_shape[0], l % grid_shape[1], m % grid_shape[2]]
                            * self.gauss_f(dist, fwhm))
                        weight += self.gauss_f(dist, fwhm)
                    else:
                        pass
        sum_list = np.sum(cyl_avg_tmp)
        cyl_avg = sum_list / weight

        return cyl_avg

    def get_average_along_axis(self, ind=2):
        """
        Get the lateral average along a certain axis direction. This function is adapted from the pymatgen vasp
        VolumetricData class

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

    def to_hdf(self, hdf, group_name="volumetric_data"):
        """
        Writes the data as a group to a HDF5 file

        Args:
            hdf (pyiron_base.generic.hdfio.ProjectHDFio): The HDF file/path to write the data to
            group_name (str): The name of the group under which the data must be stored as

        """
        with hdf.open(group_name) as hdf_vd:
            hdf_vd["TYPE"] = str(type(self))
            hdf_vd["total"] = self.total_data

    def from_hdf(self, hdf, group_name="volumetric_data"):
        """
        Recreating the VolumetricData instance by reading data from the HDF5 files

        Args:
            hdf (pyiron_base.generic.hdfio.ProjectHDFio): The HDF file/path to write the data to
            group_name (str): The name of the group under which the data must be stored as

        Returns:
            pyiron.atomistics.volumetric.generic.VolumetricData: The VolumetricData instance

        """
        with hdf.open(group_name) as hdf_vd:
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
            pos_data = np.genfromtxt(lines[6: n_atoms + 6])
            if n_atoms == 1:
                pos_data = np.array([pos_data])
            atomic_numbers = np.array(pos_data[:, 0], dtype=int)
            positions = pos_data[:, 2:]
            self._atoms = Atoms(numbers=atomic_numbers, positions=positions, cell=cell)
            end_int = n_atoms + 6 + int(np.prod(grid_shape) / 6)
            data = np.genfromtxt(lines[n_atoms+6: end_int])
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
