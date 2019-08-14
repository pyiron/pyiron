# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH " \
                "- Computational Materials Design (CM) Department"
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

    @property
    def total_data(self):
        """
        numpy.ndarray: The Nx x Ny x Nz sized array for the total data
        """
        return self._total_data

    @total_data.setter
    def total_data(self, val):
        if not (isinstance(val, (np.ndarray, list))):
            raise TypeError("Attribute total_data should be a numpy.ndarray instance or a list and "
                            "not {}".format(type(val)))
        val = np.array(val)
        shape = np.array(np.shape(val))
        if not (len(shape) == 3):
            raise ValueError("Attribute total_data should be a 3D array")
        self._total_data = val

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
            self.total_data = hdf_vd["total"]
