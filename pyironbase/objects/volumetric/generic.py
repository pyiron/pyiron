# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyironbase.objects.generic.template import PyIronObject

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class VolumetricData(PyIronObject):

    """
    A new class to handle 3-dimensional volumetric data elegantly (charge densities, electrostatic potentials etc) based
    on the numpy.ndarray instance

    Attributes:

        total_data (numpy.ndarray instance): A 3D array containing the data

    """
    def __init__(self):
        self._total_data = None

    @property
    def total_data(self):
        return self._total_data

    @total_data.setter
    def total_data(self, val):
        try:
            assert(isinstance(val, (np.ndarray, list)))
            try:
                val = np.array(val)
                shape = np.array(np.shape(val))
                assert(len(shape) == 3)
                self._total_data = val
            except AssertionError:
                raise ValueError("Attribute total_data should be a 3D array")

        except AssertionError:
            raise TypeError("Attribute total_data should be a numpy.ndarray instance or a list and "
                            "not {}".format(type(val)))

    def get_average_along_axis(self, ind=2):
        """
        Get the lateral average along a certain axis direction.

        Args:
            ind (int): Index of axis.

        Returns:
            Average total along axis
        """
        if ind == 0:
            return np.average(np.average(self._total_data, axis=1), 1)
        elif ind == 1:
            return np.average(np.average(self._total_data, axis=0), 1)
        else:
            return np.average(np.average(self._total_data, axis=0), 0)

    def to_hdf(self, hdf5, group_name="volumetric_data"):
        with hdf5.open(group_name) as hdf_vd:
            hdf_vd["TYPE"] = str(type(self))
            hdf_vd["total"] = self.total_data

    def from_hdf(self, hdf5, group_name="volumetric_data"):
        with hdf5.open(group_name) as hdf_vd:
            self.total_data = hdf_vd["total"]
