# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import numpy as np

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class Bandstructure(object):

    def __init__(self, prec=1e-5):
        self.prec = prec
        self._structure = None
        self.bmat = None
        self.point_group = None
        self._eigenvalues = None
        self._path_dict = dict()

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, val):
        self._structure = val

    @property
    def path_dict(self):
        return self._path_dict

    @path_dict.setter
    def path_dict(self, val):
        self._path_dict = val


class BandPath(object):

    def __init__(self, bs_obj, n_points=20):
        self.bs_obj = bs_obj
        self.translate_to_pylab = {"Gamma": r"$\Gamma$", "G'": r"$\Gamma^\prime$"}
        self.q_points = list()
        self.labels = bs_obj.path_dict.keys()
        self.special_points = bs_obj.path_dict.values()
        self.n_points = n_points
        self.q_dist = np.zeros(n_points)

    def _generate_points(self):
        assert(self.bs_obj.structure is not None)
        spl_distances = list()
        for i, sp in enumerate(self.special_points):
            if i < len(self.special_points) - 1:
                x1 = np.array(self.special_points[i + 1])
                x2 = np.array(sp)
                spl_distances.append(np.linalg.norm(x2 - x1))
