# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
This module is supposed to be common for both electronic and phonon band structures
"""
from __future__ import print_function
import numpy as np
from numpy import transpose as tr
from numpy.linalg import inv, norm
from pyiron_base.generic.template import PyironObject

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


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
        if not (self.bs_obj.structure is not None):
            raise AssertionError()
        spl_distances = list()
        for i, sp in enumerate(self.special_points):
            if i < len(self.special_points) - 1:
                x1 = np.array(self.special_points[i + 1])
                x2 = np.array(sp)
                spl_distances.append(np.linalg.norm(x2 - x1))


class Bandstructure(PyironObject):
    translate_to_pylab = {"Gamma": r"$\Gamma$", "G'": r"$\Gamma^\prime$"}

    def __init__(self, structure=None, prec=1e-5):
        self.prec = prec
        self._structure = None
        self.bmat = None
        self.point_group = None
        self.path_type = None
        self.num_points = None
        self.q_dist = None
        self.q_points = None
        self.q_ticks = None
        self.q_labels = None
        self.q_path = None
        self.ew_list = None
        self.ev_list = None
        self._eigenvalues = None
        self._path_dict = dict()
        if structure:
            self.structure = structure

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, val):
        self._structure = val
        self.bmat = tr(inv(val.cell))
        self.point_group = val.get_spacegroup(self.prec)["Number"]
        self._assign_path()

    @property
    def path_dict(self):
        return self._path_dict

    @path_dict.setter
    def path_dict(self, val):
        self._path_dict = val

    def _assign_path(self):
        b1, b2, b3 = self.bmat
        # print "B1",b1
        # print "B2",b2
        # print "B3",b3
        point_group = self.point_group

        special_points = {
            "Gamma": [0, 0, 0],
            "G'": b2,
            "L": 0.5 * (b1 + b2 + b3),
            "K": 1.0 / 8.0 * (3 * b1 + 6 * b2 + 3 * b3),
            "U": 1.0 / 8.0 * (2 * b1 + 5 * b2 + 5 * b3),
            "X": 0.5 * (b2 + b3),
            "W": 0.25 * b1 + 0.75 * b2 + 0.5 * b3,
            "X'": 0.5 * (b1 + 2 * b2 + b3),
            # replace by correct nomenclature
            "M": 0.5 * (b1 + b2),
            "X1": 0.5 * b1,
            "X2": 0.5 * b2,
        }

        if point_group == 225:  # fcc
            path_dict = {
                "very_short": ["L", "Gamma", "X"],
                "full": [
                    "Gamma",
                    "X",
                    "U",
                    "L",
                    "Gamma",
                    "K",
                ],  # , "U", "W", "L", "K"],
                "full_20": ["G", "X", "U", "Gamma", "L"],
                "full_CM": ["G'", "X'", "K", "Gamma", "L"],
            }
        elif point_group == 229:  # bcc
            path_dict = {
                "very_short": ["Gamma", "X"],
                "full": ["Gamma", "X", "L", "W", "Gamma"],
            }
        elif point_group == 221:  # sc
            path_dict = {
                "very_short": ["L", "Gamma", "X1"],
                "full": ["Gamma", "X1", "X", "L", "Gamma"],
            }
        elif point_group == 129:
            path_dict = {
                "very_short": ["Gamma", "X1"],
                "full": ["Gamma", "X1", "M", "X2", "Gamma"],
            }
        elif point_group == 123:  # used here for 1d system
            path_dict = {
                "very_short": ["X2", "Gamma", "X1"],
                "full": ["Gamma", "X1", "M", "X2", "Gamma"],
            }
        elif point_group == 166:
            path_dict = {
                "very_short": ["L", "Gamma", "X"],
                # "full": ["Gamma", "X", "K", "Gamma", "L"],
                # "full_20": ["Gamma", "X", "K", "Gamma", "L"]
                "full": ["Gamma", "X", "Gamma", "L"],
                "full_20": ["Gamma", "X", "Gamma", "L"],
            }
        elif point_group == 167:
            path_dict = {
                "very_short": ["L", "Gamma", "X"],
                # "full": ["Gamma", "X", "K", "Gamma", "L"],
                # "full_20": ["Gamma", "X", "K", "Gamma", "L"]
                "full": ["Gamma", "X", "Gamma", "L"],
                "full_20": ["Gamma", "X", "Gamma", "L"],
            }

        elif point_group == 186:
            special_points = {
                "Gamma": [0, 0, 0],
                "A": 0.5 * b3,
                "L": 0.5 * (b1 + b3),
                "M": 0.5 * b1,
                "K": (1.0 / 0.3) * (b1 + b2),
                "H": (1.0 / 0.3) * (b1 + b2) + 0.5 * b3,
            }
            path_dict = {"full": ["A", "L", "M", "Gamma", "A"]}

        else:  # TODO: taking fcc path for testing reasons, change this later
            path_dict = {
                "very_short": ["L", "Gamma", "X"],
                # "full": ["Gamma", "X", "K", "Gamma", "L"],
                # "full_20": ["Gamma", "X", "K", "Gamma", "L"]
                "full": ["Gamma", "X", "Gamma", "L"],
                "full_20": ["Gamma", "X", "Gamma", "L"],
            }
        self.path_dict = path_dict
        self.special_points = special_points

    def get_pathes(self):
        """
        provide dictionary with all predefined Bandstructure pathes for this structure
        """
        return self.path_dict

    def get_path(self, num_points=10, path_type="very_short"):
        if path_type in self.path_dict.keys():
            q_labels = self.path_dict[path_type]
        else:
            raise (
                "path type "
                + path_type
                + " does not exist for point group "
                + str(self.point_group)
            )

        self.num_points = num_points
        self.path_type = path_type

        q_path = np.array([self.special_points[q_s] for q_s in q_labels])
        self.q_path = q_path
        #        print "q_path: ", q_path

        # get total length of q_path
        q_length = 0.0
        q_vec_list = []
        for i, q in enumerate(q_path[1:]):
            d_q = q - q_path[i]

            q_vec_list.append(d_q)
            q_length += norm(d_q)
        delta_q = q_length / num_points

        # get q-points on path
        q_point_list = []
        q_dist_list = []
        q_dist_sum = 0.0
        q_ticks = []  # indices where q is special point (for graphical output)

        count = 0
        for i, q_vec in enumerate(q_vec_list):
            # print "q_vec: ", q_vec, q_path[i]
            q_dist = norm(q_vec)
            n_points = int(q_dist / delta_q + 0.5) + 1
            delta_q_i = q_dist / n_points

            delta = 1
            if i == len(q_vec_list) - 1:
                delta = 2

            q_ticks.append(count)
            for j in range(n_points + delta):
                q_point = q_path[i] + j * q_vec / (n_points + 1)
                q_point_list.append(q_point)
                q_dist_sum += delta_q_i
                q_dist_list.append(q_dist_sum)
                if delta == 2 and j == n_points + 1:
                    q_ticks.append(count)

                count += 1

        self.q_dist = q_dist_list
        self.q_points = q_point_list
        self.q_ticks = q_ticks
        self.q_labels = q_labels

        self.ew_list = None
        self.ev_list = []

        return q_dist_list, q_point_list, [q_labels, q_ticks]

    def set_eigenvalues(self, ew_list, ev_list=None):
        if not len(ew_list) == len(self.q_points):
            print("Failed")
            raise (
                "Number of eigenvalues inconsistent with q-path: "
                + str(len(ew_list))
                + " vs "
                + str(len(self.q_points))
            )

        self.ew_list = ew_list
        self.ev_list = ev_list

    def append_eigenvalues(self, ew, ev=None):
        if self.ew_list is None:
            self.ew_list = [ew]
            if ev is not None:
                self.ev_list = [ev]
        else:
            self.ew_list = np.append(self.ew_list, values=[ew], axis=0)
            if ev is not None:
                self.ev_list = np.append(self.ev_list, values=ev, axis=0)

    def plot(self):
        import pylab as plt

        q_ticks_int = [self.q_dist[i] for i in self.q_ticks]
        q_ticks_label = self.q_labels
        for i, q in enumerate(q_ticks_label):
            if q in self.translate_to_pylab:
                q_ticks_label[i] = self.translate_to_pylab[q]
        plt.plot(self.q_dist, self.ew_list)
        plt.xticks(q_ticks_int, q_ticks_label)
        for x in q_ticks_int:
            plt.axvline(x, color="black")
        return plt

    def to_hdf(self, hdf=None, group_name=None):
        if not group_name:
            group_name = "bandstructure"
        with hdf.open(group_name) as hdf5_band:
            hdf5_band["path_type"] = self.path_type
            hdf5_band["num_points"] = self.num_points
            hdf5_band["q_dist"] = self.q_dist
            hdf5_band["q_points"] = self.q_points
            hdf5_band["q_ticks"] = self.q_ticks
            hdf5_band["q_labels"] = self.q_labels
            hdf5_band["ew"] = self.ew_list

    def from_hdf(self, hdf=None, group_name=None):
        if not group_name:
            group_name = "bandstructure"
        with hdf.open(group_name) as hdf5_band:
            self.path_type = hdf5_band["path_type"]
            self.num_points = hdf5_band["num_points"]
            self.q_dist = hdf5_band["q_dist"]
            self.q_points = hdf5_band["q_points"]
            self.q_ticks = hdf5_band["q_ticks"]
            self.q_labels = hdf5_band["q_labels"]
            self.ew_list = hdf5_band["ew"]
