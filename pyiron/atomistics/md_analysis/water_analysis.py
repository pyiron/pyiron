# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, unicode_literals
from pyiron.atomistics.md_analysis.trajectory_analysis import TrajectoryAnalysis
import numpy as np

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class WaterTrajectory(TrajectoryAnalysis):

    def __init__(self):
        super(WaterTrajectory, self).__init__()
        self.oxygen_indices = None
        self.hydrogen_indices = None
        self.h1_indices = None
        self.h2_indices = None
        self._bond_cutoff = 1.3
        self.water_oxygen_indices = None
        self.hydroxyl_hydrogen_indices = None
        self.hydroxyl_oxygen_indices = None
        self.hydroxyl_indices = None
        self.attached_hydrogen_indices = None
        self.proton_indices = None

    @property
    def contains_water(self):

        return True

    @property
    def bond_cutoff(self):
        return self._bond_cutoff

    @bond_cutoff.setter
    def bond_cutoff(self, val):
        self._bond_cutoff = val
        self._set_trajectory()

    def _set_trajectory(self):
        super(WaterTrajectory, self)._set_trajectory()
        self.set_neighborhood()

    def set_neighborhood(self, radius=5.0):
        self.oxygen_indices = self.structure.select_index("O")
        self.hydrogen_indices = self.structure.select_index("H")
        water_oxygens = list()
        hydroxyl_hydrogen_indices = list()
        hydroxyl_oxygen_indices = list()
        hydronium_hydrogen_indices = list()
        hydronium_oxygen_indices = list()
        neighbors = self.structure.get_neighbors(radius=radius, cutoff=self.bond_cutoff)
        h1_indices = list()
        h2_indices = list()
        for i, ind_list in enumerate(np.array(neighbors.indices)[self.oxygen_indices]):
            h_ind_list = np.intersect1d(ind_list, self.hydrogen_indices)
            n_hydrogen_ions = len(h_ind_list)
            if n_hydrogen_ions == 2:
                h1_indices.append(h_ind_list[0])
                h2_indices.append(h_ind_list[1])
                water_oxygens.append(self.oxygen_indices[i])
            if n_hydrogen_ions == 1:
                hydroxyl_hydrogen_indices.append(h_ind_list[-1])
                hydroxyl_oxygen_indices.append(self.oxygen_indices[i])
            if n_hydrogen_ions == 3:
                hydronium_hydrogen_indices.append(h_ind_list)
                hydronium_oxygen_indices.append(self.oxygen_indices[i])


        self.h1_indices = np.array(h1_indices)
        self.h2_indices = np.array(h2_indices)
        self.water_oxygen_indices = np.array(water_oxygens)
        self.hydroxyl_hydrogen_indices = np.array(hydroxyl_hydrogen_indices)
        self.hydroxyl_oxygen_indices = np.array(hydroxyl_oxygen_indices)
        self.hydroxyl_indices = np.union1d(self.hydroxyl_oxygen_indices, self.hydroxyl_hydrogen_indices)
        self.attached_hydrogen_indices = np.union1d(np.union1d(self.h1_indices, self.h2_indices),
                                                    self.hydroxyl_hydrogen_indices)
        self.proton_indices = np.setdiff1d(self.hydrogen_indices, self.attached_hydrogen_indices)

    def compute_orientations(self):
        pass
