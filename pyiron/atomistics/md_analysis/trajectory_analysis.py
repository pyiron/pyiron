# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, unicode_literals
import numpy as np
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.base.job.path import JobPath
from pyiron.atomistics.structure.atoms import Atoms

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class TrajectoryAnalysis(object):

    def __init__(self):
        self._positions = None
        self.cells = None
        self.forces = None
        self.total_energies = None
        self.potential_energies = None
        self.pressures = None
        self.temperatures = None
        self._job = None
        self._trajectory_path = "output/generic"
        self.steps = None

    @property
    def job(self):
        return self._job

    @job.setter
    def job(self, val):
        if not (isinstance(val, (AtomisticGenericJob, JobPath))):
            raise AssertionError()
        if not (val.status == "finished"):
            raise AssertionError()
        self._job = val
        self._set_trajectory()

    @property
    def structure(self):
        if isinstance(self.job, JobPath):
            return Atoms().from_hdf(self.job, "output/structure")
        else:
            return self.job.structure

    def get_quantity(self, quantity):
        path = "/".join([self._trajectory_path, str(quantity)])
        return self.job[path]

    def _set_trajectory(self):
        self.positions = self.get_quantity("positions")
        self.forces = self.get_quantity("forces")
        self.cells = self.get_quantity("cells")
        self.temperatures = self.get_quantity("temperatures")
        self.total_energies = self.get_quantity("energy_tot")
        # self.potential_energies = self.get_quantity("energy_pot")
        self.pressure = self.get_quantity("pressures")


def unwrap_coordinates(positions, cell=None, is_relative=False):
    unwrapped_positions = positions.copy()
    if len(unwrapped_positions) > 1:
        if not is_relative:
            if not (cell is not None):
                raise AssertionError()
            rel_positions = np.dot(unwrapped_positions, np.linalg.inv(cell))
        else:
            rel_positions = unwrapped_positions.copy()
        d_pos = np.round(rel_positions - np.roll(rel_positions, 1, axis=0), 0)
        d_pos = np.cumsum(d_pos, axis=0)
        pos = rel_positions - d_pos
        return pos
    else:
        return unwrapped_positions
