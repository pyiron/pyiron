# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron_atomistics.job.atomistic import AtomisticGenericJob

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class GenericDFTJob(AtomisticGenericJob):

    @staticmethod
    def _get_k_mesh_by_cell(cell, kspace_per_in_ang=0.10):
        latlens = [np.linalg.norm(lat) for lat in cell]
        kmesh = np.floor(np.array([2 * np.pi / ll for ll in latlens]) / kspace_per_in_ang)
        return kmesh

    def set_kmesh_density(self, kspace_per_in_ang=0.10):
        mesh = self._get_k_mesh_by_cell(self.structure, kspace_per_in_ang)
        self.set_kpoints(mesh=mesh, scheme='MP', center_shift=None, symmetry_reduction=True, manual_kpoints=None,
                         weights=None, reciprocal=True)

    def set_kpoints(self, mesh=None, scheme='MP', center_shift=None, symmetry_reduction=True, manual_kpoints=None,
                    weights=None, reciprocal=True):
        raise NotImplementedError("The set_kpoints function is not implemented for this code.")

    def set_encut(self, encut):
        raise NotImplementedError("The set_encut function is not implemented for this code.")

    def get_encut(self):
        raise NotImplementedError("The set_encut function is not implemented for this code.")

    def set_exchange_correlation_functional(self, exchange_correlation_functional):
        raise NotImplementedError("The set_exchange_correlation_functional function is not implemented for this code.")

    def calc_static(self, electronic_steps=60, algorithm=None, retain_charge_density=False,
                    retain_electrostatic_potential=False):
        super(GenericDFTJob, self).calc_static()

    def calc_minimize(self, electronic_steps=60, ionic_steps=100, max_iter=None, pressure=None, algorithm=None,
                      retain_charge_density=False, retain_electrostatic_potential=False, ionic_energy=None,
                      ionic_forces=None, volume_only=False):
        super(GenericDFTJob, self).calc_minimize(max_iter=max_iter, pressure=pressure)

    def calc_md(self, temperature=None, n_ionic_steps=1000, n_print=1, dt=1.0, retain_charge_density=False,
                retain_electrostatic_potential=False, **kwargs):
        super(GenericDFTJob, self).calc_md(temperature=temperature, n_ionic_steps=n_ionic_steps, n_print=n_print)
