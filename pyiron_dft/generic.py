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
    def __init__(self, project, job_name):
        super(AtomisticGenericJob, self).__init__(project, job_name)
        self._generic_input['fix_symmetry'] = True

    @property
    def encut(self):
        return self.plane_wave_cutoff

    @encut.setter
    def encut(self, val):
        self.plane_wave_cutoff = val

    @property
    def xc(self):
        return self.exchange_correlation_functional

    @xc.setter
    def xc(self, val):
        self.exchange_correlation_functional = val

    @property
    def plane_wave_cutoff(self):
        raise NotImplementedError("The encut property is not implemented for this code.")

    @plane_wave_cutoff.setter
    def plane_wave_cutoff(self, val):
        raise NotImplementedError("The encut property is not implemented for this code.")

    @property
    def spin_constraints(self):
        raise NotImplementedError("The spin_constraints property is not implemented for this code.")

    @spin_constraints.setter
    def spin_constraints(self, val):
        raise NotImplementedError("The spin_constraints property is not implemented for this code.")

    @property
    def exchange_correlation_functional(self):
        raise NotImplementedError("The exchange property is not implemented for this code.")

    @exchange_correlation_functional.setter
    def exchange_correlation_functional(self, val):
        raise NotImplementedError("The exchange property is not implemented for this code.")

    @staticmethod
    def _get_k_mesh_by_cell(cell, kspace_per_in_ang=0.10):
        latlens = [np.linalg.norm(lat) for lat in cell]
        kmesh = np.floor(np.array([2 * np.pi / ll for ll in latlens]) / kspace_per_in_ang)
        return kmesh

    @property
    def fix_spin_constraint(self):
        return self._generic_input['fix_spin_constraint']

    @fix_spin_constraint.setter
    def fix_spin_constraint(self, boolean):
        assert isinstance(boolean, bool)
        self._generic_input['fix_spin_constraint'] = boolean

    @property
    def fix_symmetry(self):
        return self._generic_input['fix_symmetry']

    @fix_symmetry.setter
    def fix_symmetry(self, boolean):
        assert isinstance(boolean, bool)
        self._generic_input['fix_symmetry'] = boolean

    def set_kmesh_density(self, kspace_per_in_ang=0.10):
        mesh = self._get_k_mesh_by_cell(self.structure, kspace_per_in_ang)
        self.set_kpoints(mesh=mesh, scheme='MP', center_shift=None, symmetry_reduction=True, manual_kpoints=None,
                         weights=None, reciprocal=True)

    def set_kpoints(self, mesh=None, scheme='MP', center_shift=None, symmetry_reduction=True, manual_kpoints=None,
                    weights=None, reciprocal=True):
        raise NotImplementedError("The set_kpoints function is not implemented for this code.")

    def calc_static(self, electronic_steps=60, algorithm=None, retain_charge_density=False,
                    retain_electrostatic_potential=False):
        self._generic_input['fix_symmetry'] = True
        super(GenericDFTJob, self).calc_static()

    def calc_minimize(self, electronic_steps=60, ionic_steps=100, max_iter=None, pressure=None, algorithm=None,
                      retain_charge_density=False, retain_electrostatic_potential=False, ionic_energy=None,
                      ionic_forces=None, volume_only=False):
        self._generic_input['fix_symmetry'] = True
        super(GenericDFTJob, self).calc_minimize(max_iter=max_iter, pressure=pressure)

    def calc_md(self, temperature=None, n_ionic_steps=1000, n_print=1, dt=1.0, retain_charge_density=False,
                retain_electrostatic_potential=False, **kwargs):
        self._generic_input['fix_symmetry'] = False
        super(GenericDFTJob, self).calc_md(temperature=temperature, n_ionic_steps=n_ionic_steps, n_print=n_print)

    # Backward compatibility
    def get_encut(self):
        return self.encut

    def set_encut(self, encut):
        """
        Sets the plane-wave energy cutoff
        Args:
            encut (float): The energy cutoff in eV
        """
        self.plane_wave_cutoff = encut

    def set_exchange_correlation_functional(self, exchange_correlation_functional):
        self.exchange_correlation_functional = exchange_correlation_functional
