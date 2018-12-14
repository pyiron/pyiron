# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.base.job.interactive import InteractiveBase
from pyiron.atomistics.job.interface import AtomisticInteractiveInterface
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.atomistics.job.atomistic import AtomisticGenericJob, GenericOutput

__author__ = "Osamu Waseda, Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class GenericInteractive(AtomisticGenericJob, InteractiveBase):
    def __init__(self, project, job_name):
        super(GenericInteractive, self).__init__(project, job_name)
        self._interface = AtomisticInteractiveInterface(job=self)
        self.output = GenericInteractiveOutput(job=self)

    @property
    def interactive_enforce_structure_reset(self):
        return self._interface.interactive_enforce_structure_reset

    @interactive_enforce_structure_reset.setter
    def interactive_enforce_structure_reset(self, reset):
        self._interface.interactive_enforce_structure_reset = reset

    @property
    def initial_structure(self):
        return AtomisticGenericJob.structure.fget(self)

    @property
    def current_structure(self):
        return self.structure

    @current_structure.setter
    def current_structure(self, structure):
        self.structure = structure

    @property
    def structure(self):
        if self._interface.structure_current is not None:
            return self._interface.structure_current
        elif self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal:
            self._interface.structure_current = AtomisticGenericJob.structure.fget(self)
            return self._interface.structure_current
        else:
            return AtomisticGenericJob.structure.fget(self)

    @structure.setter
    def structure(self, structure):
        if self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal:
            # only overwrite the initial structure if it is not set already.
            if AtomisticGenericJob.structure.fget(self) is None:
                AtomisticGenericJob.structure.fset(self, structure.copy())
            self._interface.structure_current = structure
        else:
            AtomisticGenericJob.structure.fset(self, structure)

    def get_structure(self, iteration_step=-1):
        if (self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal) \
                and self.interactive_is_activated():
            # Warning: We only copy symbols, positions and cell information - no tags.
            if len(self.output.indices) != 0:
                el_lst = [el.Abbreviation for el in self.structure.species]
                return Atoms(symbols=np.array([el_lst[el] for el in self.output.indices[iteration_step]]),
                             positions=self.output.positions[iteration_step],
                             cell=self.output.cells[iteration_step])
            else:
                return None
        else:
            if self.get("output/generic/cells") is not None and len(self.get("output/generic/cells")) != 0:
                return super(GenericInteractive, self).get_structure(iteration_step=iteration_step)
            else:
                return None

    def run_if_interactive(self):
        self._interface.run_if_interactive(job=self)

    def interactive_collect(self):
        self._interface.interactive_collect(job=self)

    def run_if_interactive_non_modal(self):
        self._interface.run_if_interactive_non_modal(job=self)

    def interactive_fetch(self):
        self._interface.interactive_fetch(job=self)

    def interactive_close(self):
        self._interface.interactive_close(job=self)

    # Functions which have to be implemented by the fin
    def interactive_cells_getter(self):
        return self.initial_structure.cell

    def interactive_indices_getter(self):
        return self.current_structure.get_chemical_indices()

    def interactive_positions_getter(self):
        return self.current_structure.positions

    def interactive_steps_getter(self):
        return self._interface.interactive_steps_getter()

    def interactive_time_getter(self):
        return self.interactive_steps_getter()

    def interactive_volume_getter(self):
        return self.initial_structure.get_volume()

    def interactive_cells_setter(self, cell):
        self._interface.interactive_cells_setter(job=self, cell=cell)

    def interactive_energy_pot_getter(self):
        return self._interface.interactive_energy_pot_getter(job=self)

    def interactive_energy_tot_getter(self):
        return self._interface.interactive_energy_tot_getter(job=self)

    def interactive_forces_getter(self):
        return self._interface.interactive_forces_getter(job=self)

    def interactive_indices_setter(self, indices):
        self._interface.interactive_indices_setter(job=self, indices=indices)

    def interactive_spins_getter(self):
        return self._interface.interactive_spins_getter(job=self)

    def interactive_spin_constraints_getter(self):
        return self._interface.interactive_spin_constraints_getter(job=self)

    def interactive_magnetic_forces_getter(self):
        return self._interface.interactive_magnetic_forces_getter(job=self)

    def interactive_spin_constraints_setter(self, spins):
        self._interface.interactive_spin_constraints_setter(job=self, spins=spins)

    def interactive_open(self):
        self._interface.interactive_open(job=self)

    def interactive_positions_setter(self, positions):
        self._interface.interactive_positions_setter(job=self, positions=positions)

    def interactive_pressures_getter(self):
        return self._interface.interactive_pressures_getter(job=self)

    def interactive_stress_getter(self):
        return self._interface.interactive_stress_getter(job=self)

    def interactive_structure_setter(self, structure):
        self._interface.interactive_structure_setter(job=self, structure=structure)

    def interactive_temperatures_getter(self):
        return self._interface.interactive_temperatures_getter(job=self)

    def interactive_unwrapped_positions_getter(self):
        return self._interface.interactive_unwrapped_positions_getter(job=self)


class GenericInteractiveOutput(GenericOutput):
    def __init__(self, job):
        self._job = job

    def _key_from_cache(self, key):
        if key in self._job._interface.interactive_cache.keys() and self._job._interface.interactive_is_activated() \
                and len(self._job._interface.interactive_cache[key]) != 0:
            return self._job._interface.interactive_cache[key]
        else:
            return []

    def _lst_from_cache(self, key):
        lst = self._key_from_cache(key)
        if len(lst) != 0 and isinstance(lst[-1], list):
            return [np.array(out) for out in lst]
        else:
            return lst

    def _key_from_hdf(self, key):
        return self._job['output/interactive/' + key]

    def _key_from_property(self, key, prop):
        return_lst = self._key_from_cache(key)
        hdf5_output = self._key_from_hdf(key)
        if hdf5_output is not None:
            return_lst = hdf5_output.tolist() + return_lst
        else:
            prop_result = prop(self)
            if prop_result is not None:
                return_lst = prop(self).tolist() + return_lst
        return np.array(return_lst)

    def _lst_from_property(self, key, prop=None):
        return_lst = self._lst_from_cache(key)
        hdf5_output = self._key_from_hdf(key)
        if hdf5_output is not None and len(hdf5_output) != 0:
            if isinstance(hdf5_output[-1], list):
                return_lst = [np.array(out) for out in hdf5_output] + return_lst
            else:
                return_lst = hdf5_output.tolist() + return_lst
        elif prop is not None:
            prop_result = prop(self)
            if prop_result is not None:
                return_lst = prop_result.tolist() + return_lst
        return np.array(return_lst)

    @property
    def indices(self):
        return self._lst_from_property(key='indices')

    @property
    def cells(self):
        return self._key_from_property(key='cells', prop=GenericOutput.cells.fget)

    @property
    def energy_pot(self):
        return self._key_from_property(key='energy_pot', prop=GenericOutput.energy_pot.fget)

    @property
    def energy_tot(self):
        return self._key_from_property(key='energy_tot', prop=GenericOutput.energy_tot.fget)

    @property
    def forces(self):
        return self._lst_from_property(key='forces', prop=GenericOutput.forces.fget)

    @property
    def positions(self):
        return self._lst_from_property(key='positions', prop=GenericOutput.positions.fget)

    @property
    def pressures(self):
        return self._key_from_property(key='pressures', prop=GenericOutput.pressures.fget)

    @property
    def steps(self):
        return self._key_from_property(key='steps', prop=GenericOutput.steps.fget)

    @property
    def temperature(self):
        return self._key_from_property(key='temperature', prop=GenericOutput.temperature.fget)

    @property
    def time(self):
        return self._key_from_property(key='computation_time', prop=GenericOutput.computation_time.fget)

    @property
    def unwrapped_positions(self):
        return self._lst_from_property(key='unwrapped_positions', prop=GenericOutput.unwrapped_positions.fget)

    @property
    def volume(self):
        return self._key_from_property(key='volume', prop=GenericOutput.volume.fget)

    def __dir__(self):
        return list(set(list(self._job._interface.interactive_cache.keys()) + super(GenericOutput).__dir__()))


class InteractiveInterface(object):

    def get_cell(self):
        raise NotImplementedError

    def set_cell(self, cell):
        raise NotImplementedError

    def get_temperature(self):
        raise NotImplementedError

    def set_temperature(self, temperature):
        raise NotImplementedError

    def get_positions(self):
        raise NotImplementedError

    def set_positions(self, positions):
        raise NotImplementedError

    def get_forces(self):
        raise NotImplementedError

    def set_forces(self, forces):
        raise NotImplementedError

    def get_energy_tot(self):
        raise NotImplementedError

    def set_energy_tot(self, energy_tot):
        raise NotImplementedError

    def get_energy_pot(self):
        raise NotImplementedError

    def set_energy_pot(self, energy_pot):
        raise NotImplementedError

    def get_pressure(self):
        raise NotImplementedError

    def set_pressure(self, pressure):
        raise NotImplementedError

    def run_interactive(self):
        raise NotImplementedError
