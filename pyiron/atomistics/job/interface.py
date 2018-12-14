import numpy as np
from pyiron.base.job.interface import InteractiveInterface


class AtomisticInteractiveInterface(InteractiveInterface):
    def __init__(self):
        super(AtomisticInteractiveInterface, self).__init__()
        self._interactive_enforce_structure_reset = False
        self._interactive_grand_canonical = False
        self._interactive_fetch_completed = True
        self.structure_previous = None
        self.structure_current = None
        self.interactive_cache = {'cells': [],
                                  'energy_pot': [],
                                  'energy_tot': [],
                                  'forces': [],
                                  'positions': [],
                                  'pressures': [],
                                  'stress': [],
                                  'steps': [],
                                  'temperature': [],
                                  'indices': [],
                                  'computation_time': [],
                                  'unwrapped_positions': [],
                                  'atom_spin_constraints': [],
                                  'atom_spins': [],
                                  'magnetic_forces': [],
                                  'volume': []}

    @property
    def interactive_enforce_structure_reset(self):
        return self._interactive_enforce_structure_reset

    @interactive_enforce_structure_reset.setter
    def interactive_enforce_structure_reset(self, reset):
        if not isinstance(reset, bool):
            raise AssertionError()
        self._interactive_enforce_structure_reset = reset

    def run_if_interactive(self, job):
        job.status.running = True
        if job.structure is None:
            raise ValueError("Input structure not set. Use method set_structure()")
        if not self.interactive_is_activated():
            self.interactive_open(job=job)
        if self.structure_previous is None:
            pre_struct = job.get_structure(-1)
            if pre_struct is not None:
                self.structure_previous = pre_struct
            else:
                self.structure_previous = job.structure.copy()
        if self.structure_current is not None:
            if len(self.structure_current) != len(self.structure_previous) and not self._interactive_grand_canonical:
                raise ValueError('The number of atoms changed, this is currently not supported!')
            el_lst = list(set(list(self.structure_current.get_species_symbols()) +
                              list(self.structure_previous.get_species_symbols())))
            current_structure_index = [el_lst.index(el) for el in self.structure_current.get_chemical_symbols()]
            previous_structure_index = [el_lst.index(el) for el in self.structure_previous.get_chemical_symbols()]
            if np.array_equal(np.array(current_structure_index), np.array(previous_structure_index)) and \
                    not self._interactive_enforce_structure_reset:
                if not np.allclose(self.structure_current.cell, self.structure_previous.cell, rtol=1e-15, atol=1e-15):
                    job._logger.debug('Generic library: cell changed!')
                    self.interactive_cells_setter(job=job, cell=self.structure_current.cell)
                if not np.allclose(self.structure_current.scaled_positions,
                                   self.structure_previous.scaled_positions, rtol=1e-15, atol=1e-15):
                    job._logger.debug('Generic library: positions changed!')
                    self.interactive_positions_setter(job=job, positions=self.structure_current.positions)
                if np.any(self.structure_current.get_initial_magnetic_moments()) and \
                        not np.allclose(self.structure_current.get_initial_magnetic_moments(),
                                        self.structure_previous.get_initial_magnetic_moments()):
                    job._logger.debug('Generic library: magnetic moments changed!')
                    self.interactive_spin_constraints_setter(job=job,
                                                             spins=self.structure_current.get_initial_magnetic_moments())
            elif not self._interactive_enforce_structure_reset and \
                    len(self.structure_current) == len(self.structure_previous) and \
                    np.allclose(self.structure_current.cell, self.structure_previous.cell, rtol=1e-15, atol=1e-15) and \
                    np.allclose(self.structure_current.scaled_positions,
                                self.structure_previous.scaled_positions, rtol=1e-15, atol=1e-15) and \
                    (not np.any(self.structure_current.get_initial_magnetic_moments()) or
                     np.allclose(self.structure_current.get_initial_magnetic_moments(),
                                 self.structure_previous.get_initial_magnetic_moments())):
                job._logger.debug('Generic library: indices changed!')
                self.interactive_indices_setter(job=job, indices=self.structure_current.indices)
            else:
                job._logger.debug('Generic library: structure changed!')
                self.interactive_structure_setter(job=job, structure=self.structure_current)
            self.structure_previous = self.structure_current.copy()

    def run_if_interactive_non_modal(self, job):
        raise NotImplementedError('run_if_interactive_non_modal() is not implemented!')

    def interactive_cells_getter(self, job):
        return job.initial_structure.cell

    def interactive_collect(self, job):
        if 'cells' in self.interactive_cache.keys():
            self.interactive_cache['cells'].append(self.interactive_cells_getter(job=job))
        if 'energy_pot' in self.interactive_cache.keys():
            self.interactive_cache['energy_pot'].append(self.interactive_energy_pot_getter(job=job))
        if 'energy_tot' in self.interactive_cache.keys():
            self.interactive_cache['energy_tot'].append(self.interactive_energy_tot_getter(job=job))
        if 'forces' in self.interactive_cache.keys():
            self.interactive_cache['forces'].append(self.interactive_forces_getter(job=job))
        if 'positions' in self.interactive_cache.keys():
            self.interactive_cache['positions'].append(self.interactive_positions_getter(job=job))
        if 'pressures' in self.interactive_cache.keys():
            self.interactive_cache['pressures'].append(self.interactive_pressures_getter(job=job))
        if 'stress' in self.interactive_cache.keys():
            self.interactive_cache['stress'].append(self.interactive_stress_getter(job=job))
        if 'steps' in self.interactive_cache.keys():
            self.interactive_cache['steps'].append(self.interactive_steps_getter(job=job))
        if 'temperature' in self.interactive_cache.keys():
            self.interactive_cache['temperature'].append(self.interactive_temperatures_getter(job=job))
        if 'computation_time' in self.interactive_cache.keys():
            self.interactive_cache['computation_time'].append(self.interactive_time_getter(job=job))
        if 'indices' in self.interactive_cache.keys():
            self.interactive_cache['indices'].append(self.interactive_indices_getter(job=job))
        if 'atom_spins' in self.interactive_cache.keys():
            self.interactive_cache['atom_spins'].append(self.interactive_spins_getter(job=job))
        if 'atom_spin_constraints' in self.interactive_cache.keys():
            if job._generic_input['fix_spin_constraint']:
                self.interactive_cache['atom_spin_constraints'].append(self.interactive_spin_constraints_getter(job=job))
        if 'magnetic_forces' in self.interactive_cache.keys():
            if job._generic_input['fix_spin_constraint']:
                self.interactive_cache['magnetic_forces'].append(self.interactive_magnetic_forces_getter(job=job))
        if 'unwrapped_positions' in self.interactive_cache.keys():
            self.interactive_cache['unwrapped_positions'].append(self.interactive_unwrapped_positions_getter(job=job))
        if 'volume' in self.interactive_cache.keys():
            self.interactive_cache['volume'].append(self.interactive_volume_getter(job=job))
        if len(list(self.interactive_cache.keys())) > 0 and \
                len(self.interactive_cache[list(self.interactive_cache.keys())[0]]) \
                % self._interactive_flush_frequency == 0:
            self.interactive_flush(job=job, path="interactive")
        if job.server.run_mode.interactive_non_modal:
            self._interactive_fetch_completed = True

    def interactive_indices_getter(self, job):
        return job.current_structure.get_chemical_indices()

    def interactive_positions_getter(self, job):
        return job.current_structure.positions

    def interactive_steps_getter(self, job):
        return len(self.interactive_cache[list(self.interactive_cache.keys())[0]])

    def interactive_time_getter(self, job):
        return self.interactive_steps_getter(job=job)

    def interactive_volume_getter(self, job):
        return job.initial_structure.get_volume()

    # Functions which have to be implemented by the fin
    def interactive_cells_setter(self, job, cell):
        raise NotImplementedError('interactive_cells_getter() is not implemented!')

    def interactive_energy_pot_getter(self, job):
        raise NotImplementedError('interactive_energy_pot_getter() is not implemented!')

    def interactive_energy_tot_getter(self, job):
        raise NotImplementedError('interactive_energy_tot_getter() is not implemented!')

    def interactive_forces_getter(self, job):
        raise NotImplementedError('interactive_forces_getter() is not implemented!')

    def interactive_indices_setter(self, job, indices):
        raise NotImplementedError('interactive_indices_setter() is not implemented!')

    def interactive_spins_getter(self, job):
        raise NotImplementedError('interactive_spins_getter() is not implemented!')

    def interactive_spin_constraints_getter(self, job):
        raise NotImplementedError('interactive_spin_constraints_getter() is not implemented!')

    def interactive_magnetic_forces_getter(self, job):
        raise NotImplementedError('interactive_magnetic_forces_getter() is not implemented!')

    def interactive_spin_constraints_setter(self, job, spins):
        raise NotImplementedError('iinteractive_spin_constraints_setter() is not implemented!')

    def interactive_open(self, job):
        raise NotImplementedError('interactive_open() is not implemented!')

    def interactive_fetch(self, job):
        raise NotImplementedError('interactive_fetch() is not implemented!')

    def interactive_positions_setter(self, job, positions):
        raise NotImplementedError('interactive_positions_setter() is not implemented!')

    def interactive_pressures_getter(self, job):
        raise NotImplementedError('interactive_pressures_getter() is not implemented!')

    def interactive_stress_getter(self, job):
        raise NotImplementedError('interactive_stress_getter() is not implemented!')

    def interactive_structure_setter(self, job, structure):
        raise NotImplementedError('interactive_structure_setter() is not implemented!')

    def interactive_temperatures_getter(self, job):
        raise NotImplementedError('interactive_temperatures_getter() is not implemented!')

    def interactive_unwrapped_positions_getter(self, job):
        raise NotImplementedError('interactive_unwrapped_positions_getter() is not implemented!')
