# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from ctypes import c_double, c_int
import importlib
import numpy as np
import os
import pandas as pd
import pickle
import subprocess
import warnings

from pyiron.lammps.base import LammpsBase
from pyiron.lammps.structure import UnfoldingPrism
from pyiron.atomistics.job.interactive import GenericInteractive

__author__ = "Osamu Waseda, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2018"


class LammpsInteractive(LammpsBase, GenericInteractive):
    def __init__(self, project, job_name):
        super(LammpsInteractive, self).__init__(project, job_name)
        self._check_opened = False
        self._interactive_prism = None
        self._interactive_run_command = None
        self._interactive_grand_canonical = True
        self.interactive_cache = {'cells': [],
                                  'energy_pot': [],
                                  'energy_tot': [],
                                  'forces': [],
                                  'positions': [],
                                  'pressures': [],
                                  'steps': [],
                                  'indices': [],
                                  'temperature': [],
                                  'computation_time': [],
                                  'volume': []}

    @property
    def structure(self):
        return GenericInteractive.structure.fget(self)

    @structure.setter
    def structure(self, structure):
        GenericInteractive.structure.fset(self, structure)

    def get_structure(self, iteration_step=-1):
        return GenericInteractive.get_structure(self, iteration_step=iteration_step)

    def _interactive_lib_command(self, command):
        self._logger.debug('Lammps library: ' + command)
        self._interactive_library.command(command)

    def interactive_positions_getter(self):
        positions = np.reshape(np.array(self._interactive_library.gather_atoms("x", 1, 3)),
                               (len(self.structure), 3))
        if np.matrix.trace(self._interactive_prism.R) != 3:
            positions = np.dot(positions, self._interactive_prism.R.T)
        return positions.tolist()

    def interactive_positions_setter(self, positions):
        if np.matrix.trace(self._interactive_prism.R) != 3:
            positions = np.array(positions).reshape(-1, 3)
            positions = np.dot(positions, self._interactive_prism.R)
        positions = np.array(positions).flatten()
        if self.server.run_mode.interactive and self.server.cores == 1:
            self._interactive_library.scatter_atoms("x", 1, 3, (len(positions) * c_double)(*positions))
        else:
            self._interactive_library.scatter_atoms("x", 1, 3, positions)
        self._interactive_lib_command('change_box all remap')

    def interactive_cells_getter(self):
        cc = np.array([[self._interactive_library.get_thermo('lx'),
                        0,
                        0],
                       [self._interactive_library.get_thermo('xy'),
                        self._interactive_library.get_thermo('ly'),
                        0],
                       [self._interactive_library.get_thermo('xz'),
                        self._interactive_library.get_thermo('yz'),
                        self._interactive_library.get_thermo('lz')]])
        return self._interactive_prism.unfold_cell(cc)

    def interactive_cells_setter(self, cell):
        self._interactive_prism = UnfoldingPrism(cell)
        lx, ly, lz, xy, xz, yz = self._interactive_prism.get_lammps_prism()
        if np.matrix.trace(self._interactive_prism.R) != 3:
            warnings.warn('Warning: setting upper trangular matrix might slow down the calculation')
        if abs(xy) + abs(xz) + abs(yz) > 1.0e-6:
            if self.structure._is_scaled:
                self._interactive_lib_command(
                    'change_box all x final 0 %f y final 0 %f z final 0 %f \
                     xy final %f xz final %f yz final %f triclinic remap units box' % (lx, ly, lz, xy, xz, yz))
            else:
                self._interactive_lib_command(
                    'change_box all x final 0 %f y final 0 %f z final 0 %f \
                    xy final %f xz final %f yz final %f triclinic units box' % (lx, ly, lz, xy, xz, yz))
        else:
            if self.structure._is_scaled:
                self._interactive_lib_command(
                    'change_box all x final 0 %f y final 0 %f z final 0 %f remap units box' % (lx, ly, lz))
            else:
                self._interactive_lib_command(
                    'change_box all x final 0 %f y final 0 %f z final 0 %f units box' % (lx, ly, lz))

    def interactive_indices_setter(self, indices):
        el_struct_lst = self._structure_current.get_species_symbols()
        el_obj_lst = self._structure_current.get_species_objects()
        el_eam_lst = self.input.potential.get_element_lst()
        el_dict = {}
        for id_eam, el_eam in enumerate(el_eam_lst):
            if el_eam in el_struct_lst:
                id_el = list(el_struct_lst).index(el_eam)
                el = el_obj_lst[id_el]
                el_dict[el] = id_eam + 1
        elem_all = np.array([el_dict[self._structure_current.species[el]] for el in indices])
        if self.server.run_mode.interactive and self.server.cores == 1:
            self._interactive_library.scatter_atoms('type', 0, 1, (len(elem_all) * c_int)(*elem_all))
        else:
            self._interactive_library.scatter_atoms('type', 0, 1, elem_all)

    def interactive_volume_getter(self):
        return self._interactive_library.get_thermo('vol')

    def interactive_forces_getter(self):
        ff = np.reshape(np.array(self._interactive_library.gather_atoms("f", 1, 3)), (len(self.structure), 3))
        if np.matrix.trace(self._interactive_prism.R) != 3:
            ff = np.dot(ff, self._interactive_prism.R.T)
        return ff.tolist()

    def _interactive_lammps_input(self):
        del self.input.control['dump']
        del self.input.control['dump_modify']
        for key, value in zip(self.input.control.dataset['Parameter'], self.input.control.dataset['Value']):
            if key in ['read_data', 'units', 'dimension', 'boundary', 'atom_style', 'atom_modify', 'include', 'run',
                       'minimize']:
                continue
            else:
                self._interactive_lib_command(key + ' ' + str(value))

    def _interactive_set_potential(self):
        potential_lst = []
        if self.input.potential.files is not None:
            for potential in self.input.potential.files:
                if not os.path.exists(potential):
                    raise ValueError('Potential not found: ', potential)
                potential_lst.append([potential.split('/')[-1], potential])
        for line in self.input.potential.get_string_lst():
            if len(line) > 2:
                for potential in potential_lst:
                    if potential[0] in line:
                        line = line.replace(potential[0], potential[1])
                self._interactive_lib_command(line.split('\n')[0])

    def _executable_activate_mpi(self):
        if self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal:
            pass
        else:
            super(LammpsInteractive, self)._executable_activate_mpi()

    def _reset_interactive_run_command(self):
        df = pd.DataFrame(self.input.control.dataset)
        self._interactive_run_command = " ".join(df.T[df.index[-1]].values)

    def interactive_initialize_interface(self):
        if self.server.run_mode.interactive and self.server.cores == 1:
            lammps = getattr(importlib.import_module('lammps'), 'lammps')
            self._interactive_library = lammps()
        else:
            self._create_working_directory()
            self._interactive_library = LammpsLibrary(cores=self.server.cores, working_directory=self.working_directory)
        if not all(self.structure.pbc):
            self.input.control['boundary'] = ' '.join(['p' if coord else 'f' for coord in self.structure.pbc])
        self._reset_interactive_run_command()
        self.interactive_structure_setter(self.structure)

    def calc_minimize(self, e_tol=1e-8, f_tol=1e-8, max_iter=1000, pressure=None, n_print=100):
        if self.server.run_mode.interactive_non_modal:
            warnings.warn('calc_minimize() is not implemented for the non modal interactive mode use calc_static()!')
        super(LammpsInteractive, self).calc_minimize(e_tol=e_tol, f_tol=f_tol, max_iter=max_iter, pressure=pressure,
                                                     n_print=n_print)
        if self.interactive_is_activated() and \
                (self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal):
            self.interactive_structure_setter(self.structure)

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, time_step=1.0, n_print=100,
                temperature_damping_timescale=100.0, pressure_damping_timescale=1000.0, seed=None, tloop=None,
                initial_temperature=None, langevin=False, delta_temp=None, delta_press=None):
        super(LammpsInteractive, self).calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps,
                                               time_step=time_step, n_print=n_print,
                                               temperature_damping_timescale=temperature_damping_timescale,
                                               pressure_damping_timescale=pressure_damping_timescale, seed=seed,
                                               tloop=tloop, initial_temperature=initial_temperature,
                                               langevin=langevin, delta_temp=delta_temp, delta_press=delta_press)
        if self.interactive_is_activated() and \
                (self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal):
            self.interactive_structure_setter(self.structure)

    def run_if_interactive(self):
        if self._generic_input['calc_mode'] == 'md':
            self.input.control['run'] = self._generic_input['n_print']
            super(LammpsInteractive, self).run_if_interactive()
            self._reset_interactive_run_command()

            counter = 0
            iteration_max = int(self._generic_input['n_ionic_steps'] / self._generic_input['n_print'])
            while counter < iteration_max:
                self._interactive_lib_command(self._interactive_run_command)
                self.interactive_collect()
                counter += 1

        else:
            super(LammpsInteractive, self).run_if_interactive()
            self._interactive_lib_command(self._interactive_run_command)
            self.interactive_collect()

    def run_if_interactive_non_modal(self):
        if not self._interactive_fetch_completed:
            print('Warning: interactive_fetch being effectuated')
            self.interactive_fetch()
        super(LammpsInteractive, self).run_if_interactive()
        self._interactive_lib_command(self._interactive_run_command)
        self._interactive_fetch_completed = False

    def interactive_fetch(self):
        if self._interactive_fetch_completed and self.server.run_mode.interactive_non_modal:
            print('First run and then fetch')
        else:
            self.interactive_collect()
            self._logger.debug('interactive run - done')

    def interactive_structure_setter(self, structure):
        self._interactive_lib_command('clear')
        self._set_selective_dynamics()
        self._interactive_lib_command('units ' + self.input.control['units'])
        self._interactive_lib_command('dimension ' + str(self.input.control['dimension']))
        self._interactive_lib_command('boundary ' + self.input.control['boundary'])
        self._interactive_lib_command('atom_style ' + self.input.control['atom_style'])
        self._interactive_lib_command("atom_modify map array")
        self._interactive_prism = UnfoldingPrism(structure.cell)
        if np.matrix.trace(self._interactive_prism.R) != 3:
            warnings.warn('Warning: setting upper trangular matrix might slow down the calculation')
        xhi, yhi, zhi, xy, xz, yz = self._interactive_prism.get_lammps_prism()
        if self._interactive_prism.is_skewed():
            self._interactive_lib_command('region 1 prism' +
                                          ' 0.0 ' + str(xhi) + ' 0.0 ' + str(yhi) + ' 0.0 ' + str(zhi) +
                                          ' ' + str(xy) + ' ' + str(xz) + ' ' + str(yz) + ' units box')
        else:
            self._interactive_lib_command('region 1 block' +
                                          ' 0.0 ' + str(xhi) + ' 0.0 ' + str(yhi) + ' 0.0 ' + str(zhi) + ' units box')
        el_struct_lst = self.structure.get_species_symbols()
        el_obj_lst = self.structure.get_species_objects()
        el_eam_lst = self.input.potential.get_element_lst()
        self._interactive_lib_command('create_box ' + str(len(el_eam_lst)) + ' 1')
        el_dict = {}
        for id_eam, el_eam in enumerate(el_eam_lst):
            if el_eam in el_struct_lst:
                id_el = list(el_struct_lst).index(el_eam)
                el = el_obj_lst[id_el]
                el_dict[el] = id_eam + 1
                self._interactive_lib_command('mass {0:3d} {1:f}'.format(id_eam + 1, el.AtomicMass))
            else:
                self._interactive_lib_command('mass {0:3d} {1:f}'.format(id_eam + 1, 1.00))
        self._interactive_lib_command('create_atoms 1 random ' + str(len(structure)) + ' 12345 1')
        positions = structure.positions.flatten()
        elem_all = np.array([el_dict[el] for el in structure.get_chemical_elements()])
        if self.server.run_mode.interactive and self.server.cores == 1:
            self._interactive_library.scatter_atoms("x", 1, 3, (len(positions) * c_double)(*positions))
            self._interactive_library.scatter_atoms('type', 0, 1, (len(elem_all) * c_int)(*elem_all))
        else:
            self._interactive_library.scatter_atoms("x", 1, 3, positions)
            self._interactive_library.scatter_atoms('type', 0, 1, elem_all)
        self._interactive_lib_command('change_box all remap')
        self._interactive_lammps_input()
        self._interactive_set_potential()

    def from_hdf(self, hdf=None, group_name=None):
        """
        Recreates instance from the hdf5 file

        Args:
            hdf (str): Path to the hdf5 file
            group_name (str): Name of the group which contains the object
        """
        super(LammpsInteractive, self).from_hdf(hdf=hdf, group_name=group_name)
        self.species_from_hdf()

    def collect_output(self):
        if self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal:
            pass
        else:
            super(LammpsInteractive, self).collect_output()

    def update_potential(self):
        self._interactive_lib_command(self.potential.Config[0][0])
        self._interactive_lib_command(self.potential.Config[0][1])

    def interactive_indices_getter(self):
        return super(LammpsInteractive, self).interactive_indices_getter().tolist()

    def interactive_energy_pot_getter(self):
        return self._interactive_library.get_thermo("pe")

    def interactive_energy_tot_getter(self):
        return self._interactive_library.get_thermo("etotal")

    def interactive_steps_getter(self):
        return self._interactive_library.get_thermo("step")

    def interactive_temperatures_getter(self):
        return self._interactive_library.get_thermo("temp")

    def interactive_stress_getter(self):
        """
        This gives back an Nx3x3 array of stress/atom defined in http://lammps.sandia.gov/doc/compute_stress_atom.html
        Keep in mind that it is stress*volume in eV. Further discussion can be found on the website above.

        Returns:
            numpy.array: Nx3x3 np array of stress/atom
        """
        if not 'stress' in self.interactive_cache.keys():
            self._interactive_lib_command('compute st all stress/atom NULL')
            self._interactive_lib_command('run 0')
            self.interactive_cache['stress'] = []
        ss = np.array([self._interactive_library.extract_compute('st', 1, 2)[i][j + (j != k) * (k + 2)]
                       for i in range(len(self.structure))
                       for j in range(3)
                       for k in range(3)]).reshape(len(self.structure), 3, 3)/1.602e6
        if np.matrix.trace(self._interactive_prism.R) != 3:
            ss = np.dot(np.dot(self._interactive_prism.R, ss), self._interactive_prism.R.T)
        return ss

    def interactive_pressures_getter(self):
        pp = np.array([[self._interactive_library.get_thermo('pxx'),
                        self._interactive_library.get_thermo('pxy'),
                        self._interactive_library.get_thermo('pxz')],
                       [self._interactive_library.get_thermo('pxy'),
                        self._interactive_library.get_thermo('pyy'),
                        self._interactive_library.get_thermo('pyz')],
                       [self._interactive_library.get_thermo('pxz'),
                        self._interactive_library.get_thermo('pyz'),
                        self._interactive_library.get_thermo('pzz')]])
        if np.matrix.trace(self._interactive_prism.R) != 3:
            pp = np.dot(np.dot(self._interactive_prism.R, pp), self._interactive_prism.R.T)
        return pp / 10000  # bar -> GPa

    def interactive_close(self):
        if self.interactive_is_activated():
            self._interactive_library.close()
            with self.project_hdf5.open("output") as h5:
                if 'interactive' in h5.list_groups():
                    for key in h5['interactive'].list_nodes():
                        h5['generic/' + key] = h5['interactive/' + key]
            super(LammpsInteractive, self).interactive_close()


class LammpsLibrary(object):
    def __init__(self, cores=1, working_directory='.'):
        executable = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sub', 'lmpmpi.py')
        # print(executable)
        self._process = subprocess.Popen(['mpiexec', '-n', str(cores), 'python', executable],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         stdin=subprocess.PIPE,
                                         cwd=working_directory)

    def _send(self, command, data=None):
        """
        Send a command to the Lammps Library executable

        Args:
            command (str): command to be send to the
            data:
        """
        # print('send: ', {'c': command, 'd': data})
        pickle.dump({'c': command, 'd': data}, self._process.stdin)
        self._process.stdin.flush()

    def _receive(self):
        """
        Receive data from the Lammps library

        Returns:
            data
        """
        output = pickle.load(self._process.stdout)
        # print(output)
        return output

    def command(self, command):
        """
        Send a command to the lammps library

        Args:
            command (str):
        """
        self._send(command='command', data=command)

    def gather_atoms(self, *args):
        """
        Gather atoms from the lammps library

        Args:
            *args:

        Returns:
            np.array
        """
        self._send(command='gather_atoms', data=list(args))
        return self._receive()

    def scatter_atoms(self, *args):
        """
        Scatter atoms for the lammps library

        Args:
            *args:
        """
        self._send(command='scatter_atoms', data=list(args))

    def get_thermo(self, *args):
        """
        Get thermo from the lammps library

        Args:
            *args:

        Returns:

        """
        self._send(command='get_thermo', data=list(args))
        return self._receive()

    def extract_compute(self, *args):
        """
        Extract compute from the lammps library

        Args:
            *args:

        Returns:

        """
        self._send(command='extract_compute', data=list(args))
        return self._receive()

    def close(self):
        self._send(command='close')
        self._process.kill()

    def __del__(self):
        # print('object killed __del__')
        self.close()
