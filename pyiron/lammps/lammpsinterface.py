from ctypes import c_double, c_int
from multiprocessing import Process, Pipe
import numpy as np
import pandas as pd
import warnings

try:
    from lammps import lammps
except ImportError:
    pass
from pyiron.lammps.lammps import Lammps
from pyiron.lammps.structure import UnfoldingPrism
from pyiron.base.objects.job.interactive import GenericInteractive


class LammpsInt(GenericInteractive, Lammps):
    def __init__(self, project, job_name):
        super(LammpsInt, self).__init__(project, job_name)
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
        if self.server.run_mode.interactive_non_modal:
            self._interactive_library.scatter_atoms("x", 1, 3, positions)
        else:
            self._interactive_library.scatter_atoms("x", 1, 3, (len(positions) * c_double)(*positions))
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
            print('Warning: setting upper trangular matrix might slow down the calculation')
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
        if self.server.run_mode.interactive_non_modal:
            self._interactive_library.scatter_atoms('type', 0, 1, elem_all)
        else:
            self._interactive_library.scatter_atoms('type', 0, 1, (len(elem_all) * c_int)(*elem_all))

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
                potential_lst.append([potential.split('/')[-1], potential])
        for line in self.input.potential.get_string_lst():
            if len(line) > 2:
                for potential in potential_lst:
                    if potential[0] in line:
                        line = line.replace(potential[0], potential[1])
                self._interactive_lib_command(line.split('\n')[0])

    def _reset_interactive_run_command(self):
        df = pd.DataFrame(self.input.control.dataset)
        self._interactive_run_command = " ".join(df.T[df.index[-1]].values)

    def interactive_open(self):
        if self.server.run_mode.interactive_non_modal:
            self._interactive_library = LammpsLibrary()
        else:
            self._interactive_library = lammps()
        if not all(self.structure.pbc):
            self.input.control['boundary'] = ' '.join(['p' if coord else 'f' for coord in self.structure.pbc])
        self._reset_interactive_run_command()
        self.interactive_structure_setter(self.structure)

    def calc_minimize(self, e_tol=1e-8, f_tol=1e-8, max_iter=1000, pressure=None, n_print=1):
        if self.server.run_mode.interactive_non_modal:
            warnings.warn('calc_minimize() is not implemented for the non modal interactive mode use calc_static()!')
        super(LammpsInt, self).calc_minimize(e_tol=e_tol, f_tol=f_tol, max_iter=max_iter, pressure=pressure,
                                             n_print=n_print)

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, time_step=None, n_print=100, delta_temp=1.0,
                delta_press=None, seed=None, tloop=None, rescale_velocity=True):
        if self.server.run_mode.interactive_non_modal:
            warnings.warn('calc_md() is not implemented for the non modal interactive mode use calc_static()!')
        super(LammpsInt, self).calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps,
                                       time_step=time_step, n_print=n_print, delta_temp=delta_temp,
                                       delta_press=delta_press, seed=seed, tloop=tloop,
                                       rescale_velocity=rescale_velocity)

    def run_if_interactive(self):
        if self._generic_input['calc_mode'] == 'md':
            self.input.control['run'] = self._generic_input['n_print']
            super(LammpsInt, self).run_if_interactive()
            self._reset_interactive_run_command()

            counter = 0
            iteration_max = int(self._generic_input['n_ionic_steps'] / self._generic_input['n_print'])
            while counter < iteration_max:
                self._interactive_lib_command(self._interactive_run_command)
                self.interactive_collect()
                counter += 1

        else:
            super(LammpsInt, self).run_if_interactive()
            self._interactive_lib_command(self._interactive_run_command)
            self.interactive_collect()

    def run_if_interactive_non_modal(self):
        if not self._interactive_fetch_completed:
            print('Warning: interactive_fetch being effectuated')
            self.interactive_fetch()
        super(LammpsInt, self).run_if_interactive()
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
            print('Warning: setting upper trangular matrix might slow down the calculation')
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
        for el, pos in zip(structure.get_chemical_elements(),
                           [self._interactive_prism.pos_to_lammps(position) for position in structure.positions]):
            atom_ind = el_dict[el]
            self._interactive_lib_command('create_atoms ' + str(atom_ind)
                                          + ' single ' + str(pos[0]) + ' '
                                          + str(pos[1]) + ' ' + str(pos[2]) + ' remap yes')
            # self._interactive_lib_command('change_box all remap')
        self._interactive_lammps_input()
        self._interactive_set_potential()

    def collect_output(self):
        if self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal:
            pass
        else:
            super(LammpsInt, self).collect_output()

    def update_potential(self):
        self._interactive_lib_command(self.potential.Config[0][0])
        self._interactive_lib_command(self.potential.Config[0][1])

    def interactive_indices_getter(self):
        return super(LammpsInt, self).interactive_indices_getter().tolist()

    def interactive_energy_pot_getter(self):
        return self._interactive_library.get_thermo("pe")

    def interactive_energy_tot_getter(self):
        return self._interactive_library.get_thermo("etotal")

    def interactive_steps_getter(self):
        return self._interactive_library.get_thermo("step")

    def interactive_temperatures_getter(self):
        return self._interactive_library.get_thermo("temp")

    def interactive_stress_getter(self):
        '''
        This gives back an Nx3x3 np array of stress/atom defined in
        http://lammps.sandia.gov/doc/compute_stress_atom.html
        Keep in mind that it is stress*volume in eV.
        Further discussion can be found on the website above.
        '''
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
            super(LammpsInt, self).interactive_close()


class LammpsLibrary(object):
    def __init__(self):
        lmp_interface = lammps()
        parent_conn, child_conn = Pipe()
        lammps_process = Process(target=self.interactive_run, args=(child_conn, lmp_interface))
        lammps_process.start()
        self._interactive_library = parent_conn

    def command(self, command):
        self._interactive_library.send([self.interactive_lib_command, command])

    def gather_atoms(self, *args):
        self._interactive_library.send([self.interative_gather_atoms, *args])
        return self._interactive_library.recv()

    def scatter_atoms(self, *args):
        self._interactive_library.send([self.interactive_scatter_atoms, *args])

    def get_thermo(self, *args):
        self._interactive_library.send([self.interactive_get_thermo, *args])
        return self._interactive_library.recv()

    def extract_compute(self, *args):
        self._interactive_library.send([self.interactive_extract_compute, *args])
        return self._interactive_library.recv()

    def close(self):
        self._interactive_library.send([self.interactive_close])

    @staticmethod
    def interactive_lib_command(conn, job, command):
        job.command(command)

    @staticmethod
    def interative_gather_atoms(conn, job, *args):
        return np.array(job.gather_atoms(*args))

    @staticmethod
    def interactive_scatter_atoms(conn, job, *args):
        py_vector = args[3]
        if isinstance(py_vector[0], int):
            c_vector = (len(py_vector) * c_int)(*py_vector)
        else:
            c_vector = (len(py_vector) * c_double)(*py_vector)
        job.scatter_atoms(args[0], args[1], args[2], c_vector)

    @staticmethod
    def interactive_get_thermo(conn, job, *args):
        return np.array(job.get_thermo(*args))

    @staticmethod
    def interactive_extract_compute(conn, job, *args):
        return np.array(job.extract_compute(*args))

    @staticmethod
    def interactive_close(conn, job):
        job.close()
        conn.close()
        return 'exit'

    @staticmethod
    def interactive_run(conn, job):
        while True:
            input_info = conn.recv()
            if isinstance(input_info, list):
                input_function = input_info[0]
                input_args = input_info[1:]
                answer = input_function(conn, job, *input_args)
            else:
                answer = input_info(conn, job)
            if isinstance(answer, str) and answer == 'exit':
                break
            elif answer is not None:
                conn.send(answer)


class LammpsInt2(LammpsInt):
    def __init__(self, project, job_name):
        warnings.warn('Please use LammpsInt instead of LammpsInt2')
        super(LammpsInt2, self).__init__(project=project, job_name=job_name)
