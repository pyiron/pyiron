import ast
from ctypes import c_double, c_int
import h5py
import numpy as np
import os
import pandas as pd

try:
    from lammps import lammps
except ImportError:
    pass

from pyiron.lammps.structure import LammpsStructure, UnfoldingPrism
from pyiron.lammps.pipe import LammpsLibrary
from pyiron.base.job.interface import FileInterface
from pyiron.atomistics.job.interface import AtomisticInteractiveInterface
from pyiron.base.pyio.parser import Logstatus, extract_data_from_file
from pyiron.atomistics.md_analysis.trajectory_analysis import unwrap_coordinates


class LammpsInterface(FileInterface):
    def write_input(self, job):
        """
        Call routines that generate the code specific input files

        Returns:

        """
        lmp_structure = self._get_lammps_structure(job=job, structure=job.structure, cutoff_radius=job.cutoff_radius)
        lmp_structure.write_file(file_name="structure.inp", cwd=job.working_directory)
        if job.executable.version and 'dump_modify' in job.input.control._dataset['Parameter'] and \
                (int(job.executable.version.split('.')[0]) < 2016 or
                 (int(job.executable.version.split('.')[0]) == 2016 and
                  int(job.executable.version.split('.')[1]) < 11)):
            job.input.control['dump_modify'] = job.input.control['dump_modify'].replace(' line ', ' ')
        if not all(job.structure.pbc):
            job.input.control['boundary'] = ' '.join(['p' if coord else 'f' for coord in job.structure.pbc])
        job._set_selective_dynamics()
        job.input.control.write_file(file_name="control.inp", cwd=job.working_directory)
        job.input.potential.write_file(file_name="potential.inp", cwd=job.working_directory)
        job.input.potential.copy_pot_files(job.working_directory)

    def collect_output(self, job):
        """

        Returns:

        """
        job.input.from_hdf(job._hdf5)
        if os.path.isfile(job.job_file_name(file_name="dump.h5", cwd=job.working_directory)):
            self._collect_h5md_file(job=job, file_name="dump.h5", cwd=job.working_directory)
        else:
            self._collect_dump_file(job=job, file_name="dump.out", cwd=job.working_directory)
        self._collect_output_log(job=job, file_name="log.lammps", cwd=job.working_directory)
        final_structure = job.get_final_structure()
        with job.project_hdf5.open("output") as hdf_output:
            final_structure.to_hdf(hdf_output)

    @staticmethod
    def _get_lammps_structure(job, structure=None, cutoff_radius=None):
        lmp_structure = LammpsStructure()
        lmp_structure.potential = job.input.potential
        lmp_structure.atom_type = job.input.control["atom_style"]
        if cutoff_radius is not None:
            lmp_structure.cutoff_radius = cutoff_radius
        else:
            lmp_structure.cutoff_radius = job.cutoff_radius
        lmp_structure.el_eam_lst = job.input.potential.get_element_lst()
        if structure is not None:
            lmp_structure.structure = structure
        else:
            lmp_structure.structure = job.structure
        if not set(lmp_structure.structure.get_species_symbols()).issubset(set(lmp_structure.el_eam_lst)):
            raise ValueError('The selected potentials do not support the given combination of elements.')
        return lmp_structure

    @staticmethod
    def _collect_dump_file(job, file_name="dump.out", cwd=None):
        """
        general purpose routine to extract static from a lammps dump file

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = job.job_file_name(file_name=file_name, cwd=cwd)
        tag_dict = {"ITEM: TIMESTEP": {"arg": "0",
                                       "rows": 1,
                                       "h5": "time"},
                    # "ITEM: NUMBER OF ATOMS": {"arg": "0",
                    #                          "rows": 1,
                    #                          "h5": "number_of_atoms"},
                    "ITEM: BOX BOUNDS": {"arg": "0",
                                         "rows": 3,
                                         "h5": "cells",
                                         "func": to_amat},
                    "ITEM: ATOMS": {"arg": ":,:",
                                    "rows": len(job.structure),
                                    "splitArg": True}
                    }

        h5_dict = {"id": "id",
                   "type": "type",
                   "xsu": "coord_xs",
                   "ysu": "coord_ys",
                   "zsu": "coord_zs",
                   "f_ave[1]": "coord_xs",
                   "f_ave[2]": "coord_ys",
                   "f_ave[3]": "coord_zs",
                   "fx": "force_x",
                   "fy": "force_y",
                   "fz": "force_z",
                   }

        lammps_dict = None

        lf = Logstatus()
        lf.extract_file(file_name=file_name,
                        tag_dict=tag_dict,
                        h5_dict=h5_dict,
                        key_dict=lammps_dict)
        lf.combine_xyz('force_x', 'force_y', 'force_z', 'forces')
        lf.combine_xyz('coord_xs', 'coord_ys', 'coord_zs', 'positions')

        prism = UnfoldingPrism(job.structure.cell, digits=15)

        rel_positions = list()

        for ind, (pos, forc, cel) in enumerate(
                zip(lf.status_dict["positions"], lf.status_dict["forces"], lf.status_dict["cells"])):
            cell = cel[1]
            positions = pos[1]
            forces = forc[1]

            # rotation matrix from lammps(unfolded) cell to original cell
            rotation_lammps2orig = np.linalg.inv(prism.R)

            # convert from scaled positions to absolute in lammps cell
            positions = np.array([np.dot(cell.T, r) for r in positions])
            # rotate positions from lammps to original
            positions_atoms = np.array([np.dot(np.array(r), rotation_lammps2orig) for r in positions])

            # rotate forces from lammps to original cell
            forces_atoms = np.array([np.dot(np.array(f), rotation_lammps2orig) for f in forces])

            # unfold cell
            cell = prism.unfold_cell(cell)
            # rotate cell from unfolded lammps to original
            cell_atoms = np.array([np.dot(np.array(f), rotation_lammps2orig) for f in cell])

            lf.status_dict["positions"][ind][1] = positions_atoms

            rel_positions.append(np.dot(positions_atoms, np.linalg.inv(cell_atoms)))

            lf.status_dict["forces"][ind][1] = forces_atoms
            lf.status_dict["cells"][ind][1] = cell_atoms

        del lf.status_dict['id']
        del lf.status_dict['type']
        unwrapped_rel_pos = unwrap_coordinates(positions=np.array(rel_positions), is_relative=True)
        unwrapped_pos = list()
        # print(np.shape(unwrapped_rel_pos))
        for i, cell in enumerate(lf.status_dict["cells"]):
            unwrapped_pos.append(np.dot(np.array(unwrapped_rel_pos[i]), cell[1]))
        lf.status_dict["unwrapped_positions"] = list()
        for pos in unwrapped_pos:
            lf.status_dict["unwrapped_positions"].append([[0], pos])
        with job.project_hdf5.open("output/generic") as hdf_output:
            lf.to_hdf(hdf_output)
        return lf

    # TODO: make rotation of all vectors back to the original as in self.collect_dump_file
    @staticmethod
    def _collect_h5md_file(job, file_name="dump.h5", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = job.job_file_name(file_name=file_name, cwd=cwd)
        with h5py.File(file_name) as h5md:
            positions = [pos_i.tolist() for pos_i in h5md['/particles/all/position/value']]
            time = [time_i.tolist() for time_i in h5md['/particles/all/position/step']]
            forces = [for_i.tolist() for for_i in h5md['/particles/all/force/value']]
            # following the explanation at: http://nongnu.org/h5md/h5md.html
            cell = [np.eye(3) * np.array(cell_i.tolist()) for cell_i in h5md['/particles/all/box/edges/value']]
        with job.project_hdf5.open("output/generic") as h5_file:
            h5_file['forces'] = np.array(forces)
            h5_file['positions'] = np.array(positions)
            h5_file['time'] = np.array(time)
            h5_file['cells'] = cell

    @staticmethod
    def _collect_errors(job, file_name, cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = job.job_file_name(file_name=file_name, cwd=cwd)
        error = extract_data_from_file(file_name, tag="ERROR", num_args=1000)
        if len(error) > 0:
            error = " ".join(error[0])
            raise RuntimeError("Run time error occurred: " + str(error))
        else:
            return True

    def _collect_output_log(self, job, file_name="log.lammps", cwd=None):
        """
        general purpose routine to extract static from a lammps log file

        Args:
            file_name:
            cwd:

        Returns:

        """
        log_file = job.job_file_name(file_name=file_name, cwd=cwd)
        self._collect_errors(job, file_name)
        # log_fileName = "logfile.out"
        # self.collect_errors(log_fileName)
        attr = job.input.control.dataset["Parameter"]
        tag_dict = {"Loop time of": {"arg": "0",
                                     "type": "float",
                                     "h5": "time_loop"},

                    "Memory usage per processor": {"arg": "1",
                                                   "h5": "memory"}
                    }

        if 'minimize' in attr:
            tag_dict["Step Temp PotEng TotEng Pxx Pxy Pxz Pyy Pyz Pzz Volume"] = {"arg": ":,:",
                                                                                  "rows": "Loop",
                                                                                  "splitTag": True}

        elif 'run' in attr:
            num_iterations = ast.literal_eval(extract_data_from_file(log_file, tag="run")[0])
            num_thermo = ast.literal_eval(extract_data_from_file(log_file, tag="thermo")[0])
            num_iterations = num_iterations / num_thermo + 1

            tag_dict["thermo_style custom"] = {"arg": ":",
                                               "type": "str",
                                               "h5": "thermo_style"}
            tag_dict["$thermo_style custom"] = {"arg": ":,:",
                                                "rows": num_iterations,
                                                "splitTag": True}
            # print "tag_dict: ", tag_dict

        h5_dict = {"Step": "steps",
                   "Temp": "temperatures",
                   "PotEng": "energy_pot",
                   "TotEng": "energy_tot",
                   "Pxx": "pressure_x",
                   "Pxy": "pressure_xy",
                   "Pxz": "pressure_xz",
                   "Pyy": "pressure_y",
                   "Pyz": "pressure_yz",
                   "Pzz": "pressure_z",
                   "Volume": "volume",
                   "E_pair": "E_pair",
                   "E_mol": "E_mol"
                   }

        lammps_dict = {"step": "Step",
                       "temp": "Temp",
                       "pe": "PotEng",
                       "etotal": "TotEng",
                       "pxx": "Pxx",
                       "pxy": "Pxy",
                       "pxz": "Pxz",
                       "pyy": "Pyy",
                       "pyz": "Pyz",
                       "pzz": "Pzz",
                       "vol": "Volume"
                       }

        lf = Logstatus()
        lf.extract_file(file_name=log_file,
                        tag_dict=tag_dict,
                        h5_dict=h5_dict,
                        key_dict=lammps_dict)

        lf.store_as_vector = ['energy_tot', 'temperatures', 'steps', 'volume', 'energy_pot']
        # print ("lf_keys: ", lf.status_dict['energy_tot'])

        lf.combine_mat('pressure_x', 'pressure_xy', 'pressure_xz',
                       'pressure_y', 'pressure_yz', 'pressure_z', 'pressures')
        lf.convert_unit('pressures', 0.0001)  # bar -> GPa

        if 'minimize' not in attr:
            if 'thermo_style' in lf.status_dict.keys():
                del lf.status_dict['thermo_style']
        if 'time_loop' in lf.status_dict.keys():
            del lf.status_dict['time_loop']
        if 'memory' in lf.status_dict.keys():
            del lf.status_dict['memory']
        with job.project_hdf5.open("output/generic") as hdf_output:
            lf.to_hdf(hdf_output)


class InteractiveLammpsInterface(LammpsInterface, AtomisticInteractiveInterface):
    def __init__(self):
        super(InteractiveLammpsInterface, self).__init__()
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

    def _interactive_lib_command(self, job, command):
        job._logger.debug('Lammps library: ' + command)
        self._interactive_library.command(command)

    def interactive_positions_getter(self, job):
        positions = np.reshape(np.array(self._interactive_library.gather_atoms("x", 1, 3)),
                               (len(job.structure), 3))
        if np.matrix.trace(self._interactive_prism.R) != 3:
            positions = np.dot(positions, self._interactive_prism.R.T)
        return positions.tolist()

    def interactive_positions_setter(self, job, positions):
        if np.matrix.trace(self._interactive_prism.R) != 3:
            positions = np.array(positions).reshape(-1, 3)
            positions = np.dot(positions, self._interactive_prism.R)
        positions = np.array(positions).flatten()
        if job.server.run_mode.interactive_non_modal:
            self._interactive_library.scatter_atoms("x", 1, 3, positions)
        else:
            self._interactive_library.scatter_atoms("x", 1, 3, (len(positions) * c_double)(*positions))
        self._interactive_lib_command(job=job, command='change_box all remap')

    def interactive_cells_getter(self, job):
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

    def interactive_cells_setter(self, job, cell):
        self._interactive_prism = UnfoldingPrism(cell)
        lx, ly, lz, xy, xz, yz = self._interactive_prism.get_lammps_prism()
        if np.matrix.trace(self._interactive_prism.R) != 3:
            print('Warning: setting upper trangular matrix might slow down the calculation')
        if abs(xy) + abs(xz) + abs(yz) > 1.0e-6:
            if job.structure._is_scaled:
                self._interactive_lib_command(job=job,
                                              command='change_box all x final 0 %f y final 0 %f z final 0 %f \
                     xy final %f xz final %f yz final %f triclinic remap units box' % (lx, ly, lz, xy, xz, yz))
            else:
                self._interactive_lib_command(job=job,
                                              command='change_box all x final 0 %f y final 0 %f z final 0 %f \
                    xy final %f xz final %f yz final %f triclinic units box' % (lx, ly, lz, xy, xz, yz))
        else:
            if job.structure._is_scaled:
                self._interactive_lib_command(job=job,
                                              command='change_box all x final 0 %f y final 0 %f z final 0 %f \
                                              remap units box' % (lx, ly, lz))
            else:
                self._interactive_lib_command(job=job,
                                              command='change_box all x final 0 %f y final 0 %f z final 0 %f \
                                              units box' % (lx, ly, lz))

    def interactive_indices_setter(self, job, indices):
        el_struct_lst = self.structure_current.get_species_symbols()
        el_obj_lst = self.structure_current.get_species_objects()
        el_eam_lst = job.input.potential.get_element_lst()
        el_dict = {}
        for id_eam, el_eam in enumerate(el_eam_lst):
            if el_eam in el_struct_lst:
                id_el = list(el_struct_lst).index(el_eam)
                el = el_obj_lst[id_el]
                el_dict[el] = id_eam + 1
        elem_all = np.array([el_dict[self.structure_current.species[el]] for el in indices])
        if job.server.run_mode.interactive_non_modal:
            self._interactive_library.scatter_atoms('type', 0, 1, elem_all)
        else:
            self._interactive_library.scatter_atoms('type', 0, 1, (len(elem_all) * c_int)(*elem_all))

    def interactive_volume_getter(self, job):
        return self._interactive_library.get_thermo('vol')

    def interactive_forces_getter(self, job):
        ff = np.reshape(np.array(self._interactive_library.gather_atoms("f", 1, 3)), (len(job.structure), 3))
        if np.matrix.trace(self._interactive_prism.R) != 3:
            ff = np.dot(ff, self._interactive_prism.R.T)
        return ff.tolist()

    def _interactive_lammps_input(self, job):
        del job.input.control['dump']
        del job.input.control['dump_modify']
        for key, value in zip(job.input.control.dataset['Parameter'], job.input.control.dataset['Value']):
            if key in ['read_data', 'units', 'dimension', 'boundary', 'atom_style', 'atom_modify', 'include', 'run',
                       'minimize']:
                continue
            else:
                self._interactive_lib_command(job=job, command=key + ' ' + str(value))

    def _interactive_set_potential(self, job):
        potential_lst = []
        if job.input.potential.files is not None:
            for potential in job.input.potential.files:
                potential_lst.append([potential.split('/')[-1], potential])
        for line in job.input.potential.get_string_lst():
            if len(line) > 2:
                for potential in potential_lst:
                    if potential[0] in line:
                        line = line.replace(potential[0], potential[1])
                self._interactive_lib_command(job=job, command=line.split('\n')[0])

    def _reset_interactive_run_command(self, job):
        df = pd.DataFrame(job.input.control.dataset)
        self._interactive_run_command = " ".join(df.T[df.index[-1]].values)

    def interactive_open(self, job):
        if job.server.run_mode.interactive_non_modal:
            self._interactive_library = LammpsLibrary()
        else:
            self._interactive_library = lammps()
        if not all(job.structure.pbc):
            job.input.control['boundary'] = ' '.join(['p' if coord else 'f' for coord in job.structure.pbc])
        self._reset_interactive_run_command(job=job)
        self.interactive_structure_setter(job=job, structure=job.structure)

    def interactive_structure_setter(self, job, structure):
        self._interactive_lib_command(job=job, command='clear')
        job._set_selective_dynamics()
        self._interactive_lib_command(job=job, command='units ' + job.input.control['units'])
        self._interactive_lib_command(job=job, command='dimension ' + str(job.input.control['dimension']))
        self._interactive_lib_command(job=job, command='boundary ' + job.input.control['boundary'])
        self._interactive_lib_command(job=job, command='atom_style ' + job.input.control['atom_style'])
        self._interactive_lib_command(job=job, command="atom_modify map array")
        self._interactive_prism = UnfoldingPrism(structure.cell)
        if np.matrix.trace(self._interactive_prism.R) != 3:
            print('Warning: setting upper trangular matrix might slow down the calculation')
        xhi, yhi, zhi, xy, xz, yz = self._interactive_prism.get_lammps_prism()
        if self._interactive_prism.is_skewed():
            self._interactive_lib_command(job=job, command='region 1 prism' +
                                          ' 0.0 ' + str(xhi) + ' 0.0 ' + str(yhi) + ' 0.0 ' + str(zhi) +
                                          ' ' + str(xy) + ' ' + str(xz) + ' ' + str(yz) + ' units box')
        else:
            self._interactive_lib_command(job=job, command='region 1 block' +
                                          ' 0.0 ' + str(xhi) + ' 0.0 ' + str(yhi) + ' 0.0 ' + str(zhi) + ' units box')
        el_struct_lst = job.structure.get_species_symbols()
        el_obj_lst = job.structure.get_species_objects()
        el_eam_lst = job.input.potential.get_element_lst()
        self._interactive_lib_command(job=job, command='create_box ' + str(len(el_eam_lst)) + ' 1')
        el_dict = {}
        for id_eam, el_eam in enumerate(el_eam_lst):
            if el_eam in el_struct_lst:
                id_el = list(el_struct_lst).index(el_eam)
                el = el_obj_lst[id_el]
                el_dict[el] = id_eam + 1
                self._interactive_lib_command(job=job, command='mass {0:3d} {1:f}'.format(id_eam + 1, el.AtomicMass))
            else:
                self._interactive_lib_command(job=job, command='mass {0:3d} {1:f}'.format(id_eam + 1, 1.00))
        self._interactive_lib_command(job=job, command='create_atoms 1 random ' + str(len(structure)) + ' 12345 1')
        positions = structure.positions.flatten()
        elem_all = np.array([el_dict[el] for el in structure.get_chemical_elements()])
        if job.server.run_mode.interactive_non_modal:
            self._interactive_library.scatter_atoms("x", 1, 3, positions)
            self._interactive_library.scatter_atoms('type', 0, 1, elem_all)
        else:
            self._interactive_library.scatter_atoms("x", 1, 3, (len(positions) * c_double)(*positions))
            self._interactive_library.scatter_atoms('type', 0, 1, (len(elem_all) * c_int)(*elem_all))
        self._interactive_lib_command(job=job, command='change_box all remap')
        self._interactive_lammps_input(job=job)
        self._interactive_set_potential(job=job)

    def update_potential(self, job):
        self._interactive_lib_command(job=job, command=job.potential.Config[0][0])
        self._interactive_lib_command(job=job, command=job.potential.Config[0][1])

    def interactive_energy_pot_getter(self, job):
        return self._interactive_library.get_thermo("pe")

    def interactive_energy_tot_getter(self, job):
        return self._interactive_library.get_thermo("etotal")

    def interactive_steps_getter(self, job):
        return self._interactive_library.get_thermo("step")

    def interactive_temperatures_getter(self, job):
        return self._interactive_library.get_thermo("temp")

    def interactive_stress_getter(self, job):
        '''
        This gives back an Nx3x3 np array of stress/atom defined in
        http://lammps.sandia.gov/doc/compute_stress_atom.html
        Keep in mind that it is stress*volume in eV.
        Further discussion can be found on the website above.
        '''
        if not 'stress' in self.interactive_cache.keys():
            self._interactive_lib_command(job=job, command='compute st all stress/atom NULL')
            self._interactive_lib_command(job=job, command='run 0')
            self.interactive_cache['stress'] = []
        ss = np.array([self._interactive_library.extract_compute('st', 1, 2)[i][j + (j != k) * (k + 2)]
                       for i in range(len(job.structure))
                       for j in range(3)
                       for k in range(3)]).reshape(len(job.structure), 3, 3) / 1.602e6
        if np.matrix.trace(self._interactive_prism.R) != 3:
            ss = np.dot(np.dot(self._interactive_prism.R, ss), self._interactive_prism.R.T)
        return ss

    def interactive_pressures_getter(self, job):
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

    def run_if_interactive(self, job):
        if job._generic_input['calc_mode'] == 'md':
            job.input.control['run'] = job._generic_input['n_print']
            super(InteractiveLammpsInterface, self).run_if_interactive(job=job)
            self._reset_interactive_run_command(job=job)

            counter = 0
            iteration_max = int(job._generic_input['n_ionic_steps'] / job._generic_input['n_print'])
            while counter < iteration_max:
                self._interactive_lib_command(job=job, command=self._interactive_run_command)
                self.interactive_collect(job=job)
                counter += 1

        else:
            super(InteractiveLammpsInterface, self).run_if_interactive(job=job)
            self._interactive_lib_command(job=job, command=self._interactive_run_command)
            self.interactive_collect(job=job)

    def run_if_interactive_non_modal(self, job):
        if not self._interactive_fetch_completed:
            print('Warning: interactive_fetch being effectuated')
            self.interactive_fetch()
        super(InteractiveLammpsInterface, self).run_if_interactive(job=job)
        self._interactive_lib_command(job=job, command=self._interactive_run_command)
        self._interactive_fetch_completed = False

    def interactive_fetch(self, job):
        if self._interactive_fetch_completed and job.server.run_mode.interactive_non_modal:
            print('First run and then fetch')
        else:
            self.interactive_collect(job=job)
            job._logger.debug('interactive run - done')

    def interactive_close(self, job):
        if self.interactive_is_activated():
            self._interactive_library.close()
            with job.project_hdf5.open("output") as h5:
                if 'interactive' in h5.list_groups():
                    for key in h5['interactive'].list_nodes():
                        h5['generic/' + key] = h5['interactive/' + key]
            super(InteractiveLammpsInterface, self).interactive_close(job=job)


def to_amat(l_list):
    """

    Args:
        l_list:

    Returns:

    """
    lst = np.reshape(l_list, -1)
    if len(lst) == 9:
        xlo_bound, xhi_bound, xy, ylo_bound, yhi_bound, xz, zlo_bound, zhi_bound, yz = lst

    elif len(lst) == 6:
        xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound = lst
        xy, xz, yz = 0., 0., 0.
    else:
        raise ValueError("This format for amat not yet implemented: " + str(len(lst)))

    # > xhi_bound - xlo_bound = xhi -xlo  + MAX(0.0, xy, xz, xy + xz) - MIN(0.0, xy, xz, xy + xz)
    # > xhili = xhi -xlo   = xhi_bound - xlo_bound - MAX(0.0, xy, xz, xy + xz) + MIN(0.0, xy, xz, xy + xz)
    xhilo = (xhi_bound - xlo_bound) - max([0.0, xy, xz, xy + xz]) + min([0.0, xy, xz, xy + xz])

    # > yhilo = yhi -ylo = yhi_bound -ylo_bound - MAX(0.0, yz) + MIN(0.0, yz)
    yhilo = (yhi_bound - ylo_bound) - max([0.0, yz]) + min([0.0, yz])

    # > zhi - zlo = zhi_bound- zlo_bound
    zhilo = (zhi_bound - zlo_bound)

    cell = [[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]]
    return cell
