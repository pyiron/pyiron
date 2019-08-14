# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, unicode_literals
import os
import posixpath

import sys
import h5py
import numpy as np
import pandas as pd
import warnings
from io import StringIO

from pyiron.lammps.potential import LammpsPotentialFile, PotentialAvailable
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.base.settings.generic import Settings
from pyiron.base.pyio.parser import Logstatus, extract_data_from_file
from pyiron.lammps.control import LammpsControl
from pyiron.lammps.potential import LammpsPotential
from pyiron.lammps.structure import LammpsStructure, UnfoldingPrism
from pyiron.atomistics.md_analysis.trajectory_analysis import unwrap_coordinates

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


s = Settings()


class LammpsBase(AtomisticGenericJob):
    """
    Class to setup and run and analyze LAMMPS simulations which is a derivative of
    atomistics.job.generic.GenericJob. The functions in these modules are written in such the function names and
    attributes are very generic (get_structure(), molecular_dynamics(), version) but the functions are written to handle
    LAMMPS specific input/output.

    Args:
        project (pyiron.project.Project instance):  Specifies the project path among other attributes
        job_name (str): Name of the job

    Attributes:
        input (lammps.Input instance): Instance which handles the input
    """

    def __init__(self, project, job_name):
        super(LammpsBase, self).__init__(project, job_name)
        self.input = Input()
        self._cutoff_radius = None
        self._is_continuation = None
        self._compress_by_default = True
        s.publication_add(self.publication)

    @property
    def cutoff_radius(self):
        """

        Returns:

        """
        return self._cutoff_radius

    @cutoff_radius.setter
    def cutoff_radius(self, cutoff):
        """

        Args:
            cutoff:

        Returns:

        """
        self._cutoff_radius = cutoff

    @property
    def potential(self):
        """
        Execute view_potential() or list_potential() in order to see the pre-defined potential files

        Returns:

        """
        return self.input.potential.df

    @potential.setter
    def potential(self, potential_filename):
        """
        Execute view_potential() or list_potential() in order to see the pre-defined potential files

        Args:
            potential_filename:

        Returns:

        """
        if sys.version_info.major == 2:
            stringtypes = (str, unicode)
        else:
            stringtypes = str
        if isinstance(potential_filename, stringtypes):
            if '.lmp' in potential_filename:
                potential_filename = potential_filename.split('.lmp')[0]
            potential_db = LammpsPotentialFile()
            potential = potential_db.find_by_name(potential_filename)
        elif isinstance(potential_filename, pd.DataFrame):
            potential = potential_filename
        else:
            raise TypeError('Potentials have to be strings or pandas dataframes.')
        self.input.potential.df = potential
        for val in ["units", "atom_style", "dimension"]:
            v = self.input.potential[val]
            if v is not None:
                self.input.control[val] = v
        self.input.potential.remove_structure_block()

    @property
    def potential_available(self):
        return PotentialAvailable(list_of_potentials=self.potential_list)

    @property
    def potential_list(self):
        """
        List of interatomic potentials suitable for the current atomic structure.

        use self.potentials_view() to get more details.

        Returns:
            list: potential names
        """
        return self.list_potentials()

    @property
    def potential_view(self):
        """
        List all interatomic potentials for the current atomistic sturcture including all potential parameters.

        To quickly get only the names of the potentials you can use: self.potentials_list()

        Returns:
            pandas.Dataframe: Dataframe including all potential parameters.
        """
        return self.view_potentials()

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(LammpsBase, self).set_input_to_read_only()
        self.input.control.read_only = True
        self.input.potential.read_only = True

    def validate_ready_to_run(self):
        """

        Returns:

        """
        super(LammpsBase, self).validate_ready_to_run()
        if self.potential is None:
            raise ValueError('This job does not contain a valid potential: {}'.format(self.job_name))

    def get_potentials_for_structure(self):
        """

        Returns:

        """
        return self.list_potentials()

    def get_final_structure(self):
        """

        Returns:

        """
        warnings.warn("get_final_structure() is deprecated - please use get_structure() instead.", DeprecationWarning)
        return self.get_structure(iteration_step=-1)

    def view_potentials(self):
        """
        List all interatomic potentials for the current atomistic sturcture including all potential parameters.

        To quickly get only the names of the potentials you can use: self.potentials_list()

        Returns:
            pandas.Dataframe: Dataframe including all potential parameters.
        """
        from pyiron.lammps.potential import LammpsPotentialFile
        if not self.structure:
            raise ValueError('No structure set.')
        list_of_elements = set(self.structure.get_chemical_symbols())
        list_of_potentials = LammpsPotentialFile().find(list_of_elements)
        if list_of_potentials is not None:
            return list_of_potentials
        else:
            raise TypeError('No potentials found for this kind of structure: ', str(list_of_elements))

    def list_potentials(self):
        """
        List of interatomic potentials suitable for the current atomic structure.

        use self.potentials_view() to get more details.

        Returns:
            list: potential names
        """
        return list(self.view_potentials()['Name'].values)

    def enable_h5md(self):
        """

        Returns:

        """
        del self.input.control['dump_modify']
        del self.input.control['dump']
        self.input.control['dump'] = '1 all h5md ${dumptime} dump.h5 position force create_group yes'

    def write_input(self):
        """
        Call routines that generate the code specific input files

        Returns:

        """
        if self.structure is None:
            raise ValueError("Input structure not set. Use method set_structure()")
        lmp_structure = self._get_lammps_structure(structure=self.structure, cutoff_radius=self.cutoff_radius)
        lmp_structure.write_file(file_name="structure.inp", cwd=self.working_directory)
        version_int_lst = self._get_executable_version_number()
        if version_int_lst is not None and 'dump_modify' in self.input.control._dataset['Parameter'] and \
                (version_int_lst[0] < 2016 or (version_int_lst[0] == 2016 and version_int_lst[1] < 11)):
            self.input.control['dump_modify'] = self.input.control['dump_modify'].replace(' line ', ' ')
        if not all(self.structure.pbc):
            self.input.control['boundary'] = ' '.join(['p' if coord else 'f' for coord in self.structure.pbc])
        self._set_selective_dynamics()
        self.input.control.write_file(file_name="control.inp", cwd=self.working_directory)
        self.input.potential.write_file(file_name="potential.inp", cwd=self.working_directory)
        self.input.potential.copy_pot_files(self.working_directory)

    def _get_executable_version_number(self):
        """
        Get the version of the executable

        Returns:
            list: List of integers defining the version number
        """
        if self.executable.version:
            return [l for l in [[int(i) for i in s.split('.') if i.isdigit()]
                                for s in self.executable.version.split('_')]
                    if len(l) > 0][0]
        else:
            return None

    @property
    def publication(self):
        return {'lammps': {'lammps': {'title': 'Fast Parallel Algorithms for Short-Range Molecular Dynamics',
                                      'journal': 'Journal of Computational Physics',
                                      'volume': '117',
                                      'number': '1',
                                      'pages': '1-19',
                                      'year': '1995',
                                      'issn': '0021-9991',
                                      'doi': '10.1006/jcph.1995.1039',
                                      'url': 'http://www.sciencedirect.com/science/article/pii/S002199918571039X',
                                      'author': ['Steve Plimpton']}}}

    def collect_output(self):
        """

        Returns:

        """
        self.input.from_hdf(self._hdf5)
        if os.path.isfile(self.job_file_name(file_name="dump.h5", cwd=self.working_directory)):
            self.collect_h5md_file(file_name="dump.h5", cwd=self.working_directory)
        else:
            self.collect_dump_file(file_name="dump.out", cwd=self.working_directory)
        self.collect_output_log(file_name="log.lammps", cwd=self.working_directory)
        final_structure = self.get_structure(iteration_step=-1)
        with self.project_hdf5.open("output") as hdf_output:
            final_structure.to_hdf(hdf_output)

    def convergence_check(self):
        if self._generic_input['calc_mode'] == 'minimize':
            if self._generic_input['max_iter'] + 1 <= len(self['output/generic/energy_tot']) or \
                    len([l for l in self['log.lammps'] if 'linesearch alpha is zero' in l]) != 0:
                return False
            else:
                return True
        else:
            return True

    def collect_logfiles(self):
        """

        Returns:

        """
        return

    # TODO: make rotation of all vectors back to the original as in self.collect_dump_file
    def collect_h5md_file(self, file_name="dump.h5", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
        with h5py.File(file_name, 'r', libver='latest', swmr=True) as h5md:
            positions = [pos_i.tolist() for pos_i in h5md['/particles/all/position/value']]
            time = [time_i.tolist() for time_i in h5md['/particles/all/position/step']]
            forces = [for_i.tolist() for for_i in h5md['/particles/all/force/value']]
            # following the explanation at: http://nongnu.org/h5md/h5md.html
            cell = [np.eye(3) * np.array(cell_i.tolist()) for cell_i in h5md['/particles/all/box/edges/value']]
        with self.project_hdf5.open("output/generic") as h5_file:
            h5_file['forces'] = np.array(forces)
            h5_file['positions'] = np.array(positions)
            h5_file['time'] = np.array(time)
            h5_file['cells'] = cell

    def collect_errors(self, file_name, cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
        error = extract_data_from_file(file_name, tag="ERROR", num_args=1000)
        if len(error) > 0:
            error = " ".join(error[0])
            raise RuntimeError("Run time error occurred: " + str(error))
        else:
            return True

    def collect_output_log(self, file_name="log.lammps", cwd=None):
        """
        general purpose routine to extract static from a lammps log file

        Args:
            file_name:
            cwd:

        Returns:

        """
        self.collect_errors(file_name=file_name, cwd=cwd)
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
        with open(file_name, 'r') as f:
            f = f.readlines()
            l_start = np.where([line.startswith('Step') for line in f])[0]
            l_end = np.where([line.startswith('Loop') for line in f])[0]
            if len(l_start)>len(l_end):
                l_end = np.append(l_end, [None])
            if sys.version_info >= (3,):
                df = [pd.read_csv(StringIO('\n'.join(f[llst:llen])),
                                  delim_whitespace=True) for llst, llen in zip(l_start, l_end)]
            else:
                df = [pd.read_csv(StringIO(unicode('\n'.join(f[llst:llen]))),
                                  delim_whitespace=True) for llst, llen in zip(l_start, l_end)]
        df = df[-1]

        h5_dict = {"Step": "steps",
                   "Temp": "temperature",
                   "PotEng": "energy_pot",
                   "TotEng": "energy_tot",
                   "Volume": "volume"}

        df = df.rename(index=str, columns=h5_dict)
        pressures = np.stack((df.Pxx, df.Pxy, df.Pxz,
                              df.Pxy, df.Pyy, df.Pyz,
                              df.Pxz, df.Pyz, df.Pzz), axis=-1).reshape(-1, 3, 3)
        pressures *= 0.0001  # bar -> GPa
        df = df.drop(columns=df.columns[((df.columns.str.len() == 3) & df.columns.str.startswith('P'))])
        df['pressures'] = pressures.tolist()

        with self.project_hdf5.open("output/generic") as hdf_output:
            # This is a hack for backward comparability
            for k,v in df.items():
                hdf_output[k] = np.array(v)

    def calc_minimize(self, e_tol=0.0, f_tol=1e-2, max_iter=100000, pressure=None, n_print=100):
        """

        Args:
            e_tol:
            f_tol:
            max_iter:
            pressure:
            n_print:

        Returns:

        """
        super(LammpsBase, self).calc_minimize(e_tol=e_tol, f_tol=f_tol, max_iter=max_iter, pressure=pressure,
                                              n_print=n_print)
        self.input.control.calc_minimize(e_tol=e_tol, f_tol=f_tol, max_iter=max_iter, pressure=pressure,
                                         n_print=n_print)

    def calc_static(self):
        """

        Returns:

        """
        super(LammpsBase, self).calc_static()
        self.input.control.calc_static()

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, time_step=1.0, n_print=100,
                temperature_damping_timescale=100.0, pressure_damping_timescale=1000.0, seed=None, tloop=None, initial_temperature=None,
                langevin=False, delta_temp=None, delta_press=None):
        """
        Set an MD calculation within LAMMPS. Nosé Hoover is used by default.

        Args:
            temperature (None/float): Target temperature. If set to None, an NVE calculation is performed.
                                      It is required when the pressure is set or langevin is set
            pressure (None/float): Target pressure. If set to None, an NVE or an NVT calculation is performed.
                                   (This tag will allow for a list in the future as it is done for calc_minimize())
            n_ionic_steps (int): Number of ionic steps
            time_step (float): Step size between two steps. In fs if units==metal
            n_print (int):  Print frequency
            temperature_damping_timescale (float): The time associated with the thermostat adjusting the temperature.
                                                   (In fs. After rescaling to appropriate time units, is equivalent to
                                                   Lammps' `Tdamp`.)
            pressure_damping_timescale (float): The time associated with the barostat adjusting the temperature.
                                                (In fs. After rescaling to appropriate time units, is equivalent to
                                                Lammps' `Pdamp`.)
            seed (int):  Seed for the random number generation (required for the velocity creation)
            tloop:
            initial_temperature (None/float):  Initial temperature according to which the initial velocity field
                                               is created. If None, the initial temperature will be twice the target
                                               temperature (which would go immediately down to the target temperature
                                               as described in equipartition theorem). If 0, the velocity field is not
                                               initialized (in which case  the initial velocity given in structure will
                                               be used). If any other number is given, this value is going to be used
                                               for the initial temperature.
            langevin (bool): (True or False) Activate Langevin dynamics
            delta_temp (float): Thermostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
            delta_press (float): Barostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
        """
        if self.server.run_mode.interactive_non_modal:
            warnings.warn('calc_md() is not implemented for the non modal interactive mode use calc_static()!')
        super(LammpsBase, self).calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps,
                                        time_step=time_step, n_print=n_print,
                                        temperature_damping_timescale=temperature_damping_timescale,
                                        pressure_damping_timescale=pressure_damping_timescale,
                                        seed=seed, tloop=tloop, initial_temperature=initial_temperature,
                                        langevin=langevin)
        self.input.control.calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps,
                                   time_step=time_step, n_print=n_print, temperature_damping_timescale=temperature_damping_timescale,
                                   pressure_damping_timescale=pressure_damping_timescale,
                                   seed=seed, tloop=tloop, initial_temperature=initial_temperature, langevin=langevin,
                                   delta_temp=delta_temp, delta_press=delta_press, job_name=self.job_name)

    # define hdf5 input and output
    def to_hdf(self, hdf=None, group_name=None):
        """

        Args:
            hdf:
            group_name:

        Returns:

        """
        super(LammpsBase, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()
        self.input.to_hdf(self._hdf5)

    def from_hdf(self, hdf=None, group_name=None):  # TODO: group_name should be removed
        """

        Args:
            hdf:
            group_name:

        Returns:

        """
        super(LammpsBase, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()
        self.input.from_hdf(self._hdf5)

    def write_restart_file(self, filename="restart.out"):
        """

        Args:
            filename:

        Returns:

        """
        self.input.control.modify(write_restart=filename, append_if_not_present=True)

    def compress(self, files_to_compress=None):
        """
        Compress the output files of a job object.

        Args:
            files_to_compress (list):
        """
        if files_to_compress is None:
            files_to_compress = [f for f in list(self.list_files()) if f not in ["restart.out"]]
        super(LammpsBase, self).compress(files_to_compress=files_to_compress)

    def read_restart_file(self, filename="restart.out"):
        """

        Args:
            filename:

        Returns:

        """
        self._is_continuation = True
        self.input.control.set(read_restart=filename)
        self.input.control['reset_timestep'] = 0
        self.input.control.remove_keys(['dimension', 'read_data', 'boundary', 'atom_style', 'velocity'])

    def collect_dump_file(self, file_name="dump.out", cwd=None):
        """
        general purpose routine to extract static from a lammps dump file

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
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
                                    "rows": len(self.structure),
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

        prism = UnfoldingPrism(self.structure.cell, digits=15)

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
        with self.project_hdf5.open("output/generic") as hdf_output:
            lf.to_hdf(hdf_output)
        return lf

    # Outdated functions:
    def set_potential(self, file_name):
        """

        Args:
            file_name:

        Returns:

        """
        print('This function is outdated use the potential setter instead!')
        self.potential = file_name

    def next(self, snapshot=-1, job_name=None, job_type=None):
        """
        Restart a new job created from an existing Lammps calculation.
        Args:
            project (pyiron.project.Project instance): Project instance at which the new job should be created
            snapshot (int): Snapshot of the calculations which would be the initial structure of the new job
            job_name (str): Job name
            job_type (str): Job type. If not specified a Lammps job type is assumed

        Returns:
            new_ham (lammps.lammps.Lammps instance): New job
        """
        return super(LammpsBase, self).restart(snapshot=snapshot, job_name=job_name, job_type=job_type)

    def restart(self, snapshot=-1, job_name=None, job_type=None):
        """
        Restart a new job created from an existing Lammps calculation.
        Args:
            project (pyiron.project.Project instance): Project instance at which the new job should be created
            snapshot (int): Snapshot of the calculations which would be the initial structure of the new job
            job_name (str): Job name
            job_type (str): Job type. If not specified a Lammps job type is assumed

        Returns:
            new_ham (lammps.lammps.Lammps instance): New job
        """
        new_ham = super(LammpsBase, self).restart(snapshot=snapshot, job_name=job_name, job_type=job_type)
        if new_ham.__name__ == self.__name__:
            new_ham.potential = self.potential
            if os.path.isfile(os.path.join(self.working_directory, "restart.out")):
                new_ham.read_restart_file(filename="restart.out")
                new_ham.restart_file_list.append(posixpath.join(self.working_directory, "restart.out"))
        return new_ham

    def _get_lammps_structure(self, structure=None, cutoff_radius=None):
        lmp_structure = LammpsStructure()
        lmp_structure.potential = self.input.potential
        lmp_structure.atom_type = self.input.control["atom_style"]
        if cutoff_radius is not None:
            lmp_structure.cutoff_radius = cutoff_radius
        else:
            lmp_structure.cutoff_radius = self.cutoff_radius
        lmp_structure.el_eam_lst = self.input.potential.get_element_lst()
        if structure is not None:
            lmp_structure.structure = structure
        else:
            lmp_structure.structure = self.structure
        if not set(lmp_structure.structure.get_species_symbols()).issubset(set(lmp_structure.el_eam_lst)):
            raise ValueError('The selected potentials do not support the given combination of elements.')
        return lmp_structure

    def _set_selective_dynamics(self):
        if 'selective_dynamics' in self.structure._tag_list.keys():
            if self.structure.selective_dynamics._default is None:
                self.structure.selective_dynamics._default = [True, True, True]
            sel_dyn = np.logical_not(self.structure.selective_dynamics.list())
            # Enter loop only if constraints present
            if len(np.argwhere(np.any(sel_dyn, axis=1)).flatten()) != 0:
                all_indices = np.arange(len(self.structure), dtype=int)
                constraint_xyz = np.argwhere(np.all(sel_dyn, axis=1)).flatten()
                not_constrained_xyz = np.setdiff1d(all_indices, constraint_xyz)
                # LAMMPS starts counting from 1
                constraint_xyz += 1
                ind_x = np.argwhere(sel_dyn[not_constrained_xyz, 0]).flatten()
                ind_y = np.argwhere(sel_dyn[not_constrained_xyz, 1]).flatten()
                ind_z = np.argwhere(sel_dyn[not_constrained_xyz, 2]).flatten()
                constraint_xy = not_constrained_xyz[np.intersect1d(ind_x, ind_y)] + 1
                constraint_yz = not_constrained_xyz[np.intersect1d(ind_y, ind_z)] + 1
                constraint_zx = not_constrained_xyz[np.intersect1d(ind_z, ind_x)] + 1
                constraint_x = not_constrained_xyz[np.setdiff1d(np.setdiff1d(ind_x, ind_y), ind_z)] + 1
                constraint_y = not_constrained_xyz[np.setdiff1d(np.setdiff1d(ind_y, ind_z), ind_x)] + 1
                constraint_z = not_constrained_xyz[np.setdiff1d(np.setdiff1d(ind_z, ind_x), ind_y)] + 1
                if len(constraint_xyz) > 0:
                    self.input.control['group___constraintxyz'] = 'id ' + ' '.join([str(ind) for ind in constraint_xyz])
                    self.input.control['fix___constraintxyz'] = 'constraintxyz setforce 0.0 0.0 0.0'
                    if self._generic_input['calc_mode'] == 'md':
                        self.input.control['velocity___constraintxyz'] = 'set 0.0 0.0 0.0'
                if len(constraint_xy) > 0:
                    self.input.control['group___constraintxy'] = 'id ' + ' '.join([str(ind) for ind in constraint_xy])
                    self.input.control['fix___constraintxy'] = 'constraintxy setforce 0.0 0.0 NULL'
                    if self._generic_input['calc_mode'] == 'md':
                        self.input.control['velocity___constraintxy'] = 'set 0.0 0.0 NULL'
                if len(constraint_yz) > 0:
                    self.input.control['group___constraintyz'] = 'id ' + ' '.join([str(ind) for ind in constraint_yz])
                    self.input.control['fix___constraintyz'] = 'constraintyz setforce NULL 0.0 0.0'
                    if self._generic_input['calc_mode'] == 'md':
                        self.input.control['velocity___constraintyz'] = 'set NULL 0.0 0.0'
                if len(constraint_zx) > 0:
                    self.input.control['group___constraintxz'] = 'id ' + ' '.join([str(ind) for ind in constraint_zx])
                    self.input.control['fix___constraintxz'] = 'constraintxz setforce 0.0 NULL 0.0'
                    if self._generic_input['calc_mode'] == 'md':
                        self.input.control['velocity___constraintxz'] = 'set 0.0 NULL 0.0'
                if len(constraint_x) > 0:
                    self.input.control['group___constraintx'] = 'id ' + ' '.join([str(ind) for ind in constraint_x])
                    self.input.control['fix___constraintx'] = 'constraintx setforce 0.0 NULL NULL'
                    if self._generic_input['calc_mode'] == 'md':
                        self.input.control['velocity___constraintx'] = 'set 0.0 NULL NULL'
                if len(constraint_y) > 0:
                    self.input.control['group___constrainty'] = 'id ' + ' '.join([str(ind) for ind in constraint_y])
                    self.input.control['fix___constrainty'] = 'constrainty setforce NULL 0.0 NULL'
                    if self._generic_input['calc_mode'] == 'md':
                        self.input.control['velocity___constrainty'] = 'set NULL 0.0 NULL'
                if len(constraint_z) > 0:
                    self.input.control['group___constraintz'] = 'id ' + ' '.join([str(ind) for ind in constraint_z])
                    self.input.control['fix___constraintz'] = 'constraintz setforce NULL NULL 0.0'
                    if self._generic_input['calc_mode'] == 'md':
                        self.input.control['velocity___constraintz'] = 'set NULL NULL 0.0'


class Input:
    def __init__(self):
        self.control = LammpsControl()
        self.potential = LammpsPotential()

    def to_hdf(self, hdf5):
        """

        Args:
            hdf5:

        Returns:

        """
        with hdf5.open("input") as hdf5_input:
            self.control.to_hdf(hdf5_input)
            self.potential.to_hdf(hdf5_input)

    def from_hdf(self, hdf5):
        """

        Args:
            hdf5:

        Returns:

        """
        with hdf5.open("input") as hdf5_input:
            self.control.from_hdf(hdf5_input)
            self.potential.from_hdf(hdf5_input)


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
