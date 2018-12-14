# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function

import os
import posixpath
import numpy as np
import pandas as pd
import warnings

from pyiron.lammps.potential import LammpsPotentialFile
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.atomistics.job.interactive import GenericInteractive
from pyiron.base.settings.generic import Settings
from pyiron.lammps.interface import InteractiveLammpsInterface
from pyiron.lammps.control import LammpsControl
from pyiron.lammps.potential import LammpsPotential

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

try:
    from lammps import PyLammps
except ImportError:
    pass

s = Settings()


class Lammps(GenericInteractive, AtomisticGenericJob):
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
        super(Lammps, self).__init__(project, job_name)
        self.__name__ = "Lammps"
        self.__version__ = None  # Reset the version number to the executable is set automatically
        self._executable_activate()
        self.input = Input()
        self._cutoff_radius = None
        self._is_continuation = None
        self._interface = InteractiveLammpsInterface(job=self)

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
        if isinstance(potential_filename, str):
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

    def validate_ready_to_run(self):
        """
        
        Returns:

        """
        super(Lammps, self).validate_ready_to_run()
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
        return self.get_structure()

    def get_structure(self, iteration_step=-1):
        """
        
        Args:
            iteration_step:

        Returns:

        """
        if not (self.structure is not None):
            raise AssertionError()
        snapshot = self.structure.copy()
        snapshot.cell = self.get("output/generic/cells")[iteration_step]
        snapshot.positions = self.get("output/generic/positions")[iteration_step]
        return snapshot

    def view_potentials(self):
        """

        Returns:

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

        Returns:
            list:
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
        self._interface.write_input(job=self)

    def collect_output(self):
        """
        
        Returns:

        """
        if self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal:
            pass
        else:
            self._interface.collect_output(job=self)

    def convergence_check(self):
        if self._generic_input['calc_mode'] == 'minimize':
            if self._generic_input['max_iter']+1 <= len(self['output/generic/energy_tot']) or \
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

    def calc_minimize(self, e_tol=0.0, f_tol=1e-8, max_iter=100000, pressure=None, n_print=100):
        """
        
        Args:
            e_tol:
            f_tol:
            max_iter:
            pressure:
            n_print:

        Returns:

        """
        if self.server.run_mode.interactive_non_modal:
            warnings.warn('calc_minimize() is not implemented for the non modal interactive mode use calc_static()!')
        super(Lammps, self).calc_minimize(e_tol=e_tol, f_tol=f_tol, max_iter=max_iter, pressure=pressure, n_print=n_print)
        self.input.control.calc_minimize(e_tol=e_tol, f_tol=f_tol, max_iter=max_iter, pressure=pressure, n_print=n_print)

    def calc_static(self):
        """
        
        Returns:

        """
        super(Lammps, self).calc_static()
        self.input.control.calc_static()

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, dt=None, time_step=None, n_print=100,
                delta_temp=100.0, delta_press=None, seed=None, tloop=None, rescale_velocity=True, langevin=False):
        """
        
        Args:
            temperature:
            pressure:
            n_ionic_steps:
            dt:
            time_step:
            n_print:
            delta_temp:
            delta_press:
            seed:
            tloop:
            rescale_velocity:

        """
        if self.server.run_mode.interactive_non_modal:
            warnings.warn('calc_md() is not implemented for the non modal interactive mode use calc_static()!')
        if dt is not None:
            time_step = dt
        super(Lammps, self).calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps, 
                                    time_step=time_step, n_print=n_print, delta_temp=delta_temp, delta_press=delta_press,
                                    seed=seed, tloop=tloop, rescale_velocity=rescale_velocity, langevin=langevin)
        self.input.control.calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps,
                                   time_step=time_step, n_print=n_print, delta_temp=delta_temp, delta_press=delta_press,
                                   seed=seed, tloop=tloop, rescale_velocity=rescale_velocity, langevin=langevin)

    # define hdf5 input and output
    def to_hdf(self, hdf=None, group_name=None):
        """
        
        Args:
            hdf:
            group_name:

        Returns:

        """
        super(Lammps, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()
        self.input.to_hdf(self._hdf5)

    def from_hdf(self, hdf=None, group_name=None):  # TODO: group_name should be removed
        """
        
        Args:
            hdf:
            group_name:

        Returns:

        """
        super(Lammps, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()
        self.input.from_hdf(self._hdf5)

    def write_restart_file(self, filename="restart.out"):
        """
        
        Args:
            filename:

        Returns:

        """
        self.input.control.modify(write_restart=filename, append_if_not_present=True)

    def read_restart_file(self, filename="restart.out"):
        """
        
        Args:
            filename:

        Returns:

        """
        self._is_continuation = True
        self.input.control.set(read_restart=filename)
        self.input.control.remove_keys(['dimension', 'read_data', 'boundary', 'atom_style', 'velocity'])

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
        return super(Lammps, self).restart(snapshot=snapshot, job_name=job_name, job_type=job_type)

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
        new_ham = super(Lammps, self).restart(snapshot=snapshot, job_name=job_name, job_type=job_type)
        if new_ham.__name__ == self.__name__:
            new_ham.potential = self.potential
            if os.path.isfile(os.path.join(self.working_directory, "restart.out")):
                new_ham.read_restart_file(filename="restart.out")
                new_ham.restart_file_list.append(posixpath.join(self.working_directory, "restart.out"))
        return new_ham

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


class LammpsInt(Lammps):
    def __init__(self, project, job_name):
        warnings.warn('Please use Lammps instead of LammpsInt')
        super(LammpsInt, self).__init__(project=project, job_name=job_name)


class LammpsInt2(LammpsInt):
    def __init__(self, project, job_name):
        warnings.warn('Please use Lammps instead of LammpsInt2')
        super(LammpsInt2, self).__init__(project=project, job_name=job_name)
