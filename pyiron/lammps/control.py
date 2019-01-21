# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
import numpy as np
from pyiron.base.generic.parameters import GenericParameters

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class LammpsControl(GenericParameters):
    def __init__(self, input_file_name=None, **qwargs):
        super(LammpsControl, self).__init__(input_file_name=input_file_name,
                                            table_name="control_inp",
                                            comment_char="#")
        self._init_block_dict()

    @property
    def dataset(self):
        return self._dataset

    def _init_block_dict(self):
        block_dict = OrderedDict()
        block_dict['read_restart'] = 'read_restart'
        block_dict['structure'] = ('units', 'dimension', 'boundary', 'atom_style', 'read_data', 'include')
        block_dict['potential'] = ('group', 'set', 'pair_style', 'pair_coeff', 'bond_style', 'bond_coeff',
                                   'angle_style', 'angle_coeff', 'kspace_style', 'neighbor')
        block_dict['compute'] = 'compute'
        block_dict['control'] = ('fix', 'min_style', 'min_modify', 'neigh_modify',
                                 'velocity', 'dump', 'dump_modify', 'thermo_style', 'thermo_modify', 'thermo',
                                 'timestep', 'dielectric', 'fix_modify')
        block_dict['run'] = ('run', 'minimize')
        block_dict['write_restart'] = 'write_restart'
        self._block_dict = block_dict

    def load_default(self, file_content=None):
        if file_content is None:
            file_content = ('units         metal\n'+
                            'dimension     3\n'+
                            'boundary      p p p\n'+
                            'atom_style    atomic\n'+
                            'read_data     structure.inp\n'+
                            'include       potential.inp\n'+
                            'fix           ensemble all nve\n'+
                            'variable      dumptime equal 100\n'+
                            'dump          1 all custom ${dumptime} dump.out id type xsu ysu zsu fx fy fz\n'+
                            'dump_modify   1 sort id format line "%d %d %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g"\n'+
                            'thermo_style  custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol\n'+
                            'thermo_modify format  float %20.15g\n'+
                            'thermo        100\n'+
                            'run           0\n')
        self.load_string(file_content)

    def calc_minimize(self, e_tol=0.0, f_tol=1e-8, max_iter=100000, pressure=None, n_print=100):
        max_evaluations = 10 * max_iter
        if pressure is not None:
            if type(pressure) == float or type(pressure) == int:
                pressure = pressure*np.ones(3)
            str_press = ''
            for press, str_axis in zip(pressure, [' x ', ' y ', ' z ']):
                if press is not None:
                    str_press += str_axis+str(press*1.0e4)
            if len(str_press) == 0:
                raise ValueError('Pressure values cannot be three times None')
            elif len(str_press)>1:
                str_press += ' couple none'
            self.set(fix___ensemble=r'all box/relax' + str_press)
        else:
            self.remove_keys(["fix"])
        self.set(minimize=str(e_tol) + ' ' + str(f_tol) + ' ' + str(max_iter) + " " + str(max_evaluations))
        self.remove_keys(['run', 'velocity'])
        self.modify(variable___dumptime___equal=n_print, thermo=n_print)

    def calc_static(self):
        self.set(run='0')
        self.remove_keys(['minimize', 'velocity'])

    def set_initial_velocity(self, temperature, seed=None, gaussian=False, append_value=False, zero_lin_momentum=True,
                             zero_rot_momentum=True):
        """
        Create initial velocities via velocity all create. More information can be found on LAMMPS website:
        https://lammps.sandia.gov/doc/velocity.html

        Args:
            temperature: (int or float)
            seed: (int) Seed for the initial random number generator
            gaussian: (True/False) Create velocity according to the Gaussian distribution (otherwise uniform)
            append_value: (True/False) Add the velocity values to the current velocities (probably not functional now)
            zero_lin_momentum: (True/False) Cancel the total linear momentum
            zero_rot_momentum: (True/False) Cancel the total angular momentum
        """

        if seed is None:
            seed = np.random.randint(99999)
        arg = ''
        if gaussian:
            arg = ' dist gaussian'
        if append_value:
            arg += ' sum yes'
        if not zero_lin_momentum:
            arg += ' mom no'
        if not zero_rot_momentum:
            arg += ' rot no'
        self.modify(velocity='all create ' + str(temperature) + ' ' + str(seed) + arg,
                    append_if_not_present=True)

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, time_step=1.0, n_print=100,
                delta_temp=100.0, delta_press=1000.0, seed=None, tloop=None, initial_temperature=None, langevin=False):
        """
        Set an MD calculation within LAMMPS. Nosé Hoover is used by default

        Args:
            temperature (None/float): Target temperature. If set to None, an NVE calculation is performed.
                                      It is required when the pressure is set or langevin is set
            pressure (None/float): Target pressure. If set to None, an NVE or an NVT calculation is performed.
                                   (This tag will allow for a list in the future as it is done for calc_minimize())
            n_ionic_steps (int): Number of ionic steps
            time_step (float): Step size between two steps. In fs if units==metal
            n_print (int):  Print frequency
            delta_temp (float):  Temperature damping factor (cf. https://lammps.sandia.gov/doc/fix_nh.html)
            delta_press (float): Pressure damping factor (cf. https://lammps.sandia.gov/doc/fix_nh.html)
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
        """

        if time_step is not None:
            # time_step in fs
            if self['units'] == 'metal':
                self['timestep'] = time_step * 1e-3  # lammps in ps
            elif self['units'] in ['si', 'cgs']:
                self['timestep'] = time_step * 1e-12
            elif self['units'] in ['real', 'electron']:
                self['timestep'] = time_step
            else:
                raise NotImplementedError()

        if initial_temperature is None and temperature is not None:
            initial_temperature = 2*temperature

        if seed is None:
            seed = np.random.randint(99999)
        if pressure is not None:
            pressure = float(pressure)*1.0e4  # TODO; why needed?
            if delta_press is None:
                delta_press = delta_temp*10
            if temperature is None or temperature == 0.0:
                raise ValueError('Target temperature for fix nvt/npt/nph cannot be 0.0')
            if langevin:
                fix_ensemble_str = 'all nph aniso {0} {1} {2}'.format(str(pressure), str(pressure), str(delta_press))
                self.modify(fix___langevin='all langevin {0} {1} {2} {3}'.format(str(temperature), str(temperature), str(delta_temp), str(seed)),
                            append_if_not_present=True)
            else:
                fix_ensemble_str = 'all npt temp {0} {1} {2} aniso {3} {4} {5}'.format(str(temperature), str(temperature),
                                                                                       str(delta_temp), str(pressure), str(pressure),
                                                                                       str(delta_press))
        elif temperature is not None:
            temperature = float(temperature)  # TODO; why needed?
            if temperature == 0.0:
                raise ValueError('Target temperature for fix nvt/npt/nph cannot be 0.0')
            if langevin:
                fix_ensemble_str = 'all nve'
                self.modify(fix___langevin='all langevin {0} {1} {2} {3}'.format(str(temperature), str(temperature), str(delta_temp), str(seed)),
                            append_if_not_present=True)
            else:
                fix_ensemble_str = 'all nvt temp {0} {1} {2}'.format(str(temperature), str(temperature), str(delta_temp))
        else:
            fix_ensemble_str = 'all nve'
            initial_temperature = 0
        if tloop is not None:
            fix_ensemble_str += " tloop " + str(tloop)
        self.remove_keys(["minimize"])
        self.modify(fix___ensemble=fix_ensemble_str,
                    variable=' dumptime equal {} '.format(n_print),
                    thermo=int(n_print),
                    run=int(n_ionic_steps),
                    append_if_not_present=True)
        if initial_temperature>0:
            self.set_initial_velocity(initial_temperature, gaussian=True)
