# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
import hashlib
import numpy as np
import warnings
from pyiron.base.generic.parameters import GenericParameters
import scipy.constants as spc

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
                                 'timestep', 'dielectric', 'fix_modify', 'reset_timestep')
        block_dict['run'] = ('run', 'minimize')
        block_dict['write_restart'] = 'write_restart'
        self._block_dict = block_dict

    def load_default(self, file_content=None):
        if file_content is None:
            file_content = ('units               metal\n'+
                            'dimension           3\n'+
                            'boundary            p p p\n'+
                            'atom_style          atomic\n'+
                            'read_data           structure.inp\n'+
                            'include             potential.inp\n'+
                            'fix___ensemble      all nve\n'+
                            'variable___dumptime equal 100\n'+
                            'dump___1            all custom ${dumptime} dump.out id type xsu ysu zsu fx fy fz\n'+
                            'dump_modify___1     sort id format line "%d %d %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g"\n'+
                            'thermo_style        custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol\n'+
                            'thermo_modify       format float %20.15g\n'+
                            'thermo              100\n'+
                            'run                 0\n')
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
            self.remove_keys(["fix___nve"])
        self.set(minimize=str(e_tol) + ' ' + str(f_tol) + ' ' + str(max_iter) + " " + str(max_evaluations))
        self.remove_keys(['run', 'velocity'])
        self.modify(variable___dumptime='equal '+str(n_print), thermo=n_print)

    def calc_static(self):
        self.set(run='0')
        self.remove_keys(['minimize', 'velocity'])

    def set_initial_velocity(self, temperature, seed=None, gaussian=False, append_value=False, zero_lin_momentum=True,
                             zero_rot_momentum=True, job_name=''):
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
            job_name: (str) job name to generate seed
        """

        if seed is None:
            seed = self.generate_seed_from_job(job_name=job_name, seed=1)
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

    @staticmethod
    def generate_seed_from_job(job_name='', seed=0):
        """
        Generate a unique seed from the job name.

        Args:
            job_name (str): job_name of the current job to generate the seed
            seed (int): integer to access part of the seed

        Returns:
            int: random seed generated based on the hash
        """
        return int(str(int(hashlib.sha256(job_name.encode()).hexdigest(), 16))[5 * seed:5 * seed + 5])

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, time_step=1.0, n_print=100,
                temperature_damping_timescale=100.0, pressure_damping_timescale=1000.0, seed=None, tloop=None,
                initial_temperature=None, langevin=False, delta_temp=None, delta_press=None, job_name=''):
        """
        Set an MD calculation within LAMMPS. Nosé Hoover is used by default.

        Args:
            temperature (None/float): Target temperature. If set to None, an NVE calculation is performed.
                                      It is required when the pressure is set or langevin is set
            pressure (None/float/numpy.ndarray/list): Target pressure. If set to None, an NVE or an NVT calculation is
                performed. A length-3 list or array may be given to specify x-, y- and z-components individually. In
                this case, floats and `None` may be mixed to allow relaxation only in particular directions.
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
            job_name (str): Job name of the job to generate a unique random seed.
        """
        # Conversion factors for transfroming pyiron units to Lammps units
        fs_to_ps = spc.femto / spc.pico
        fs_to_s = spc.femto / 1.
        GPa_to_bar = spc.giga * 1. / spc.bar
        GPa_to_Pa = spc.giga
        GPa_to_barye = spc.giga * 1. / (1.e-6 * spc.bar)  # Lammps is in "barye"
        GPa_to_atm = spc.giga * 1. / spc.atm
        lammps_unit_conversions = {'metal': {'time': fs_to_ps,
                                             'pressure': GPa_to_bar},
                                   'si': {'time': fs_to_s,
                                          'pressure': GPa_to_Pa},
                                   'cgs': {'time': fs_to_s,
                                           'pressure': GPa_to_barye},
                                   'real': {'time': 1,
                                            'pressure': GPa_to_atm},
                                   'electron': {'time': 1,
                                                'pressure': GPa_to_Pa}}
        time_units = lammps_unit_conversions[self['units']]['time']
        pressure_units = lammps_unit_conversions[self['units']]['pressure']
        # No need for temperature conversion; pyiron and all available Lammps units are both in Kelvin
        # (well, except unitless Lennard-Jones units...)
        if self['units'] == 'lj':
            raise NotImplementedError

        # Transform time
        if time_step is not None:
            try:
                self['timestep'] = time_step * time_units
            except KeyError:
                raise NotImplementedError()

        # Transform thermostat strength
        if delta_temp is not None:
            warnings.warn("WARNING: `delta_temp` is deprecated, please use `temperature_damping_timescale`.")
            temperature_damping_timescale = delta_temp
        else:
            temperature_damping_timescale *= time_units

        # Transform barostat strength
        if delta_press is not None:
            warnings.warn("WARNING: `delta_press` is deprecated, please use `pressure_damping_timescale`.")
            pressure_damping_timescale = delta_press
        else:
            pressure_damping_timescale *= time_units

        # Apply initial overheating (default uses the theorem of equipartition of energy between KE and PE)
        if initial_temperature is None and temperature is not None:
            initial_temperature = 2 * temperature

        if seed is None:
            seed = self.generate_seed_from_job(job_name=job_name)

        # Set thermodynamic ensemble
        if pressure is not None:  # NPT
            if not hasattr(pressure, '__len__'):
                pressure = pressure * np.ones(3)
            else:
                pressure = np.array(pressure)

            if sum(pressure != None) == 0:
                raise ValueError('Pressure cannot be three times None')

            if len(pressure) != 3:
                raise ValueError('Pressure must be a float or a 3d vector')

            if temperature is None or temperature == 0.0:
                raise ValueError('Target temperature for fix nvt/npt/nph cannot be 0')

            pressure[pressure != None] *= pressure_units

            pressure_string = ''
            for coord, value in zip(['x', 'y', 'z'], pressure):
                if value is not None:
                    pressure_string += ' {0} {1} {1} {2}'.format(coord, str(value), str(pressure_damping_timescale))

            if langevin:  # NPT(Langevin)
                fix_ensemble_str = 'all nph' + pressure_string
                self.modify(fix___langevin='all langevin {0} {1} {2} {3} zero yes'.format(str(temperature),
                                                                                          str(temperature),
                                                                                          str(temperature_damping_timescale),
                                                                                          str(seed)),
                            append_if_not_present=True)
            else:  # NPT(Nose-Hoover)
                fix_ensemble_str = 'all npt temp {0} {1} {2}'.format(str(temperature),
                                                                     str(temperature),
                                                                     str(temperature_damping_timescale))
                fix_ensemble_str += pressure_string
        elif temperature is not None:  # NVT
            if temperature == 0.0:
                raise ValueError('Target temperature for fix nvt/npt/nph cannot be 0.0')

            if langevin:  # NVT(Langevin)
                fix_ensemble_str = 'all nve'
                self.modify(fix___langevin='all langevin {0} {1} {2} {3} zero yes'.format(str(temperature),
                                                                                          str(temperature),
                                                                                          str(temperature_damping_timescale),
                                                                                          str(seed)),
                            append_if_not_present=True)
            else:  # NVT(Nose-Hoover)
                fix_ensemble_str = 'all nvt temp {0} {1} {2}'.format(str(temperature),
                                                                     str(temperature),
                                                                     str(temperature_damping_timescale))
        else:  # NVE
            if langevin:
                warnings.warn("Temperature not set; Langevin ignored.")
            fix_ensemble_str = 'all nve'
            initial_temperature = 0

        if tloop is not None:
            fix_ensemble_str += " tloop " + str(tloop)

        self.remove_keys(["minimize"])
        self.modify(fix___ensemble=fix_ensemble_str,
                    variable___dumptime=' equal {} '.format(n_print),
                    thermo=int(n_print),
                    run=int(n_ionic_steps),
                    append_if_not_present=True)

        if initial_temperature > 0:
            self.set_initial_velocity(initial_temperature, gaussian=True, job_name=job_name)
