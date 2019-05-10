# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function

from copy import copy
import numpy as np

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class ThermoBulk(object):
    """
    Class should provide all tools to compute bulk thermodynamic quantities. Central quantity is the Free Energy F(V,T).
    ToDo: Make it a (light weight) pyiron object (introduce a new tool rather than job object). 
       
    Args:
        project: 
        name: 
            
    """
    eV_to_J_per_mol = 1.60217662e-19 * 6.022e23
    kB = 1 / 8.6173303e-5

    def __init__(self, project=None, name=None):
        # only for compatibility with pyiron objects
        self._project = project
        self._name = name

        self._volumes = None
        self._temperatures = None
        self._energies = None
        self._entropy = None
        self._pressure = None
        self._num_atoms = None

        self._fit_order = 3

    def copy(self):
        """
        
        Returns:

        """
        cls = self.__class__
        result = cls.__new__(cls)
        result.__init__()
        result.__dict__['_volumes'] = copy(self._volumes)
        result.__dict__['_temperatures'] = copy(self._temperatures)
        result.__dict__['_energies'] = copy(self._energies)
        result.__dict__['_fit_order'] = self._fit_order
        return result

    def _reset_energy(self):
        """
        
        Returns:

        """
        if self._volumes is not None:
            if self._temperatures is not None:
                self._energies = np.zeros((len(self._temperatures), len(self._volumes)))
                # self.energies = 0

    @property
    def num_atoms(self):
        """
        
        Returns:

        """
        if self._num_atoms is None:
            return 1   # normalize per cell if number of atoms unknown
        return self._num_atoms

    @num_atoms.setter
    def num_atoms(self, num):
        """
        
        Args:
            num: 

        Returns:

        """
        self._num_atoms = num

    @property
    def _coeff(self):
        """
        
        Returns:

        """
        return np.polyfit(self._volumes, self._energies.T, deg=self._fit_order)

    @property
    def temperatures(self):
        """
        
        Returns:

        """
        return self._temperatures

    @property
    def _d_temp(self):
        """
        
        Returns:

        """
        return self.temperatures[1] - self.temperatures[0]

    @property
    def _d_vol(self):
        """
        
        Returns:

        """
        return self.volumes[1] - self.volumes[0]

    @temperatures.setter
    def temperatures(self, temp_lst):
        """
        
        Args:
            temp_lst: 

        Returns:

        """
        if not hasattr(temp_lst, '__len__'):
            raise ValueError('Requires list as input parameter')
        len_temp = -1
        if self._temperatures is not None:
            len_temp = len(self._temperatures)
        self._temperatures = np.array(temp_lst)
        if len(temp_lst) != len_temp:
            self._reset_energy()

    @property
    def volumes(self):
        """
        
        Returns:

        """
        return self._volumes

    @volumes.setter
    def volumes(self, volume_lst):
        """
        
        Args:
            volume_lst: 

        Returns:

        """
        if not hasattr(volume_lst, '__len__'):
            raise ValueError('Requires list as input parameter')
        len_vol = -1
        if self._volumes is not None:
            len_vol = len(self._volumes)
        self._volumes = np.array(volume_lst)
        if len(volume_lst) != len_vol:
            self._reset_energy()

    @property
    def entropy(self):
        """
        
        Returns:

        """
        if self._entropy is None:
            self._compute_thermo()
        return self._entropy

    @property
    def pressure(self):
        """
        
        Returns:

        """
        if self._pressure is None:
            self._compute_thermo()
        return self._pressure

    @property
    def energies(self):
        """
        
        Returns:

        """
        return self._energies

    @energies.setter
    def energies(self, erg_lst):
        """
        
        Args:
            erg_lst: 

        Returns:

        """
        if np.ndim(erg_lst) == 2:
            self._energies = erg_lst
        elif np.ndim(erg_lst) == 1:
            if len(erg_lst) == len(self.volumes):
                self._energies = np.tile(erg_lst, (len(self.temperatures), 1))
            else:
                raise ValueError()
        else:
            self._energies = np.ones((len(self.volumes), len(self.temperatures))) * erg_lst

    def set_temperatures(self, temperature_min=0, temperature_max=1500, temperature_steps=50):
        """
        
        Args:
            temperature_min: 
            temperature_max: 
            temperature_steps: 

        Returns:

        """
        self.temperatures = np.linspace(temperature_min, temperature_max, temperature_steps)

    def set_volumes(self, volume_min, volume_max=None, volume_steps=10):
        """
        
        Args:
            volume_min: 
            volume_max: 
            volume_steps: 

        Returns:

        """
        if volume_max is None:
            volume_max = 1.1 * volume_min
        self.volumes = np.linspace(volume_min, volume_max, volume_steps)

    def meshgrid(self):
        """
        
        Returns:

        """
        return np.meshgrid(self.volumes, self.temperatures)

    def get_minimum_energy_path(self, pressure=None):
        """
        
        Args:
            pressure: 

        Returns:

        """
        if pressure is not None:
            raise NotImplemented()
        v_min_lst = []
        for c in self._coeff.T:
            v_min = np.roots(np.polyder(c, 1))
            p_der2 = np.polyder(c, 2)
            p_val2 = np.polyval(p_der2, v_min)
            v_m_lst = v_min[p_val2 > 0]
            if len(v_m_lst) > 0:
                v_min_lst.append(v_m_lst[0])
            else:
                v_min_lst.append(np.nan)
        return np.array(v_min_lst)

    def get_free_energy(self, vol, pressure=None):
        """
        
        Args:
            vol: 
            pressure: 

        Returns:

        """
        if not pressure:
            return np.polyval(self._coeff, vol)
        else:
            raise NotImplementedError()

    def interpolate_volume(self, volumes, fit_order=None):
        """
        
        Args:
            volumes: 
            fit_order: 

        Returns:

        """
        if fit_order is not None:
            self._fit_order = fit_order
        new = self.copy()
        new.volumes = volumes
        new.energies = np.array([np.polyval(self._coeff, v) for v in volumes]).T
        return new

    def _compute_thermo(self):
        """
        
        Returns:

        """
        self._entropy, self._pressure = np.gradient(-self.energies, self._d_temp, self._d_vol)

    def get_free_energy_p(self):
        """
        
        Returns:

        """
        coeff = np.polyfit(self._volumes, self.energies.T, deg=self._fit_order)
        return np.polyval(coeff, self.get_minimum_energy_path())

    def get_entropy_p(self):
        """
        
        Returns:

        """
        s_coeff = np.polyfit(self._volumes, self.entropy.T, deg=self._fit_order)
        return np.polyval(s_coeff, self.get_minimum_energy_path())

    def get_entropy_v(self):
        """
        
        Returns:

        """
        eq_volume = self.volumes[0]
        s_coeff = np.polyfit(self.volumes, self.entropy.T, deg=self._fit_order)
        const_v = eq_volume * np.ones(len(s_coeff.T))
        return np.polyval(s_coeff, const_v)

    def plot_free_energy(self):
        """
        
        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        plt.plot(self.temperatures, self.get_free_energy_p()/self.num_atoms)
        plt.xlabel('Temperature [K]')
        plt.ylabel('Free energy [eV]')

    def plot_entropy(self):
        """
        
        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        plt.plot(self.temperatures, self.eV_to_J_per_mol / self.num_atoms * self.get_entropy_p(), label='S$_p$')
        plt.plot(self.temperatures, self.eV_to_J_per_mol / self.num_atoms * self.get_entropy_v(), label='S$_V$')
        plt.legend()
        plt.xlabel('Temperature [K]')
        plt.ylabel('Entropy [J K$^{-1}$ mol-atoms$^{-1}$]')

    def plot_heat_capacity(self, to_kB=True):
        """
        
        Args:
            to_kB: 

        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        if to_kB:
            units = self.kB / self.num_atoms
            plt.ylabel('Heat capacity [kB]')
        else:
            units = self.eV_to_J_per_mol
            plt.ylabel('Heat capacity [J K$^{-1}$ mol-atoms$^{-1}$]')
        temps = self.temperatures[:-2]
        c_p = temps * np.gradient(self.get_entropy_p(), self._d_temp)[:-2]
        c_v = temps * np.gradient(self.get_entropy_v(), self._d_temp)[:-2]
        plt.plot(temps, units * c_p, label='c$_p$')
        plt.plot(temps, units * c_v, label='c$_v$')
        plt.legend(loc='lower right')
        plt.xlabel('Temperature [K]')

    def contour_pressure(self):
        """
        
        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        x, y = self.meshgrid()
        p_coeff = np.polyfit(self.volumes, self.pressure.T, deg=self._fit_order)
        p_grid = np.array([np.polyval(p_coeff, v) for v in self._volumes]).T
        plt.contourf(x, y, p_grid)
        plt.plot(self.get_minimum_energy_path(), self.temperatures)
        plt.xlabel('Volume [$\AA^3$]')
        plt.ylabel('Temperature [K]')

    def contour_entropy(self):
        """
        
        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        s_coeff = np.polyfit(self.volumes, self.entropy.T, deg=self._fit_order)
        s_grid = np.array([np.polyval(s_coeff, v) for v in self.volumes]).T
        x, y = self.meshgrid()
        plt.contourf(x, y, s_grid)
        plt.plot(self.get_minimum_energy_path(), self.temperatures)
        plt.xlabel('Volume [$\AA^3$]')
        plt.ylabel('Temperature [K]')

    def plot_contourf(self, ax=None, show_min_erg_path=False):
        """
        
        Args:
            ax: 
            show_min_erg_path: 

        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        x, y = self.meshgrid()
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax.contourf(x, y, self.energies)
        if show_min_erg_path:
            plt.plot(self.get_minimum_energy_path(), self.temperatures, 'w--')
        plt.xlabel('Volume [$\AA^3$]')
        plt.ylabel('Temperature [K]')
        return ax

    def plot_min_energy_path(self, *args, ax=None, **qwargs):
        """
        
        Args:
            *args: 
            ax: 
            **qwargs: 

        Returns:

        """
        try:
            import pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 1)
            ax.xlabel('Volume [$\AA^3$]')
            ax.ylabel('Temperature [K]')
        ax.plot(self.get_minimum_energy_path(), self.temperatures, *args, **qwargs)
        return ax
