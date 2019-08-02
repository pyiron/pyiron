# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function

import numpy as np
import scipy.integrate
import scipy.optimize as spy
import scipy.constants
import warnings
from pyiron.atomistics.master.parallel import AtomisticParallelMaster
from pyiron.base.master.parallel import JobGenerator

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"



eV_div_A3_to_GPa = 1e21 / scipy.constants.physical_constants['joule-electron volt relationship'][0]


def _debye_kernel(xi):
    return xi ** 3 / (np.exp(xi) - 1)


def debye_integral(x):
    return scipy.integrate.quad(_debye_kernel, 0, x)[0]


def debye_function(x):
    if hasattr(x, '__len__'):
        return np.array([3 / xx ** 3 * debye_integral(xx) for xx in x])
    return 3 / x ** 3 * debye_integral(x)


# https://gitlab.com/ase/ase/blob/master/ase/eos.py
def birchmurnaghan_energy(V, E0, B0, BP, V0):
    'BirchMurnaghan equation from PRB 70, 224107'
    eta = (V0 / V) ** (1 / 3)
    return E0 + 9 * B0 * V0 / 16 * (eta ** 2 - 1) ** 2 * (6 + BP * (eta ** 2 - 1) - 4 * eta ** 2)


def vinet_energy(V, E0, B0, BP, V0):
    'Vinet equation from PRB 70, 224107'
    eta = (V / V0) ** (1 / 3)
    return (E0 + 2 * B0 * V0 / (BP - 1) ** 2 * (
            2 - (5 + 3 * BP * (eta - 1) - 3 * eta) * np.exp(-3 * (BP - 1) * (eta - 1) / 2)))


def murnaghan(V, E0, B0, BP, V0):
    'From PRB 28,5480 (1983'
    E = E0 + B0 * V / BP * (((V0 / V) ** BP) / (BP - 1) + 1) - V0 * B0 / (BP - 1)
    return E


def birch(V, E0, B0, BP, V0):
    """
    From Intermetallic compounds: Principles and Practice, Vol. I: Principles
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
    paper downloaded from Web

    case where n=0
    """
    E = (E0 +
         9 / 8 * B0 * V0 * ((V0 / V) ** (2 / 3) - 1) ** 2 +
         9 / 16 * B0 * V0 * (BP - 4) * ((V0 / V) ** (2 / 3) - 1) ** 3)
    return E


def pouriertarantola(V, E0, B0, BP, V0):
    'Pourier-Tarantola equation from PRB 70, 224107'
    eta = (V / V0) ** (1 / 3)
    squiggle = -3 * np.log(eta)

    E = E0 + B0 * V0 * squiggle ** 2 / 6 * (3 + squiggle * (BP - 2))
    return E


def fitfunction(parameters, vol, fittype='vinet'):
    """
    Fit the energy volume curve

    Args:
        parameters (list): [E0, B0, BP, V0] list of fit parameters
        vol (float/numpy.dnarray): single volume or a vector of volumes as numpy array
        fittype (str): on of the following ['birch', 'birchmurnaghan', 'murnaghan', 'pouriertarantola', 'vinet']

    Returns:
        (float/numpy.dnarray): single energy as float or a vector of energies as numpy array
    """
    [E0, b0, bp, V0] = parameters
    # Unit correction
    B0 = b0 / eV_div_A3_to_GPa
    BP = bp
    V = vol
    if fittype.lower() == 'birchmurnaghan':
        return birchmurnaghan_energy(V, E0, B0, BP, V0)
    elif fittype.lower() == 'vinet':
        return vinet_energy(V, E0, B0, BP, V0)
    elif fittype.lower() == 'murnaghan':
        return murnaghan(V, E0, B0, BP, V0)
    elif fittype.lower() == 'pouriertarantola':
        return pouriertarantola(V, E0, B0, BP, V0)
    elif fittype.lower() == 'birch':
        return birch(V, E0, B0, BP, V0)
    else:
        raise ValueError


def fit_leastsq(p0, datax, datay, fittype='vinet'):
    """
    Least square fit

    Args:
        p0 (list): [E0, B0, BP, V0] list of fit parameters
        datax (float/numpy.dnarray): volumes to fit
        datay (float/numpy.dnarray): energies corresponding to the volumes
        fittype (str): on of the following ['birch', 'birchmurnaghan', 'murnaghan', 'pouriertarantola', 'vinet']

    Returns:
        list: [E0, B0, BP, V0], [E0_err, B0_err, BP_err, V0_err]
    """
    # http://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i

    errfunc = lambda p, x, y, fittype: fitfunction(p, x, fittype) - y

    pfit, pcov, infodict, errmsg, success = \
        spy.leastsq(errfunc, p0, args=(datax, datay, fittype), full_output=1, epsfcn=0.0001)

    if (len(datay) > len(p0)) and pcov is not None:
        s_sq = (errfunc(pfit, datax, datay, fittype)**2).sum()/(len(datay)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = np.inf

    error = []
    for i in range(len(pfit)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
    pfit_leastsq = pfit
    perr_leastsq = np.array(error)
    return pfit_leastsq, perr_leastsq


class DebyeModel(object):
    """
    Calculate Thermodynamic Properties based on the Murnaghan output
    """
    def __init__(self, murnaghan, num_steps=50):
        self._murnaghan = murnaghan

        # self._atoms_per_cell = len(murnaghan.structure)
        self._v_min = None
        self._v_max = None
        self._num_steps = None

        self._volume = None
        self._init_volume()

        self.num_steps = num_steps
        self._fit_volume = None
        self._debye_T = None

    def _init_volume(self):
        vol = self._murnaghan['output/volume']
        self._v_min, self._v_max = np.min(vol), np.max(vol)

    def _set_volume(self):
        if self._v_min and self._v_max and self._num_steps:
            self._volume = np.linspace(self._v_min, self._v_max, self._num_steps)
            self._reset()
            # print ('set_volume: ', self._num_steps)

    @property
    def num_steps(self):
        return self._num_steps

    @num_steps.setter
    def num_steps(self, val):
        self._num_steps = val
        self._set_volume()

    @property
    def volume(self):
        if self._volume is None:
            self._init_volume()
            self._set_volume()
        return self._volume

    @volume.setter
    def volume(self, volume_lst):
        self._volume = volume_lst
        self._v_min = np.min(volume_lst)
        self._v_max = np.max(volume_lst)
        self._reset()

    def _reset(self):
        self._debye_T = None

    def polynomial(self, poly_fit=None, volumes=None):
        if poly_fit is None:
            self._murnaghan.fit_polynomial()  # TODO: include polyfit in output
            poly_fit = self._murnaghan.fit_dict['poly_fit']
        p_fit = np.poly1d(poly_fit)
        if volumes is None:
            return p_fit(self.volume)
        return p_fit(volumes)

    @property
    def debye_temperature(self):
        if self._debye_T is not None:
            return self._debye_T

        GPaTokBar = 10
        Ang3_to_Bohr3 = scipy.constants.angstrom**3/scipy.constants.physical_constants['Bohr radius'][0]**3
        convert = 67.48  # conversion factor, Moruzzi Eq. (4)
        empirical = 0.617  # empirical factor, Moruzzi Eq. (6)
        gamma_low, gamma_high = 1, 2/3  # low/high T gamma

        out = self._murnaghan['output']
        V0 = out['equilibrium_volume']
        B0 = out['equilibrium_bulk_modulus']
        Bp = out['equilibrium_b_prime']

        vol = self.volume

        mass = set(self._murnaghan.structure.get_masses())
        if len(mass) > 1:
            raise NotImplementedError('Debye temperature only for single species systems!')
        mass = list(mass)[0]

        r0 = (3 * V0 * Ang3_to_Bohr3 / (4 * np.pi)) ** (1. / 3.)
        debye_zero = empirical * convert * np.sqrt(r0 * B0 * GPaTokBar / mass)
        # print('r0, B0, Bp, mass, V0', r0, B0, Bp, mass, V0)
        # print('gamma_low, gamma_high: ', gamma_low, gamma_high)
        # print('debye_zero, V0: ', debye_zero, V0)
        if vol is None:
            print('WARNING: vol: ', vol)

        debye_low = debye_zero * (V0 / vol) ** (-gamma_low + 0.5 * (1 + Bp))
        debye_high = debye_zero * (V0 / vol) ** (-gamma_high + 0.5 * (1 + Bp))

        self._debye_T = (debye_low, debye_high)
        return self._debye_T

    def energy_vib(self, T, debye_T=None, low_T_limit=True):
        kB = 0.086173422 / 1000  # eV/K
        if debye_T is None:
            if low_T_limit:
                debye_T = self.debye_temperature[0]  # low
            else:
                debye_T = self.debye_temperature[1]  # high
        if hasattr(debye_T, '__len__'):
            val = [9. / 8. * kB * d_T + T * kB * (3 * np.log(1 - np.exp(-d_T / T)) - debye_function(d_T / T))
                   for d_T in debye_T]
            val = np.array(val)
        else:
            val = 9. / 8. * kB * debye_T + T * kB * (3 * np.log(1 - np.exp(-debye_T / T)) - debye_function(debye_T / T))
        atoms_per_cell = len(self._murnaghan.structure)
        return atoms_per_cell * val


class MurnaghanJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        parameter_lst = []
        for strain in np.linspace(1 - self._job.input['vol_range'],
                                  1 + self._job.input['vol_range'],
                                  self._job.input['num_points']):
            basis = self._job.ref_job.structure.copy()
            basis.set_cell(basis.cell * strain ** (1. / 3.), scale_atoms=True)
            parameter_lst.append([np.round(strain, 7), basis])
        return parameter_lst

    @staticmethod
    def job_name(parameter):
        return "strain_" + str(parameter[0]).replace('.', '_')

    def modify_job(self, job, parameter):
        job.structure = parameter[1]
        return job


class EnergyVolumeFit(object):
    """
    Fit energy volume curves

    Args:
        volume_lst (list/numpy.dnarray): vector of volumes
        energy_lst (list/numpy.dnarray): vector of energies

    Attributes:

        .. attribute:: volume_lst

            vector of volumes

        .. attribute:: energy_lst

            vector of energies

        .. attribute:: fit_dict

            dictionary of fit parameters
    """
    def __init__(self, volume_lst=None, energy_lst=None):
        self._volume_lst = volume_lst
        self._energy_lst = energy_lst
        self._fit_dict = None

    @property
    def volume_lst(self):
        return self._volume_lst

    @volume_lst.setter
    def volume_lst(self, vol_lst):
        self._volume_lst = vol_lst

    @property
    def energy_lst(self):
        return self._energy_lst

    @energy_lst.setter
    def energy_lst(self, eng_lst):
        self._energy_lst = eng_lst

    @property
    def fit_dict(self):
        return self._fit_dict

    def _get_volume_and_energy_lst(self, volume_lst=None, energy_lst=None):
        """
        Internal function to get the vector of volumes and the vector of energies

        Args:
            volume_lst (list/numpy.dnarray/None): vector of volumes
            energy_lst (list/numpy.dnarray/None): vector of energies

        Returns:
            list: vector of volumes and vector of energies
        """
        if volume_lst is None:
            if self._volume_lst is None:
                raise ValueError('Volume list not set.')
            volume_lst = self._volume_lst
        if energy_lst is None:
            if self._energy_lst is None:
                raise ValueError('Volume list not set.')
            energy_lst = self._energy_lst
        return volume_lst, energy_lst

    def fit_eos_general_intern(self, fittype='birchmurnaghan'):
        self._fit_dict = self.fit_eos_general(volume_lst=self._volume_lst, energy_lst=self._energy_lst, fittype=fittype)

    def fit_eos_general(self, volume_lst=None, energy_lst=None, fittype='birchmurnaghan'):
        """
        Fit on of the equations of state

        Args:
            volume_lst (list/numpy.dnarray/None): vector of volumes
            energy_lst (list/numpy.dnarray/None): vector of energies
            fittype (str): on of the following ['birch', 'birchmurnaghan', 'murnaghan', 'pouriertarantola', 'vinet']

        Returns:
            dict: dictionary with fit results
        """
        volume_lst, energy_lst = self._get_volume_and_energy_lst(volume_lst=volume_lst, energy_lst=energy_lst)
        fit_dict = {}
        pfit_leastsq, perr_leastsq = self._fit_leastsq(volume_lst=volume_lst, energy_lst=energy_lst, fittype=fittype)
        fit_dict["fit_type"] = fittype
        fit_dict["volume_eq"] = pfit_leastsq[3]
        fit_dict["energy_eq"] = pfit_leastsq[0]
        fit_dict["bulkmodul_eq"] = pfit_leastsq[1]
        fit_dict["b_prime_eq"] = pfit_leastsq[2]
        fit_dict["least_square_error"] = perr_leastsq  # [e0, b0, bP, v0]

        return fit_dict

    def fit_polynomial(self, volume_lst=None, energy_lst=None, fit_order=3):
        """
        Fit a polynomial

        Args:
            volume_lst (list/numpy.dnarray/None): vector of volumes
            energy_lst (list/numpy.dnarray/None): vector of energies
            fit_order (int): Degree of the polynomial

        Returns:
            dict: dictionary with fit results
        """
        volume_lst, energy_lst = self._get_volume_and_energy_lst(volume_lst=volume_lst, energy_lst=energy_lst)
        fit_dict = {}

        # compute a polynomial fit
        z = np.polyfit(volume_lst, energy_lst, fit_order)
        p_fit = np.poly1d(z)
        fit_dict["poly_fit"] = z

        # get equilibrium lattice constant
        # search for the local minimum with the lowest energy
        p_deriv_1 = np.polyder(p_fit, 1)
        roots = np.roots(p_deriv_1)

        # volume_eq_lst = np.array([np.real(r) for r in roots if np.abs(np.imag(r)) < 1e-10])
        volume_eq_lst = np.array([np.real(r) for r in roots if (abs(np.imag(r)) < 1e-10 and
                                                                r>=min(volume_lst) and
                                                                r<=max(volume_lst))])

        e_eq_lst = p_fit(volume_eq_lst)
        arg = np.argsort(e_eq_lst)
        # print ("v_eq:", arg, e_eq_lst)
        if len(e_eq_lst) == 0:
            return None
        e_eq = e_eq_lst[arg][0]
        volume_eq = volume_eq_lst[arg][0]

        # get bulk modulus at equ. lattice const.
        p_2deriv = np.polyder(p_fit, 2)
        p_3deriv = np.polyder(p_fit, 3)
        a2 = p_2deriv(volume_eq)
        a3 = p_3deriv(volume_eq)

        b_prime = -(volume_eq * a3 / a2 + 1)

        fit_dict["fit_type"] = "polynomial"
        fit_dict["fit_order"] = fit_order
        fit_dict["volume_eq"] = volume_eq
        fit_dict["energy_eq"] = e_eq
        fit_dict["bulkmodul_eq"] = eV_div_A3_to_GPa * volume_eq * a2
        fit_dict["b_prime_eq"] = b_prime
        fit_dict["least_square_error"] = self.get_error(volume_lst, energy_lst, p_fit)
        return fit_dict

    def _fit_leastsq(self, volume_lst, energy_lst, fittype='birchmurnaghan'):
        """
        Internal helper function for the least square fit

        Args:
            volume_lst (list/numpy.dnarray/None): vector of volumes
            energy_lst (list/numpy.dnarray/None): vector of energies
            fittype (str): on of the following ['birch', 'birchmurnaghan', 'murnaghan', 'pouriertarantola', 'vinet']

        Returns:
            list: [E0, B0, BP, V0], [E0_err, B0_err, BP_err, V0_err]
        """
        try:
            import matplotlib.pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        vol_lst = np.array(volume_lst).flatten()
        eng_lst = np.array(energy_lst).flatten()
        a, b, c = plt.polyfit(vol_lst, eng_lst, 2)
        v0 = -b / (2 * a)
        ev_angs_to_gpa = 1e21 / scipy.constants.physical_constants['joule-electron volt relationship'][0]
        pfit_leastsq, perr_leastsq = self._fit_leastsq_funct([a * v0 ** 2 + b * v0 + c, 2 * a * v0 * ev_angs_to_gpa, 4, v0],
                                                              vol_lst, eng_lst, fitfunction, fittype)
        return pfit_leastsq, perr_leastsq  # [e0, b0, bP, v0]

    @staticmethod
    def _fit_leastsq_funct(p0, datax, datay, function, fittype):
        """
        Internal least square fit function

        Args:
            p0 (list): [E0, B0, BP, V0] list of fit parameters
            datax (float/numpy.dnarray): volumes to fit
            datay (float/numpy.dnarray): energies corresponding to the volumes
            fittype (str): on of the following ['birch', 'birchmurnaghan', 'murnaghan', 'pouriertarantola', 'vinet']

        Returns:
            list: [E0, B0, BP, V0], [E0_err, B0_err, BP_err, V0_err]
        """
        # http://stackoverflow.com/questions/14581358/getting-standard-errors-on-fitted-parameters-using-the-optimize-leastsq-method-i

        errfunc = lambda p, x, y, fittype: function(p, x, fittype) - y

        pfit, pcov, infodict, errmsg, success = \
            spy.leastsq(errfunc, p0, args=(datax, datay, fittype), full_output=1, epsfcn=0.0001)

        if (len(datay) > len(p0)) and pcov is not None:
            s_sq = (errfunc(pfit, datax, datay, fittype) ** 2).sum() / (len(datay) - len(p0))
            pcov = pcov * s_sq
        else:
            pcov = np.inf

        error = []
        for i in range(len(pfit)):
            try:
                error.append(np.absolute(pcov[i][i]) ** 0.5)
            except:
                error.append(0.00)
        pfit_leastsq = pfit
        perr_leastsq = np.array(error)
        return pfit_leastsq, perr_leastsq

    @staticmethod
    def get_error(x_lst, y_lst, p_fit):
        """

        Args:
            x_lst:
            y_lst:
            p_fit:

        Returns:
            numpy.dnarray
        """
        y_fit_lst = np.array(p_fit(x_lst))
        error_lst = (y_lst - y_fit_lst) ** 2
        return np.mean(error_lst)

    def fit_energy(self, volume_lst):
        """
        Gives the energy value for the corresponding energy volume fit defined in the fit dictionary.

        Args:
            volume_lst: list of volumes

        Returns:
            list of energies

        """
        if not self._fit_dict:
            return ValueError("parameter 'fit_dict' has to be defined!")
        v = volume_lst
        e0 = self._fit_dict["energy_eq"]
        b0 = self._fit_dict["bulkmodul_eq"] / 160.21766208
        b_p = self._fit_dict["b_prime_eq"]
        v0 = self._fit_dict["volume_eq"]
        if self._fit_dict['fit_type'] == 'birch':
            return self.birch(v, e0, b0, b_p, v0)
        elif self._fit_dict['fit_type'] == 'birchmurnaghan':
            return self.birchmurnaghan_energy(v, e0, b0, b_p, v0)
        elif self._fit_dict['fit_type'] == 'murnaghan':
            return self.murnaghan(v, e0, b0, b_p, v0)
        elif self._fit_dict['fit_type'] == 'pouriertarantola':
            return self.pouriertarantola(v, e0, b0, b_p, v0)
        else:
            return self.vinet_energy(v, e0, b0, b_p, v0)

    @staticmethod
    def birchmurnaghan_energy(V, E0, B0, BP, V0):
        'BirchMurnaghan equation from PRB 70, 224107'
        return birchmurnaghan_energy(V, E0, B0, BP, V0)

    @staticmethod
    def vinet_energy(V, E0, B0, BP, V0):
        'Vinet equation from PRB 70, 224107'
        return vinet_energy(V, E0, B0, BP, V0)

    @staticmethod
    def murnaghan(V, E0, B0, BP, V0):
        'From PRB 28,5480 (1983'
        return murnaghan(V, E0, B0, BP, V0)

    @staticmethod
    def birch(V, E0, B0, BP, V0):
        """
        From Intermetallic compounds: Principles and Practice, Vol. I: Principles
        Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
        paper downloaded from Web

        case where n=0
        """
        return birch(V, E0, B0, BP, V0)

    @staticmethod
    def pouriertarantola(V, E0, B0, BP, V0):
        return pouriertarantola(V, E0, B0, BP, V0)


# ToDo: not all abstract methods implemented
class Murnaghan(AtomisticParallelMaster):
    def __init__(self, project, job_name='murnaghan'):
        """

        Args:
            project:
            job_name:
        """
        super(Murnaghan, self).__init__(project, job_name)
        self.__name__ = 'Murnaghan'
        self.__version__ = '0.3.0'

        # print ("h5_path: ", self.project_hdf5._h5_path)

        # define default input
        self.input['num_points'] = (11, 'number of sample points')
        self.input['fit_type'] = ("polynomial", "['polynomial', 'birch', 'birchmurnaghan', 'murnaghan', 'pouriertarantola', 'vinet']")
        self.input['fit_order'] = (3, 'order of the fit polynom')
        self.input['vol_range'] = (0.1, 'relative volume variation around volume defined by ref_ham')

        self.debye_model = DebyeModel(self)
        self.fit_module = EnergyVolumeFit()

        self.fit_dict = None
        self._debye_T = None
        self._job_generator = MurnaghanJobGenerator(self)

    @property
    def fit(self):
        return self.debye_model

    @property
    def equilibrium_volume(self):
        return self.fit_dict["volume_eq"]

    @property
    def equilibrium_energy(self):
        return self.fit_dict["energy_eq"]

    def fit_polynomial(self, fit_order=3, vol_erg_dic=None):
        return self.poly_fit(fit_order=fit_order, vol_erg_dic=vol_erg_dic)

    def fit_murnaghan(self, vol_erg_dic=None):
        return self._fit_eos_general(vol_erg_dic=vol_erg_dic, fittype='murnaghan')

    def fit_birch_murnaghan(self, vol_erg_dic=None):
        return self._fit_eos_general(vol_erg_dic=vol_erg_dic, fittype='birchmurnaghan')

    def fit_vinet(self, vol_erg_dic=None):
        return self._fit_eos_general(vol_erg_dic=vol_erg_dic, fittype='vinet')

    def _fit_eos_general(self, vol_erg_dic=None, fittype='birchmurnaghan'):
        self._set_fit_module(vol_erg_dic=vol_erg_dic)
        fit_dict = self.fit_module.fit_eos_general(fittype=fittype)
        self.input['fit_type'] = fit_dict["fit_type"]
        self.input['fit_order'] = 0
        with self.project_hdf5.open('input') as hdf5_input:
            self.input.to_hdf(hdf5_input)
        with self.project_hdf5.open("output") as hdf5:
            hdf5["equilibrium_energy"] = fit_dict["energy_eq"]
            hdf5["equilibrium_volume"] = fit_dict["volume_eq"]
            hdf5["equilibrium_bulk_modulus"] = fit_dict["bulkmodul_eq"]
            hdf5["equilibrium_b_prime"] = fit_dict["b_prime_eq"]

        self.fit_dict = fit_dict
        return fit_dict

    def _fit_leastsq(self, volume_lst, energy_lst, fittype='birchmurnaghan'):
        return self.fit_module._fit_leastsq(volume_lst=volume_lst, energy_lst=energy_lst, fittype=fittype)

    def _set_fit_module(self, vol_erg_dic=None):
        if vol_erg_dic is not None:
            if "volume" in vol_erg_dic.keys() and "energy" in vol_erg_dic.keys():
                self.fit_module = EnergyVolumeFit(volume_lst=vol_erg_dic["volume"], energy_lst=vol_erg_dic["energy"])
            else:
                raise KeyError
        else:
            df = self.output_to_pandas()
            self.fit_module = EnergyVolumeFit(volume_lst=df["volume"].values, energy_lst=df["energy"].values)

    def poly_fit(self, fit_order=3, vol_erg_dic=None):
        self._set_fit_module(vol_erg_dic=vol_erg_dic)
        fit_dict = self.fit_module.fit_polynomial(fit_order=fit_order)
        if fit_dict is None:
            self._logger.warning("Minimum could not be found!")
        else:
            self.input['fit_type'] = fit_dict["fit_type"]
            self.input['fit_order'] = fit_dict["fit_order"]
            with self.project_hdf5.open('input') as hdf5_input:
                self.input.to_hdf(hdf5_input)
            with self.project_hdf5.open("output") as hdf5:
                hdf5["equilibrium_energy"] = fit_dict["energy_eq"]
                hdf5["equilibrium_volume"] = fit_dict["volume_eq"]
                hdf5["equilibrium_bulk_modulus"] = fit_dict["bulkmodul_eq"]
                hdf5["equilibrium_b_prime"] = fit_dict["b_prime_eq"]

            with self._hdf5.open("output") as hdf5:
                self.get_structure(iteration_step=-1).to_hdf(hdf5)

            self.fit_dict = fit_dict
        return fit_dict

    def list_structures(self):
        if self.ref_job.structure is not None:
            return [parameter[1] for parameter in self._job_generator.parameter_list]
        else:
            return []

    def collect_output(self):
        if self.server.run_mode.interactive:
            ham = self.project_hdf5.inspect(self.child_ids[0])
            erg_lst = ham["output/generic/energy_tot"]
            vol_lst = ham["output/generic/volume"]
            arg_lst = np.argsort(vol_lst)

            self._output["volume"] = vol_lst[arg_lst]
            self._output["energy"] = erg_lst[arg_lst]
        else:
            erg_lst, vol_lst, err_lst, id_lst = [], [], [], []
            for job_id in self.child_ids:
                ham = self.project_hdf5.inspect(job_id)
                print('job_id: ', job_id, ham.status)
                energy = ham["output/generic/energy_tot"][-1]
                volume = ham["output/generic/volume"][-1]
                erg_lst.append(np.mean(energy))
                err_lst.append(np.var(energy))
                vol_lst.append(volume)
                id_lst.append(job_id)
            vol_lst = np.array(vol_lst)
            erg_lst = np.array(erg_lst)
            err_lst = np.array(err_lst)
            id_lst = np.array(id_lst)
            arg_lst = np.argsort(vol_lst)

            self._output["volume"] = vol_lst[arg_lst]
            self._output["energy"] = erg_lst[arg_lst]
            self._output["error"] = err_lst[arg_lst]
            self._output["id"] = id_lst[arg_lst]

        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val
        if self.input['fit_type'] == "polynomial":
            self.fit_polynomial(fit_order=self.input['fit_order'])
        else:
            self._fit_eos_general(fittype=self.input['fit_type'])

    def plot(self, num_steps=100, plt_show=True):
        try:
            import matplotlib.pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        if not self.fit_dict:
            if self.input['fit_type'] == "polynomial":
                self.fit_polynomial(fit_order=self.input['fit_order'])
            else:
                self._fit_eos_general(fittype=self.input['fit_type'])
        df = self.output_to_pandas()
        vol_lst, erg_lst = df["volume"].values, df["energy"].values
        x_i = np.linspace(np.min(vol_lst), np.max(vol_lst), num_steps)
        color = 'blue'

        if self.fit_dict is not None:
            if self.input['fit_type'] == "polynomial":
                p_fit = np.poly1d(self.fit_dict["poly_fit"])
                least_square_error = self.fit_module.get_error(vol_lst, erg_lst, p_fit)
                plt.title("Murnaghan: error: " + str(least_square_error))
                plt.plot(x_i, p_fit(x_i), '-', label=self.input['fit_type'], color=color, linewidth=3)
            else:
                V0 = self.fit_dict["volume_eq"]
                E0 = self.fit_dict["energy_eq"]
                B0 = self.fit_dict["bulkmodul_eq"]
                BP = self.fit_dict["b_prime_eq"]
                if self.input['fit_type'].lower() == 'birchmurnaghan':
                    eng_fit_lst = birchmurnaghan_energy(x_i, E0, B0 / eV_div_A3_to_GPa, BP, V0)
                elif self.input['fit_type'].lower() == 'vinet':
                    eng_fit_lst = vinet_energy(x_i, E0, B0 / eV_div_A3_to_GPa, BP, V0)
                elif self.input['fit_type'].lower() == 'murnaghan':
                    eng_fit_lst = murnaghan(x_i, E0, B0 / eV_div_A3_to_GPa, BP, V0)
                elif self.input['fit_type'].lower() == 'pouriertarantola':
                    eng_fit_lst = pouriertarantola(x_i, E0, B0 / eV_div_A3_to_GPa, BP, V0)
                elif self.input['fit_type'].lower() == 'birch':
                    eng_fit_lst = birch(x_i, E0, B0, BP, V0)
                else:
                    raise ValueError
                plt.plot(x_i, eng_fit_lst, '-', label=self.input['fit_type'], color=color, linewidth=3)

        plt.plot(vol_lst, erg_lst, 'x', color=color, markersize=20)
        plt.legend()
        plt.xlabel("Volume ($\AA^3$)")
        plt.ylabel("energy (eV)")
        if plt_show:
            plt.show()

    def get_structure(self, iteration_step=-1):
        """

        Returns: Structure with equilibrium volume

        """
        if not (self.structure is not None):
            raise AssertionError()
        if iteration_step == -1:
            snapshot = self.structure.copy()
            old_vol = snapshot.get_volume()
            new_vol = self["output/equilibrium_volume"]
            k = (new_vol / old_vol) ** (1. / 3.)
            new_cell = snapshot.cell * k
            snapshot.set_cell(new_cell, scale_atoms=True)
            return snapshot
        elif iteration_step == 0:
            return self.structure
        else:
            raise ValueError('iteration_step should be either 0 or -1.')
