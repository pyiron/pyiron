# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import numpy as np

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class Dos(object):

    """
    The DOS class stores all information to store and retrieve the total and resolved density of states from an
    electronic structure calculation.

    Args:
        n_bins (int): Number of histogram bins required to calculate the DOS
        es_obj: The pyiron.objects.waves.core.ElectronicStructure instance for which the DOS has to be computed
        eigenvalues (list/numpy.ndarray): If es-obj is None, the eigenvalues could be specified as a list

    """

    def __init__(self, n_bins=100, es_obj=None, eigenvalues=None, bin_density=None):
        self.orbital_dict = {"s": [0], "p": [1, 2, 3], "d": [4, 5, 6, 7, 8]}
        self.n_bins = n_bins
        self.es_obj = es_obj
        dos_min = np.min(self.es_obj.eigenvalues)
        dos_max = np.max(self.es_obj.eigenvalues)
        if bin_density is not None:
            n_bins = int((dos_max - dos_min) * bin_density)
        if es_obj is not None:
            self.t_dos, self.energies = np.histogram(
                self.es_obj.eigenvalues, bins=int(n_bins), density=True
            )
        else:
            self.t_dos, self.energies = np.histogram(
                eigenvalues, bins=int(n_bins), density=True
            )
        self.energies = self.energies[1:] - (
            (self.energies[1] - self.energies[0]) / 2.0
        )

    def plot_total_dos(self, **kwargs):
        """
        Plots the total DOS

        Args:
            **kwargs: Variables for matplotlib.pylab.plot customization (linewidth, linestyle, etc.)

        Returns:
            matplotlib.pylab.plot
        """
        try:
            import matplotlib.pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        fig = plt.figure(1, figsize=(6, 4))
        ax1 = fig.add_subplot(111)
        ax1.set_xlabel("E (eV)", fontsize=14)
        ax1.set_ylabel("DOS", fontsize=14)
        plt.fill_between(self.energies, self.t_dos, **kwargs)
        return plt

    def plot_orbital_resolved_dos(self, **kwargs):
        """
        Plots the orbital resolved DOS

        Args:
            **kwargs: Variable for matplotlib.pylab.plot customization (linewidth, linestyle, etc.)

        Returns:
            matplotlib.pylab.plot
        """
        try:
            import matplotlib.pylab as plt
        except ImportError:
            import matplotlib.pyplot as plt
        if not (self.es_obj.grand_dos_matrix is not None):
            raise NoResolvedDosError(
                "Can not plot the orbital resolved dos since resolved dos values are not"
                " available"
            )
        plot = self.plot_total_dos()
        for key, val in self.orbital_dict.items():
            r_dos = self.get_orbital_resolved_dos(val)
            plt.plot(self.energies, r_dos, label=key, **kwargs)
        plot.legend()
        return plot

    def get_spin_resolved_dos(self, spin_indices):
        """
        Gives the dos contribution of a given indices of spin as arranged in the
        pyiron.objects.waves.ElectronicStructure instance.

        Args:
            spin_indices (list/numpy.ndarray): The index/indices of the spins for which the dos contribution is required

        Returns:
            numpy.ndarray: The required dos

        """
        if not (self.es_obj.grand_dos_matrix is not None):
            raise NoResolvedDosError(
                "Can not get the spin resolved dos since resolved dos values are not"
                " available"
            )

        grand_sum = np.sum(self.es_obj.grand_dos_matrix)
        tot_val = self.es_obj.grand_dos_matrix.copy() / grand_sum
        n_spin, n_kpts, n_bands, _, _ = np.shape(tot_val)
        k = 0
        b = 0
        r_dos = np.zeros_like(self.t_dos)
        w_dos = np.zeros_like(self.t_dos)
        for i, e in enumerate(self.es_obj.eigenvalues):
            weight = np.sum(tot_val[spin_indices, k, b, :, :])
            weight_sum = np.sum(tot_val[:, k, b, :, :])
            if b < n_bands - 1:
                b += 1
            else:
                b = 0
                k += 1
            index = len(self.energies[self.energies < e]) - 1
            if index >= 0:
                r_dos[index] = r_dos[index] + weight
                w_dos[index] = w_dos[index] + weight_sum
            else:
                r_dos[0] = r_dos[0] + weight
                w_dos[0] = w_dos[0] + weight_sum
        ind_0 = np.argwhere(w_dos < 1e-8).flatten()
        ind_1 = np.argwhere(w_dos >= 1e-8).flatten()
        r_dos[ind_1] /= w_dos[ind_1]
        r_dos[ind_0] = 0.0
        return r_dos * self.t_dos

    def get_spatially_resolved_dos(self, atom_indices, spin_indices=0):
        """
        Gives the dos contribution of a given indices of atoms as arranged in the
        pyiron.objects.waves.ElectronicStructure instance.

        Args:
            atom_indices (list/numpy.ndarray): The index/indices of the atoms for which the dos contribution is required
            spin_indices (list/numpy.ndarray): The index/indices of the spins for which the dos contribution is required

        Returns:
            numpy.ndarray: The required dos

        """
        if not (self.es_obj.grand_dos_matrix is not None):
            raise NoResolvedDosError(
                "Can not get the spatially resolved dos since resolved dos values are not"
                " available"
            )
        grand_sum = np.sum(self.es_obj.grand_dos_matrix)
        tot_val = self.es_obj.grand_dos_matrix.copy() / grand_sum
        _, n_kpts, n_bands, _, _ = np.shape(tot_val)
        k = 0
        b = 0
        r_dos = np.zeros_like(self.t_dos)
        w_dos = np.zeros_like(self.t_dos)
        for i, e in enumerate(self.es_obj.eigenvalues):
            weight = np.sum(tot_val[spin_indices, k, b, atom_indices, :])
            weight_sum = np.sum(tot_val[spin_indices, k, b, :, :])
            if b < n_bands - 1:
                b += 1
            else:
                b = 0
                k += 1
            index = len(self.energies[self.energies < e]) - 1
            if index >= 0:
                r_dos[index] = r_dos[index] + weight
                w_dos[index] = w_dos[index] + weight_sum
            else:
                r_dos[0] = r_dos[0] + weight
                w_dos[0] = w_dos[0] + weight_sum
        ind_0 = np.argwhere(w_dos < 1e-8).flatten()
        ind_1 = np.argwhere(w_dos >= 1e-8).flatten()
        r_dos[ind_1] /= w_dos[ind_1]
        r_dos[ind_0] = 0.0
        return r_dos * self.t_dos

    def get_orbital_resolved_dos(self, orbital_indices, spin_indices=0):
        """
        Gives the dos contribution of a given indices of orbitals as arranged in the
        pyiron.objects.waves.ElectronicStructure instance.

        Args:
            orbital_indices (list/numpy.ndarray): The index/indices of the orbitals for which the dos contribution is required
            spin_indices (list/numpy.ndarray): The index/indices of the spins for which the dos contribution is required

        Returns:
            numpy.ndaray: The required dos

        """
        if not (self.es_obj.grand_dos_matrix is not None):
            raise NoResolvedDosError(
                "Can not get the orbital resolved dos since resolved dos values are not"
                " available"
            )
        grand_sum = np.sum(self.es_obj.grand_dos_matrix)
        tot_val = self.es_obj.grand_dos_matrix.copy() / grand_sum
        _, n_kpts, n_bands, _, _ = np.shape(tot_val)
        k = 0
        b = 0
        r_dos = np.zeros_like(self.t_dos)
        w_dos = np.zeros_like(self.t_dos)
        for i, e in enumerate(self.es_obj.eigenvalues):
            weight = np.sum(tot_val[spin_indices, k, b, :, orbital_indices])
            weight_sum = np.sum(tot_val[spin_indices, k, b, :, :])
            if b < n_bands - 1:
                b += 1
            else:
                b = 0
                k += 1
            index = len(self.energies[self.energies < e]) - 1
            if index >= 0:
                r_dos[index] = r_dos[index] + weight
                w_dos[index] = w_dos[index] + weight_sum
            else:
                r_dos[0] = r_dos[0] + weight
                w_dos[0] = w_dos[0] + weight_sum
        ind_0 = np.argwhere(w_dos < 1e-8).flatten()
        ind_1 = np.argwhere(w_dos >= 1e-8).flatten()
        r_dos[ind_1] /= w_dos[ind_1]
        r_dos[ind_0] = 0.0
        return r_dos * self.t_dos

    def get_spatial_orbital_resolved_dos(
        self, atom_indices, orbital_indices, spin_indices=0
    ):
        """
        Gives the dos contribution of a given indices of atoms as well as orbitals as arranged in the
        pyiron.objects.waves.ElectronicStructure instance.

        Args:
            atom_indices (list/numpy.ndarray): The index/indices of the atoms for which the dos contribution is required
            orbital_indices (list/numpy.ndarray): The index/indices of the orbitals for which the dos contribution is required
            spin_indices (list/numpy.ndarray): The index/indices of the spins for which the dos contribution is required

        Returns:
            numpy.ndaray: The required dos
        """
        if not (self.es_obj.grand_dos_matrix is not None):
            raise NoResolvedDosError(
                "Can not get the resolved dos since resolved dos values are not"
                " available"
            )
        grand_sum = np.sum(self.es_obj.grand_dos_matrix)
        tot_val = self.es_obj.grand_dos_matrix.copy() / grand_sum
        _, n_kpts, n_bands, _, _ = np.shape(tot_val)
        k = 0
        b = 0
        r_dos = np.zeros_like(self.t_dos)
        w_dos = np.zeros_like(self.t_dos)
        for i, e in enumerate(self.es_obj.eigenvalues):
            weight = np.sum(
                [
                    np.sum(tot_val[spin_indices, k, b, atom_indices, o])
                    for o in orbital_indices
                ]
            )
            weight_sum = np.sum(tot_val[spin_indices, k, b, :, :])
            if b < n_bands - 1:
                b += 1
            else:
                b = 0
                k += 1
            index = len(self.energies[self.energies < e]) - 1
            if index >= 0:
                r_dos[index] = r_dos[index] + weight
                w_dos[index] = w_dos[index] + weight_sum
            else:
                r_dos[0] = r_dos[0] + weight
                w_dos[0] = w_dos[0] + weight_sum
        ind_0 = np.argwhere(w_dos < 1e-8).flatten()
        ind_1 = np.argwhere(w_dos >= 1e-8).flatten()
        r_dos[ind_1] /= w_dos[ind_1]
        r_dos[ind_0] = 0.0
        return r_dos * self.t_dos


class NoResolvedDosError(Exception):
    """
    Raised when information on the resolved dos in unavailable
    """

    pass
