# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import division, print_function
from ase.atoms import Atoms as ASEAtoms, get_distances as ase_get_distances, Atom as ASEAtom
import ast
from copy import copy
import numpy as np
from pyiron_base import Settings

__author__ = "Joerg Neugebauer, Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Osamu Waseda"
__email__ = "waseda@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


class Neighbors(object):
    """
    Class for storage of the neighbor information for a given atom based on the KDtree algorithm
    """

    def __init__(self, ref_structure, tolerance=2):
        self.distances = None
        self.vecs = None
        self.indices = None
        self._shells = None
        self._tolerance = tolerance
        self._ref_structure = ref_structure

    @property
    def shells(self):
        if self._shells is None:
            if self.distances is None:
                return None
            self._shells = []
            for dist in self.distances:
                self._shells.append(np.unique(np.round(dist, decimals=self._tolerance), return_inverse=True)[1]+1)
            if isinstance(self.distances, np.ndarray):
                self._shells = np.array(self._shells)
            return self._shells

    def get_shell_dict(self, max_shell=2):
        """

        Args:
            max_shell (int): maximum number of shells

        Returns:

        """

        shells = self.shells[0]
        dist = self.distances[0]
        shell_dict = {}
        for i_shell in set(shells):
            if i_shell > max_shell:
                break
            shell_dict[i_shell] = np.mean(dist[shells == i_shell])
            # print ("shells: ", i_shell, shell_dict[i_shell])
        if not (max(shell_dict.keys()) == max_shell):
            raise AssertionError()
        return shell_dict

    def get_global_shells(self, decimals=5):
        """
        Set shell indices based on all distances available in the system instead of
        setting them according to the local distances (in contrast to shells defined
        as an attribute in this class)

        Args:
            decimals (int): decimals in np.round for rounding up distances

        Returns:
            shells (numpy.ndarray): shell indices (cf. shells)
        """
        if self.distances is None:
            raise ValueError('neighbors not set')
        distances = np.unique(np.round(a=self.distances, decimals=decimals))
        shells = self.distances[:,:,np.newaxis]-distances[np.newaxis,np.newaxis,:]
        shells = np.absolute(shells).argmin(axis=-1)+1
        return shells

    def get_shell_matrix(
        self, shell=None, restraint_matrix=None
    ):
        """

        Args:
            shell (int/None): shell number. If None, all shells are returned
            restraint_matrix: NxN matrix with True or False, where False will remove the entries.
                              If an integer is given the sum of the chemical indices corresponding to the number will
                              be set to True and the rest to False

        Returns:
            NxN matrix with 1 for the pairs of atoms in the given shell

        """
        if shell is not None and shell<=0:
            raise ValueError("Parameter 'shell' must be an integer greater than 0")
        neigh_list = self.get_neighbors(
            num_neighbors=num_neighbors, id_list=id_list, tolerance=tolerance
        )
        Natom = len(self)
        if shell is None:
            shell_lst = np.unique(neigh_list.shells)
        else:
            shell_lst = np.array([shell]).flatten()
        if restraint_matrix is None:
            restraint_matrix = np.ones((Natom, Natom)) == 1
        elif type(restraint_matrix) == list and len(restraint_matrix) == 2:
            restraint_matrix = np.outer(
                1 * (self.get_chemical_symbols() == restraint_matrix[0]),
                1 * (self.get_chemical_symbols() == restraint_matrix[1]),
            )
            restraint_matrix = (restraint_matrix + restraint_matrix.transpose()) > 0
        shell_matrix_lst = []
        for shell in shell_lst:
            shell_matrix = np.zeros((Natom, Natom))
            for ii, ss in enumerate(neigh_list.shells):
                unique, counts = np.unique(
                    neigh_list.indices[ii][ss == np.array(shell)], return_counts=True
                )
                shell_matrix[ii][unique] = counts
            shell_matrix[restraint_matrix == False] = 0
            shell_matrix_lst.append(shell_matrix)
        if len(shell_matrix_lst)==1:
            return shell_matrix_lst[0]
        else:
            return shell_matrix_lst



