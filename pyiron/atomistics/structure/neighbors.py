# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

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
        """
            Returns the cell numbers of each atom according to the distances
        """
        if self._shells is None:
            if self.distances is None:
                return None
            self._shells = []
            for dist in self.distances:
                self._shells.append(np.unique(np.round(dist, decimals=self._tolerance), return_inverse=True)[1]+1)
            if isinstance(self.distances, np.ndarray):
                self._shells = np.array(self._shells)
            return self._shells

    def update_vectors(self):
        """
            Update vecs and distances with the same indices
        """
        if np.max(np.absolute(self.vecs)) > 0.49*np.min(np.linalg.norm(self._ref_structure.cell, axis=-1)):
            raise AssertionError('Largest distance value is larger than half the box -> rerun get_neighbors')
        myself = np.ones_like(self.indices)*np.arange(len(self.indices))[:, np.newaxis]
        vecs = self._ref_structure.get_distances(
            myself.flatten(), self.indices.flatten(), mic=np.all(self._ref_structure.pbc), vector=True
        )
        self.vecs = vecs.reshape(self.vecs.shape)
        self.distances = np.linalg.norm(self.vecs, axis=-1)

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
        self, shell_numbers=None, restraint_matrix=None
    ):
        """

        Args:
            shell_numbers (int/None): shell number. If None, all shells are returned
            restraint_matrix: NxN matrix with True or False, where False will remove the entries.
                              If an integer is given the sum of the chemical indices corresponding to the number will
                              be set to True and the rest to False

        Returns:
            NxN matrix with 1 for the pairs of atoms in the given shell

        """
        if shell_numbers is not None and shell_numbers<=0:
            raise ValueError("Parameter 'shell' must be an integer greater than 0")
        Natom = len(self._ref_structure)
        if shell_numbers is None:
            shell_lst = np.unique(self.shells)
        else:
            shell_lst = np.array([shell_numbers]).flatten()
        if restraint_matrix is None:
            restraint_matrix = np.ones((Natom, Natom)) == 1
        elif isinstance(restraint_matrix, list) and len(restraint_matrix) == 2:
            restraint_matrix = np.outer(
                1 * (self._ref_structure.get_chemical_symbols() == restraint_matrix[0]),
                1 * (self._ref_structure.get_chemical_symbols() == restraint_matrix[1]),
            )
            restraint_matrix = (restraint_matrix + restraint_matrix.transpose()) > 0
        shell_matrix_lst = []
        for shell in shell_lst:
            shell_matrix = np.zeros((Natom, Natom))
            for ii, ss in enumerate(self.shells):
                unique, counts = np.unique(
                    self.indices[ii][ss == np.array(shell)], return_counts=True
                )
                shell_matrix[ii][unique] = counts
            shell_matrix[restraint_matrix == False] = 0
            shell_matrix_lst.append(shell_matrix)
        if len(shell_matrix_lst)==1:
            return shell_matrix_lst[0]
        else:
            return shell_matrix_lst

    def find_neighbors_by_vector(self, vector, deviation=False):
        """
        Args:
            vector (list/np.ndarray): vector by which positions are translated (and neighbors are searched)
            deviation (bool): whether to return distance between the expect positions and real positions

        Returns:
            np.ndarray: list of id's for the specified translation

        Example:
            a_0 = 2.832
            structure = pr.create_structure('Fe', 'bcc', a_0)
            id_list = structure.find_neighbors_by_vector([0, 0, a_0])
            # In this example, you get a list of neighbor atom id's at z+=a_0 for each atom.
            # This is particularly powerful for SSA when the magnetic structure has to be translated
            # in each direction.
        """

        dist = np.linalg.norm(self.vecs-np.array(vector), axis=-1)
        if deviation:
            return self.indices[np.arange(len(dist)), np.argmin(dist, axis=-1)], np.min(dist, axis=-1)
        return self.indices[np.arange(len(dist)), np.argmin(dist, axis=-1)]



