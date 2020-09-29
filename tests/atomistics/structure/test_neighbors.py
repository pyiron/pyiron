# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
from pyiron.atomistics.structure.atoms import Atoms, CrystalStructure


class TestAtoms(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        pass

    @classmethod
    def setUpClass(cls):
        pass

    def test_get_global_shells(self):
        structure = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc').repeat(2)
        neigh = structure.get_neighbors()
        self.assertTrue(np.array_equal(neigh.shells, neigh.get_global_shells()))
        structure += Atoms(elements='C', positions=[[0, 0, 0.5*4]])
        neigh = structure.get_neighbors()
        self.assertFalse(np.array_equal(neigh.shells, neigh.get_global_shells()))
        structure = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc').repeat(2)
        neigh = structure.get_neighbors()
        shells = neigh.get_global_shells()
        structure.positions += 0.01*(np.random.random((len(structure), 3))-0.5)
        neigh = structure.get_neighbors()
        self.assertTrue(np.array_equal(shells, neigh.get_global_shells(cluster_by_vecs=True, cluster_by_distances=True)))
        neigh.reset_clusters()
        self.assertTrue(np.array_equal(shells, neigh.get_global_shells(cluster_by_vecs=True)))
        self.assertFalse(np.array_equal(shells, neigh.get_global_shells()))

    def test_get_local_shells(self):
        structure = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc').repeat(2)
        neigh = structure.get_neighbors()
        shells = neigh.get_local_shells()
        structure.positions += 0.01*(np.random.random((len(structure), 3))-0.5)
        neigh = structure.get_neighbors()
        self.assertTrue(np.array_equal(shells, neigh.get_local_shells(cluster_by_vecs=True, cluster_by_distances=True)))
        neigh.reset_clusters()
        self.assertTrue(np.array_equal(shells, neigh.get_local_shells(cluster_by_vecs=True)))
        self.assertFalse(np.array_equal(shells, neigh.get_local_shells()))

    def test_get_shell_matrix(self):
        structure = CrystalStructure(elements='Fe', lattice_constants=2.83, bravais_basis='bcc').repeat(2)
        structure[0] = 'Ni'
        neigh = structure.get_neighbors(num_neighbors=8)
        mat = neigh.get_shell_matrix()
        self.assertEqual(mat[0].sum(), 8*len(structure))
        mat = neigh.get_shell_matrix(chemical_symbols=['Fe', 'Ni'])
        self.assertEqual(mat[0].sum(), 16)

if __name__ == "__main__":
    unittest.main()
