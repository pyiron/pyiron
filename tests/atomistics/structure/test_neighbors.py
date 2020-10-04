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
        mat = neigh.get_shell_matrix(chemical_pair=['Fe', 'Ni'])
        self.assertEqual(mat[0].sum(), 16)
        mat = neigh.get_shell_matrix(chemical_pair=['Ni', 'Ni'])
        self.assertEqual(mat[0].sum(), 0)

    def test_cluster_analysis(self):
        basis = CrystalStructure("Al", bravais_basis="fcc", lattice_constants=4.2).repeat(10)
        neigh = basis.get_neighbors(num_neighbors=100)
        key, counts = neigh.cluster_analysis(id_list=[0,1], return_cluster_sizes=True)
        self.assertTrue(np.array_equal(key[1], [0,1]))
        self.assertEqual(counts[0], 2)
        key, counts = neigh.cluster_analysis(id_list=[0,int(len(basis)/2)], return_cluster_sizes=True)
        self.assertTrue(np.array_equal(key[1], [0]))
        self.assertEqual(counts[0], 1)

    def test_get_bonds(self):
        basis = CrystalStructure("Al", bravais_basis="fcc", lattice_constants=4.2).repeat(5)
        neigh = basis.get_neighbors(num_neighbors=20)
        bonds = neigh.get_bonds()
        self.assertTrue(np.array_equal(np.sort(bonds[0]['Al'][0]),
                        np.sort(neigh.indices[0, neigh.shells[0]==1])))

if __name__ == "__main__":
    unittest.main()
