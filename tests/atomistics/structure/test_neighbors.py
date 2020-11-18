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

    def test_allow_ragged(self):
        struct = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc').repeat(10)
        del struct[0]
        neigh = struct.get_neighbors_by_distance(cutoff_radius=3)
        self.assertTrue(neigh.allow_ragged)
        self.assertTrue(isinstance(neigh.indices, list))
        indices = neigh.indices.copy()
        with self.assertRaises(ValueError):
            neigh.allow_ragged = 'yes'
        neigh.allow_ragged = False
        self.assertTrue(isinstance(neigh.indices, np.ndarray))
        self.assertGreater(len(neigh.indices[0]), len(indices[0]))
        with self.assertRaises(IndexError):
            struct.positions[neigh.indices] = struct.positions[neigh.indices]
        neigh.allow_ragged = True
        self.assertTrue(np.array_equal(neigh.indices[0], indices[0]))
        neigh = struct.get_neighbors(cutoff_radius=3, num_neighbors=None)
        self.assertFalse(neigh.allow_ragged)
        self.assertTrue(isinstance(neigh.indices, np.ndarray))

    def test_getters(self):
        struct = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc')
        neigh = struct.get_neighbors()
        self.assertTrue(np.array_equal(neigh.distances, neigh.get_distances()))
        self.assertTrue(np.array_equal(neigh.vecs, neigh.get_vectors()))
        self.assertTrue(np.array_equal(neigh.indices, neigh.get_indices()))

    def test_getter_and_ragged(self):
        struct = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc').repeat(2)
        del struct[0]
        neigh = struct.get_neighbors_by_distance(cutoff_radius=3)
        vecs = neigh.get_vectors(allow_ragged=False)
        distances = neigh.get_distances(allow_ragged=False)
        indices = neigh.get_indices(allow_ragged=False)
        neigh.allow_ragged = False
        self.assertTrue(np.array_equal(neigh.distances, distances))
        self.assertTrue(np.array_equal(neigh.indices, indices))
        self.assertTrue(np.array_equal(neigh.vecs, vecs))

    def test_get_neighborhood_single(self):
        struct = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc')
        neigh = struct.get_neighborhood(np.random.random(3), cutoff_radius=3)
        distances = neigh.distances.copy()
        neigh.allow_ragged = False
        self.assertTrue(np.array_equal(neigh.distances, distances))
        neigh.allow_ragged = True
        self.assertTrue(np.array_equal(neigh.distances, distances))
        neigh.allow_ragged = False
        self.assertTrue(np.array_equal(neigh.distances, distances))

    def test_get_neighborhood(self):
        struct = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc')
        struct.positions += 0.01*(2*np.random.random(struct.positions.shape)-1)
        struct = struct.center_coordinates_in_unit_cell()
        positions = np.random.random((2, 3)).dot(struct.cell)
        neigh = struct.get_neighborhood(positions)
        positions = np.random.random((2, 3)).dot(struct.cell)
        new_neigh = neigh.get_neighborhood(positions)
        self.assertTrue(np.array_equal(new_neigh.distances, neigh.get_distances(positions)))
        self.assertTrue(np.array_equal(new_neigh.vecs, neigh.get_vectors(positions)))
        self.assertTrue(np.array_equal(new_neigh.indices, neigh.get_indices(positions)))

    def test_get_neighbors(self):
        struct = CrystalStructure(elements='Fe', lattice_constants=2.85, bravais_basis='bcc').repeat(10)
        cell = struct.cell.copy()
        cell += np.random.random((3,3))-0.5
        struct.positions += np.random.random((len(struct), 3))-0.5
        struct.set_cell(cell, scale_atoms=True)
        neigh = struct.get_neighbors()
        self.assertAlmostEqual(np.absolute(neigh.distances-np.linalg.norm(neigh.vecs, axis=-1)).max(), 0)
        myself = np.ones_like(neigh.indices)
        myself = myself*np.arange(len(myself))[:,np.newaxis]
        dist = struct.get_distances(myself.flatten(), neigh.indices.flatten(), mic=True)
        self.assertAlmostEqual(np.absolute(dist-neigh.distances.flatten()).max(), 0)
        vecs = struct.get_distances(myself.flatten(), neigh.indices.flatten(), mic=True, vector=True)
        self.assertAlmostEqual(np.absolute(vecs-neigh.vecs.reshape(-1, 3)).max(), 0)
        dist = struct.get_scaled_positions()
        dist = dist[:,np.newaxis,:]-dist[np.newaxis,:,:]
        dist -= np.rint(dist)
        dist = np.einsum('nmi,ij->nmj', dist, struct.cell)
        dist = np.linalg.norm(dist, axis=-1).flatten()
        dist = dist[dist>0]
        self.assertAlmostEqual(neigh.distances.min(), dist.min())

    def test_pbc_false(self):
        struct = CrystalStructure(elements='Fe', lattice_constants=2.85, bravais_basis='bcc').repeat(10)
        struct.pbc = False
        cell = struct.cell.copy()
        cell += np.random.random((3,3))-0.5
        struct.set_cell(cell, scale_atoms=True)
        neigh = struct.get_neighbors()
        self.assertAlmostEqual(np.absolute(neigh.distances-np.linalg.norm(neigh.vecs, axis=-1)).max(), 0)
        myself = np.ones_like(neigh.indices)
        myself = myself*np.arange(len(myself))[:,np.newaxis]
        dist = np.linalg.norm(struct.positions[myself]-struct.positions[neigh.indices], axis=-1)
        self.assertAlmostEqual(np.absolute(dist-neigh.distances).max(), 0)

    def test_fe_large(self):
        struct = CrystalStructure(elements='Fe', lattice_constants=2.85, bravais_basis='bcc').repeat(10)
        neigh = struct.get_neighbors()
        self.assertAlmostEqual(np.absolute(neigh.distances-np.linalg.norm(neigh.vecs, axis=-1)).max(), 0)
        self.assertAlmostEqual(neigh.vecs[neigh.shells==1].sum(), 0)
        self.assertAlmostEqual(neigh.vecs[0, neigh.shells[0]==1].sum(), 0)

    def test_fe_small(self):
        struct = CrystalStructure(elements='Fe', lattice_constants=2.85, bravais_basis='bcc')
        neigh = struct.get_neighbors()
        self.assertAlmostEqual(neigh.vecs[neigh.shells==1].sum(), 0)

    def test_al_large(self):
        struct = CrystalStructure(elements='Al', lattice_constants=4.04, bravais_basis='fcc').repeat(10)
        neigh = struct.get_neighbors()
        self.assertAlmostEqual(np.absolute(neigh.distances-np.linalg.norm(neigh.vecs, axis=-1)).max(), 0)
        self.assertAlmostEqual(neigh.vecs[neigh.shells==1].sum(), 0)
        self.assertAlmostEqual(neigh.vecs[0, neigh.shells[0]==1].sum(), 0)

    def test_al_small(self):
        struct = CrystalStructure(elements='Al', lattice_constants=4.04, bravais_basis='fcc')
        neigh = struct.get_neighbors()
        self.assertAlmostEqual(np.absolute(neigh.distances-np.linalg.norm(neigh.vecs, axis=-1)).max(), 0)
        self.assertAlmostEqual(neigh.vecs[neigh.shells==1].sum(), 0)
        with self.assertRaises(ValueError):
            struct.get_neighbors(num_neighbors=0)

    def test_wrapped_positions(self):
        structure = CrystalStructure(elements='Al', lattice_constants=4, bravais_basis='fcc').repeat(2)
        neigh = structure.get_neighbors()
        distances = neigh.distances
        new_positions = structure.positions+structure.cell.diagonal()
        self.assertFalse(np.all(np.isclose(distances, neigh.get_distances(new_positions, num_neighbors=13)[:,1:])))
        neigh.wrap_positions = True
        self.assertTrue(np.all(np.isclose(distances, neigh.get_distances(new_positions, num_neighbors=13)[:,1:])))

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
        structure.center_coordinates_in_unit_cell()
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

    def test_find_neighbors_by_vector(self):
        basis = Atoms(symbols=2*["Fe"],
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)],
                      cell=np.identity(3),
                      pbc=True)
        neigh = basis.get_neighbors(num_neighbors=14)
        id_lst, dist = neigh.find_neighbors_by_vector([0, 0, 1],
                                                      deviation=True)
        self.assertEqual(len(np.unique(np.unique(id_lst, return_counts=True)[1])), 1)
        self.assertLess(np.linalg.norm(dist), 1.0e-4)
        id_lst = neigh.find_neighbors_by_vector([0, 0, 0])
        self.assertTrue(np.array_equal(id_lst, np.arange(len(basis))))

if __name__ == "__main__":
    unittest.main()
