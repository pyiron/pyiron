import unittest
import sys
import numpy as np
import os
from pyiron_atomistics.structure.atom import Atom
from pyiron_atomistics.structure.atoms import Atoms, CrystalStructure
from pyiron_atomistics.structure.sparse_list import SparseList
from pyiron_atomistics.structure.periodic_table import PeriodicTable, ChemicalElement
from pyiron_base.objects.generic.hdfio import FileHDFio


class TestAtoms(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        if sys.version_info[0] >= 3:
            file_location = os.path.dirname(os.path.abspath(__file__))
            if os.path.isfile(os.path.join(file_location, "/static/pyiron_atomistics/test_hdf")):
                os.remove(os.path.join(file_location, "/static/pyiron_atomistics/test_hdf"))

    def setUp(self):
        pass
        self.CO2 = Atoms("CO2", positions=[[0, 0, 0], [0, 0, 1.5], [0, 1.5, 0]])
        C = Atom('C').element
        self.C3 = Atoms([C, C, C], positions=[[0, 0, 0], [0, 0, 2], [0, 2, 0]])
        self.C2 = Atoms(2 * [Atom('C')])

    def test__init__(self):
        pos, cell = generate_fcc_lattice()
        pse = PeriodicTable()
        el = pse.element("Al")
        basis = Atoms()
        self.assertIsInstance(basis, Atoms)
        self.assertIsInstance(basis.info, dict)
        self.assertIsInstance(basis.arrays, dict)
        self.assertIsInstance(basis.adsorbate_info, dict)
        self.assertIsInstance(basis.units, dict)
        self.assertIsInstance(basis.pbc, (bool, list, np.ndarray))
        self.assertIsInstance(basis.indices, np.ndarray)
        self.assertIsNone(basis._internal_positions)
        self.assertIsNone(basis.positions)
        self.assertIsNone(basis.scaled_positions)
        self.assertIsInstance(basis.species, list)
        self.assertIsInstance(basis.elements, np.ndarray)
        self.assertIsNone(basis.cell)
        basis = Atoms(symbols='Al', positions=pos, cell=cell)
        self.assertIsInstance(basis, Atoms)
        self.assertEqual(basis.get_spacegroup()["Number"], 225)
        basis = Atoms(elements='Al', positions=pos, cell=cell)
        self.assertIsInstance(basis, Atoms)
        basis = Atoms(elements=['Al'], positions=pos, cell=cell)
        self.assertIsInstance(basis, Atoms)
        self.assertRaises(ValueError, Atoms, symbols="Pt", elements='Al', positions=pos, cell=cell)
        basis = Atoms(numbers=[13], positions=pos, cell=cell)
        self.assertEqual(basis.get_majority_species()[1], "Al")
        basis = Atoms(species=[el], indices=[0], positions=pos, cell=cell)
        self.assertEqual(basis.get_majority_species()[1], "Al")
        self.assertIsInstance(basis, Atoms)
        self.assertIsInstance(basis.info, dict)
        self.assertIsInstance(basis.arrays, dict)
        self.assertIsInstance(basis.adsorbate_info, dict)
        self.assertIsInstance(basis.units, dict)
        self.assertIsInstance(basis.pbc, (bool, list, np.ndarray))
        self.assertIsInstance(basis.indices, np.ndarray)
        self.assertIsInstance(basis.species, list)
        self.assertIsInstance(basis.cell, np.ndarray)
        self.assertIsInstance(basis._internal_positions, np.ndarray)
        self.assertIsInstance(basis.positions, np.ndarray)
        self.assertIsInstance(basis.scaled_positions, np.ndarray)
        self.assertIsInstance(basis.elements, np.ndarray)

    def test_set_species(self):
        pos, cell = generate_fcc_lattice()
        pse = PeriodicTable()
        el = pse.element("Pt")
        basis = Atoms(symbols='Al', positions=pos, cell=cell)
        self.assertEqual(basis.get_chemical_formula(), "Al")
        basis.set_species([el])
        self.assertEqual(basis.get_chemical_formula(), "Pt")
        self.assertTrue("Al" not in [sp.Abbreviation] for sp in basis._species_to_index_dict.keys())
        self.assertTrue("Pt" in [sp.Abbreviation] for sp in basis._species_to_index_dict.keys())

    def test_new_array(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols='Al', positions=pos, cell=cell)
        basis.set_repeat([10, 10, 10])
        spins = np.ones(len(basis))
        basis.new_array(name="spins", a=spins)
        self.assertTrue(np.array_equal(basis.arrays['spins'], spins))

    def test_set_array(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols='Al', positions=pos, cell=cell)
        basis.set_repeat([10, 10, 10])
        spins = np.ones(len(basis), dtype=float)
        basis.set_array(name="spins", a=2*spins, dtype=int)
        self.assertTrue(np.array_equal(basis.arrays['spins'], 2 * spins))

    def test_get_array(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols='Al', positions=pos, cell=cell)
        basis.set_repeat([10, 10, 10])
        spins = np.ones(len(basis), dtype=float)
        basis.set_array(name="spins", a=2*spins, dtype=int)
        self.assertTrue(np.array_equal(basis.arrays['spins'], 2 * spins))
        self.assertTrue(np.array_equal(basis.get_array(name="spins"), 2 * spins))

    def test_add_tags(self):
        self.CO2.add_tag(test_tag="a")
        self.assertIsInstance(self.CO2.test_tag, SparseList)
        self.assertEqual(self.CO2.test_tag[0], "a")
        self.assertEqual(self.CO2.test_tag[0], self.CO2.test_tag[2])
        self.assertIsInstance(self.CO2.test_tag.list(), list)
        self.CO2.add_tag(selective_dynamics=[True, True, True])
        self.CO2.selective_dynamics[1] = [True, False, True]
        self.assertEqual(self.CO2.selective_dynamics[1], [True, False, True])
        self.assertIsInstance(self.CO2.selective_dynamics.list(), list)

    def test_get_tags(self):
        self.CO2.add_tag(test_tag="a")
        self.assertIsInstance(self.CO2.test_tag, SparseList)
        self.assertIsInstance(self.CO2.get_tags(), type(dict().keys()))

    def test_get_pbc(self):
        self.assertTrue(np.array_equal(self.CO2.pbc, self.CO2.get_pbc()))
        self.assertEqual(len(self.CO2.get_pbc()), 3)

    def test_set_pbc(self):
        self.CO2.set_pbc(value=[True, True, False])
        self.assertTrue(np.array_equal(self.CO2.pbc, self.CO2.get_pbc()))
        self.assertTrue(np.array_equal([True, True, False], self.CO2.get_pbc()))
        self.CO2.set_pbc(value=False)
        self.assertTrue(np.array_equal([False, False, False], self.CO2.get_pbc()))
        self.assertTrue(np.array_equal(self.CO2.pbc, self.CO2.get_pbc()))

    def test_chemical_element(self):
        self.assertIsInstance(self.CO2.convert_element('C'), ChemicalElement)
        self.assertEqual(len(self.CO2.species), 2)

    def test_copy(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols='Al', positions=pos, cell=cell)
        basis_copy = basis.copy()
        self.assertEqual(basis, basis_copy)
        basis_copy[:] = "Pt"
        self.assertNotEqual(basis, basis_copy)

    def test_numbers_to_elements(self):
        num_list = [1, 12, 13, 6]
        self.assertTrue(np.array_equal([el.Abbreviation for el in self.CO2.numbers_to_elements(num_list)],
                                       ['H', 'Mg', 'Al', 'C']))

    def test_to_hdf(self):
        if sys.version_info[0] >= 3:
            filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static/pyiron_atomistics/test_hdf")
            abs_filename = os.path.abspath(filename)
            hdf_obj = FileHDFio(abs_filename)
            pos, cell = generate_fcc_lattice()
            basis = Atoms(symbols='Al', positions=pos, cell=cell)
            basis.set_repeat([2, 2, 2])
            basis.to_hdf(hdf_obj, "test_structure")
            self.assertTrue(np.array_equal(hdf_obj["test_structure/positions"], basis.positions))
            basis_new = Atoms().from_hdf(hdf_obj, "test_structure")
            self.assertEqual(basis, basis_new)

    def test_from_hdf(self):
        if sys.version_info[0] >= 3:
            filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static/pyiron_atomistics/test_hdf")
            abs_filename = os.path.abspath(filename)
            hdf_obj = FileHDFio(abs_filename)
            pos, cell = generate_fcc_lattice()
            basis_store = Atoms(symbols='Al', positions=pos, cell=cell)
            basis_store.set_repeat([2, 2, 2])
            basis_store.to_hdf(hdf_obj, "simple_structure")
            basis = Atoms().from_hdf(hdf_obj, group_name="simple_structure")
            self.assertEqual(len(basis), 8)
            self.assertEqual(basis.get_majority_species()[1], "Al")
            self.assertEqual(basis.get_spacegroup()['Number'], 225)

    def create_Fe_bcc(self):
        self.pse = PeriodicTable()
        self.pse.add_element("Fe", "Fe_up", spin="up", pseudo_name='GGA')
        self.pse.add_element("Fe", "Fe_down", spin="down", pseudo_name='GGA')
        Fe_up = self.pse.element("Fe_up")
        Fe_down = self.pse.element("Fe_down")
        self.Fe_bcc = Atoms([Fe_up, Fe_down], scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]], cell=np.identity(3))
        self.Fe_bcc.add_tag("group")
        self.Fe_bcc.group[:] = 0

    def test_convert_formula(self):
        self.assertEqual(self.CO2.convert_formula('C'), ['C'])
        self.assertEqual(self.CO2.convert_formula('C3'), ['C', 'C', 'C'])
        self.assertEqual(self.CO2.convert_formula('CO2'), ['C', 'O', 'O'])
        self.assertEqual(self.CO2.convert_formula('CO2Fe'), ['C', 'O', 'O', 'Fe'])
        self.assertEqual(self.CO2.convert_formula('CO2FeF21'), ['C', 'O', 'O', 'Fe', 'F', 'F'])

    def test__getitem__(self):
        self.assertEqual(self.CO2[0].symbol, 'C')
        self.assertEqual(self.C3[2].position.tolist(), [0, 2, 0])
        self.assertTrue((self.C3[1:].positions == np.array([[0, 0, 2], [0, 2, 0]])).all())
        short_basis = self.CO2[0]
        self.assertIsInstance(short_basis, Atom)
        short_basis = self.CO2[[0]]
        self.assertIsInstance(short_basis, Atoms)
        self.assertEqual(short_basis.indices[0], 0)
        self.assertEqual(len(short_basis.species), 1)
        short_basis = self.CO2[[2]]
        self.assertIsInstance(short_basis, Atoms)
        self.assertEqual(short_basis.indices[0], 0)
        self.assertEqual(len(short_basis.species), 1)
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        basis_O.positions += [0., 0., 0.5]
        basis = basis_Mg + basis_O
        basis.center_coordinates_in_unit_cell()
        basis.set_repeat([3, 3, 3])
        mg_indices = basis.select_index("Mg")
        o_indices = basis.select_index("O")
        basis_new = basis[mg_indices] + basis[o_indices]
        self.assertEqual(len(basis_new._tag_list), len(basis[mg_indices]) + len(basis[o_indices]))
        self.assertEqual(basis_new.get_spacegroup()["Number"], 225)

    def test_positions(self):
        self.assertEqual(self.CO2[1:].positions[1:].tolist(), [[0.0, 1.5, 0.0]])
        self.CO2.positions[1][0] = 5.
        self.assertEqual(self.CO2.positions[1].tolist(), [5.0, 0, 1.5])

    def test_set_positions(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols='Al', positions=pos, cell=cell)
        basis.set_positions(np.array([[2.5, 2.5, 2.5]]))
        self.assertTrue(np.array_equal(basis.positions, [[2.5, 2.5, 2.5]]))

    def test_set_scaled_positions(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols='Al', positions=pos, cell=cell, a=4.2)
        basis.set_scaled_positions(np.array([[0.5, 0.5, 0.5]]))
        self.assertTrue(np.array_equal(basis.scaled_positions, [[0.5, 0.5, 0.5]]))
        self.assertTrue(np.array_equal(basis.positions, np.dot([[0.5, 0.5, 0.5]], basis.cell)))

    def test_cell(self):
        CO = Atoms("CO",
                   positions=[[0, 0, 0], [0, 0, 2]],
                   cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   pbc=[True, True, True])
        self.assertTrue((CO.get_cell() == np.identity(3)).all())
        self.assertTrue((CO.cell == np.identity(3)).all())
        CO.cell[2][2] = 10.
        self.assertTrue(CO.cell[2, 2] == 10.)

    def test_add(self):
        COX = self.C2 + Atom("O", position=[0, 0, -2])
        COX += Atom("O", position=[0, 0, -4])
        COX += COX
        n_objects = len(set(COX.get_species_objects()))
        n_species = len(set(COX.get_chemical_elements()))
        self.assertEqual(n_objects, n_species)

    def test_pbc(self):
        CO = Atoms("CO",
                   positions=[[0, 0, 0], [0, 0, 2]],
                   cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   pbc=[True, True, True])
        self.assertTrue((CO.pbc == np.array([True, True, True])).all())
        CO.set_pbc((True, True, False))

    def test_get_masses_DOF(self):
        self.assertEqual(len(self.CO2.get_masses_dof()), len(self.CO2.positions.flatten()))

    def test_get_parent_basis(self):
        periodic_table = PeriodicTable()
        periodic_table.add_element(parent_element="O", new_element="O_up")
        O_up = periodic_table.element("O_up")

        O_basis = Atoms([O_up], cell=10.0 * np.eye(3), scaled_positions=[[0.5, 0.5, 0.5]])
        O_simple = Atoms(["O"], cell=10.0 * np.eye(3), scaled_positions=[[0.5, 0.5, 0.5]])
        O_parent = O_basis.get_parent_basis()
        self.assertNotEqual(O_basis, O_parent)
        self.assertEqual(O_simple, O_parent)
        self.assertEqual(O_parent[0].symbol, "O")
        periodic_table.add_element(parent_element="O", new_element="O_down")
        O_down = periodic_table.element("O_down")
        O_basis = Atoms([O_up, O_down], cell=10.0 * np.eye(3), scaled_positions=[[0.5, 0.5, 0.5], [0, 0, 0]])
        O_simple = Atoms(["O", "O"], cell=10.0 * np.eye(3), scaled_positions=[[0.5, 0.5, 0.5]])
        O_parent = O_basis.get_parent_basis()
        self.assertNotEqual(O_basis, O_parent)
        self.assertEqual(O_simple, O_parent)
        self.assertEqual(O_parent.get_chemical_formula(), "O2")
        self.assertEqual(len(O_basis.species), 2)
        self.assertEqual(len(O_simple.species), 1)
        self.assertEqual(len(O_parent.species), 1)

    def test_profiling(self):
        num = 1000
        C100 = Atoms(num * ["C"], positions=[(0, 0, 0) for _ in range(num)])
        self.assertEqual(len(C100), num)

    def test_Au(self):
        a = 4.05  # Gold lattice constant
        b = a / 2.
        fcc = Atoms(['Au'],
                    cell=[(0, b, b), (b, 0, b), (b, b, 0)],
                    pbc=True)
        # print fcc
        # print "volume: ", fcc.get_volume()

    def test_set_absolute(self):
        a = 4.05  # Gold lattice constant
        b = a / 2.
        positions = np.array([(0.5, 0.4, 0.)])
        fcc = Atoms(symbols=['Au'],
                    scaled_positions=positions,
                    cell=[(0, b, b), (b, 0, b), (b, b, 0)],
                    pbc=True)
        # fcc.set_absolute()
        # print fcc.positions
        # fcc.set_relative()
        self.assertTrue(np.linalg.norm(fcc.scaled_positions - positions) < 1e-10)

    def test_repeat(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        basis_O.positions += [0., 0., 0.5]
        basis = basis_Mg + basis_O
        basis.center_coordinates_in_unit_cell()
        basis.set_repeat([3, 3, 3])
        self.assertEqual(basis.get_spacegroup()["Number"], 225)

    def test_boundary(self):
        cell = 2.2 * np.identity(3)
        NaCl = Atoms('NaCl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        NaCl.set_repeat([3, 3, 3])
        # NaCl.plot3d()
        NaCl_bound = NaCl.get_boundary_region(0.2)
        # NaCl_bound.plot3d()

    def test_get_distance(self):
        cell = 2.2 * np.identity(3)
        NaCl = Atoms('NaCl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        self.assertAlmostEqual(NaCl.get_distance(0, 1), 2.2*0.5*np.sqrt(3))
        self.assertAlmostEqual(NaCl.get_distance(0, [0, 0, 0.5]), 0.5)
        self.assertAlmostEqual(NaCl.get_distance([0, 0, 0], [0, 0, 0.5]), 0.5)

    def test_get_neighbors(self):
        cell = 2.2 * np.identity(3)
        NaCl = Atoms('NaCl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        # NaCl.repeat([3, 3, 3])
        # NaCl.positions = [(1,1,1)]
        boundary = NaCl.get_boundary_region(3.5)
        extended_cell = NaCl + boundary
        # extended_cell.plot3d()
        nbr_dict = NaCl.get_neighbors(num_neighbors=12, t_vec=True)
        # print nbr_dict.distances
        # print [set(s) for s in nbr_dict.shells]

    def test_center_coordinates(self):
        cell = 2.2 * np.identity(3)
        NaCl = Atoms('NaCl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        NaCl.set_repeat([3, 3, 3])
        NaCl.positions += [2.2, 2.2, 2.2]
        NaCl.center_coordinates_in_unit_cell(origin=-0.5)
        self.assertTrue(-0.5 < np.min(NaCl.scaled_positions))
        self.assertTrue(np.max(NaCl.scaled_positions < 0.5))
        NaCl.center_coordinates_in_unit_cell(origin=0.)
        self.assertTrue(0 <= np.min(NaCl.positions))
        self.assertTrue(np.max(NaCl.scaled_positions < 1))

    def test_get_shells(self):
        dim = 3
        cell = 2.2 * np.identity(dim)
        Al_sc = Atoms('AlAl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([3, 3, 3])
        self.assertEqual(np.round(Al_sc.get_shells()[2], 6), 2.2)

    def test_cluster_analysis(self):
        import random
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms(elements=['Al', 'Al'], scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([4, 4, 4])
        radius = Al_sc.get_shell_radius()
        neighbors = Al_sc.get_neighbors(radius=radius, num_neighbors=100, t_vec=False, exclude_self=True)

        c_Zn = 0.1
        pse = PeriodicTable()
        Zn = pse.element("Zn")
        random.seed(123456)
        for _ in range(1):
            Zn_ind = random.sample(range(len(Al_sc)), int(c_Zn * len(Al_sc)))
            # for i_Zn in Zn_ind:
            #     Al_sc.elements[i_Zn] = Zn

            cluster = Al_sc.cluster_analysis(Zn_ind, neighbors)
            cluster_len = np.sort([len(v) for k, v in cluster.items()])
            # print np.histogram(cluster_len), np.sum(cluster_len), len(Zn_ind)
            # for key, value in cluster.items():
            #     el = pse.Element((key % 100) + 1)
            #     for i_el in value:
            #         Al_sc.elements[i_el] = el
            # Al_sc.plot3d()

    def test_get_bonds(self):
        dim = 3
        cell = 2.62 * np.identity(dim)
        d1, d2 = 0.6, 0.6
        H2O = Atoms('H2O', scaled_positions=[(d1, d2, 0), (d1, -d2, 0), (0, 0, 0)], cell=cell)
        H2O.set_repeat([1, 1, 3])
        # H2O.plot3d(show_bonds=True) #, bond_stretch=2)
        # print H2O.get_bonds(radius=2.)[0]
        # print np.sum(H2O.get_masses())/H2O.get_volume()

    # def test_get_symmetry(self):
    #     cell = 2.2 * np.identity(3)
    #     Al_sc = Atoms('AlAl', positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
    #     Al_sc.repeat([2, 2, 2])
    #     self.assertEqual(len(Al_sc.get_symmetry()['translations']), 768)

    def test_get_parent_symbols(self):
        self.assertTrue(np.array_equal(self.CO2.get_parent_symbols(), ["C", "O", "O"]))
        self.assertTrue(np.array_equal(self.CO2.get_parent_symbols(), self.CO2.get_chemical_symbols()))
        cell = np.eye(3) * 10.0
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis = Atoms([o_up], scaled_positions=[[0.27, 0.27, 0.27]], cell=cell)
        self.assertTrue(np.array_equal(basis.get_parent_symbols(), ["O"]))
        self.assertFalse(np.array_equal(basis.get_parent_symbols(), basis.get_chemical_symbols()))

    def test_get_chemical_symbols(self):
        self.assertTrue(np.array_equal(self.CO2.get_chemical_symbols(), ["C", "O", "O"]))
        cell = np.eye(3) * 10.0
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis = Atoms([o_up], scaled_positions=[[0.27, 0.27, 0.27]], cell=cell)
        self.assertTrue(np.array_equal(basis.get_chemical_symbols(), ["O_up"]))

    def test_get_symmetry_dataset(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms('AlAl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([2, 2, 2])
        self.assertEqual(Al_sc.get_symmetry_dataset()['number'], 229)

    def test_get_space_group(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms('AlAl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        self.assertEqual(Al_sc.get_spacegroup()['InternationalTableSymbol'], 'Im-3m')
        self.assertEqual(Al_sc.get_spacegroup()['Number'], 229)
        cell = 4.2 * (0.5 * np.ones((3, 3)) - 0.5 * np.eye(3))
        Al_fcc = Atoms('Al', scaled_positions=[(0, 0, 0)], cell=cell)
        self.assertEqual(Al_fcc.get_spacegroup()['InternationalTableSymbol'], 'Fm-3m')
        self.assertEqual(Al_fcc.get_spacegroup()['Number'], 225)
        a = 3.18
        c = 1.623 * a
        cell = np.eye(3)
        cell[0, 0] = a
        cell[2, 2] = c
        cell[1, 0] = -a/2.
        cell[1, 1] = np.sqrt(3) * a / 2.
        pos = np.array([[0., 0., 0.], [1./3., 2./3., 1./2.]])
        Mg_hcp = Atoms('Mg2', scaled_positions=pos, cell=cell)
        self.assertEqual(Mg_hcp.get_spacegroup()['Number'], 194)
        cell = np.eye(3)
        cell[0, 0] = a
        cell[2, 2] = c
        cell[1, 1] = np.sqrt(3) * a
        pos = np.array([[0., 0., 0.], [0.5, 0.5, 0.], [0.5, 0.16666667, 0.5], [0., 0.66666667, 0.5]])
        Mg_hcp = Atoms('Mg4', scaled_positions=pos, cell=cell)
        self.assertEqual(Mg_hcp.get_spacegroup()['Number'], 194)

    def test_get_primitive_cell(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms('AlFe', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([2, 2, 2])
        primitive_cell = Al_sc.get_primitive_cell()
        self.assertEqual(primitive_cell.get_spacegroup()['Number'], 221)

    def test_get_ir_reciprocal_mesh(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms('AlAl', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        self.assertEqual(len(Al_sc.get_ir_reciprocal_mesh([3, 3, 3])[0]), 27)

    def test_get_number_species_atoms(self):
        self.assertEqual(list(self.CO2.get_number_species_atoms().values()), [1, 2])

    def test_get_chemical_formula(self):
        self.assertEqual(self.CO2.get_chemical_formula(), "CO2")

    def test_get_equivalent_atoms(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms('AlFe', scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([2, 2, 2])

    def test_center(self):
        old_pos = self.CO2.positions.copy()
        self.CO2.center(vacuum=5)
        new_array = old_pos + 5 * np.ones(3)
        self.assertTrue(np.array_equal(self.CO2.positions, new_array))

    def test_get_positions(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        self.assertTrue(np.array_equal(basis_Mg.positions, basis_Mg.get_positions()))

    def test_get_scaled_positions(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        self.assertTrue(np.array_equal(basis_Mg.scaled_positions, basis_Mg.get_scaled_positions()))

    def test_occupy_lattice(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        basis_O.scaled_positions += [0., 0., 0.5]
        basis = basis_Mg + basis_O
        basis.center_coordinates_in_unit_cell()
        orig_basis = basis.copy()
        self.assertEqual(basis.get_chemical_formula(), "Mg4O4")
        Mg_indices = basis.select_index("Mg")
        O_indices = basis.select_index("O")
        basis.occupy_lattice(Na=Mg_indices)
        self.assertEqual(basis.get_chemical_formula(), "Na4O4")
        basis.occupy_lattice(Cl=O_indices)
        self.assertEqual(basis.get_chemical_formula(), "Cl4Na4")
        self.assertTrue(np.array_equal(basis.select_index("Na"), Mg_indices))
        self.assertTrue(np.array_equal(basis.select_index("Cl"), O_indices))
        orig_basis.set_repeat([2, 2, 2])
        Mg_indices = orig_basis.select_index("Mg")
        O_indices = orig_basis.select_index("O")
        orig_basis.occupy_lattice(Cl=O_indices, Na=Mg_indices)
        self.assertEqual(orig_basis.get_chemical_formula(), "Cl32Na32")
        orig_basis.occupy_lattice(H=O_indices[0])
        self.assertEqual(orig_basis.get_chemical_formula(), "Cl31HNa32")

    def test_select_index(self):
        self.assertTrue(np.array_equal(self.CO2.select_index("C"), [0]))
        self.assertTrue(np.array_equal(self.CO2.select_index("O"), [1, 2]))

    def test_parent_index(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        basis_O.positions += [0., 0., 0.5]
        basis = basis_Mg + basis_O
        basis.center_coordinates_in_unit_cell()
        basis.set_repeat([2, 2, 2])
        o_indices = basis.select_index("O")
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis[o_indices] = o_up
        self.assertTrue(np.array_equal(o_indices, basis.select_index(o_up)))
        self.assertEqual(len(basis.select_index("O")), 0)
        self.assertTrue(np.array_equal(o_indices, basis.select_parent_index("O")))

    def test__eq__(self):
        test_basis = self.CO2.copy()
        self.assertEqual(test_basis, self.CO2)
        test_basis.positions[2] += 0.0
        self.assertEqual(test_basis, self.CO2)
        self.assertNotEqual(self.C2, self.CO2)

    def test__add__(self):
        cell = np.eye(3) * 10.0
        basis_0 = Atoms(["O"], scaled_positions=[[0.5, 0.5, 0.5]], cell=cell)
        basis_1 = Atoms(["H"], scaled_positions=[[0.75, 0.75, 0.75]], cell=cell)
        basis_2 = Atoms(["H"], scaled_positions=[[0.25, 0.25, 0.25]], cell=cell)
        basis_3 = Atoms(["H", "O", "N"], scaled_positions=[[0.35, 0.35, 0.35], [0., 0., 0.], [0., 0., 0.1]], cell=cell)

        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis_4 = Atoms([o_up], scaled_positions=[[0.27, 0.27, 0.27]], cell=np.eye(3) * 20.0)
        b = basis_0 + basis_1
        self.assertEqual(b.get_chemical_formula(), "HO")
        b = basis_0 + basis_1 + basis_2
        self.assertEqual(b.get_chemical_formula(), "H2O")
        b += basis_2
        self.assertEqual(b.get_chemical_formula(), "H3O")
        b = basis_0 + basis_1 + basis_2 + basis_3
        self.assertEqual(b.get_chemical_formula(), "H3NO2")
        self.assertTrue(np.array_equal(b.scaled_positions[b.select_index("N")], [[0., 0., 0.1]]))
        self.assertTrue(np.allclose(b.scaled_positions[b.select_index("H")], [[0.75, 0.75, 0.75], [0.25, 0.25, 0.25],
                                                                              [0.35, 0.35, 0.35]]))
        self.assertTrue(np.allclose(b.scaled_positions[b.select_index("O")], [[0.5, 0.5, 0.5], [0., 0., 0.]]))
        b.set_repeat([2, 2, 2])
        self.assertEqual(b.get_chemical_formula(), "H24N8O16")
        b += basis_4
        self.assertEqual(b.get_chemical_formula(), "H24N8O16O_up")
        self.assertTrue(np.allclose(b.scaled_positions[b.select_index(o_up)], [[0.27, 0.27, 0.27]]))
        COX = self.C2 + Atom("O", position=[0, 0, -2])
        COX += Atom("O", position=[0, 0, -4])
        COX += COX
        n_objects = len(set(COX.get_species_objects()))
        n_species = len(set(COX.get_chemical_elements()))
        self.assertEqual(n_objects, n_species)
        self.assertEqual(n_objects, 2)
        self.assertEqual(n_species, 2)
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        # basis_O.set_relative()
        basis_O.scaled_positions += [0., 0., 0.5]
        basis = basis_Mg + basis_O
        self.assertEqual(len(basis._tag_list), len(basis_Mg._tag_list) + len(basis_O._tag_list))
        basis.center_coordinates_in_unit_cell()
        self.assertEqual(basis.get_spacegroup()["Number"], 225)

    def test__delitem__(self):
        cell = np.eye(3) * 10.0
        basis_0 = Atoms(["O"], scaled_positions=[[0.5, 0.5, 0.5]], cell=cell)
        basis_1 = Atoms(["H"], scaled_positions=[[0.75, 0.75, 0.75]], cell=cell)
        basis_2 = Atoms(["H"], scaled_positions=[[0.25, 0.25, 0.25]], cell=cell)
        basis_3 = Atoms(["H", "O", "N"], scaled_positions=[[0.35, 0.35, 0.35], [0., 0., 0.], [0., 0., 0.1]], cell=cell)

        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis_4 = Atoms([o_up], scaled_positions=[[0.27, 0.27, 0.27]], cell=cell)
        b = basis_0 + basis_1 + basis_2 + basis_3 + basis_4
        O_indices = b.select_index("O")
        self.assertEqual(len(b), 7)
        self.assertEqual(len(b.indices), 7)
        self.assertEqual(len(b.species), 4)
        b.__delitem__(O_indices[0])
        self.assertEqual(b.get_chemical_formula(), "H3NOO_up")
        self.assertEqual(len(b), 6)
        self.assertEqual(len(b.indices), 6)
        self.assertEqual(len(b._tag_list), 6)
        self.assertEqual(len(b.species), 4)
        O_indices = b.select_index("O")
        b.__delitem__(O_indices)
        self.assertEqual(b.get_chemical_formula(), "H3NO_up")
        self.assertEqual(len(b), 5)
        self.assertEqual(len(b.indices), 5)
        self.assertEqual(len(b.species), 3)
        self.assertEqual(np.max(b.indices), 2)
        N_indices = b.select_index("N")
        b.__delitem__(N_indices)
        self.assertEqual(b.get_chemical_formula(), "H3O_up")
        self.assertEqual(len(b), 4)
        self.assertEqual(len(b.indices), 4)
        self.assertEqual(len(b.species), 2)
        self.assertEqual(np.max(b.indices), 1)
        O_indices = b.select_index(o_up)
        b.__delitem__(O_indices)
        self.assertEqual(b.get_chemical_formula(), "H3")
        self.assertEqual(len(b), 3)
        self.assertEqual(len(b.indices), 3)
        self.assertEqual(len(b.species), 1)
        self.assertEqual(np.max(b.indices), 0)

    def test__setitem__(self):
        basis = self.CO2.copy()
        basis[0] = 'H'
        basis[1] = 'H'
        self.assertEqual(basis.get_chemical_formula(), "H2O")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis = self.CO2.copy()
        basis[0] = 'H'
        basis[np.int64(0)] = 'H'
        self.assertEqual(basis.get_chemical_formula(), "HO2")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis[0] = 'O'
        self.assertEqual(basis.get_chemical_formula(), "O3")
        self.assertEqual(len(basis.species), 1)
        self.assertEqual(len(basis.get_species_symbols()), 1)
        basis = self.CO2.copy()
        basis[[2]] = 'N'
        self.assertEqual(basis.get_chemical_formula(), "CNO")
        self.assertEqual(len(basis.species), 3)
        self.assertEqual(len(basis.get_species_symbols()), 3)
        basis = self.CO2.copy()
        basis[[0]] = 'H'
        basis[np.array([0])] = 'H'
        self.assertEqual(basis.get_chemical_formula(), "HO2")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)

        basis = self.CO2.copy()
        basis[[0]] = 'N'
        self.assertEqual(basis.get_chemical_formula(), "NO2")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis[[0]] = 'O'
        self.assertEqual(basis.get_chemical_formula(), "O3")
        self.assertEqual(len(basis.species), 1)
        self.assertEqual(len(basis.get_species_symbols()), 1)
        basis[[0, 2]] = 'H'
        self.assertEqual(basis.get_chemical_formula(), "H2O")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis[[0, 2]] = o_up
        self.assertEqual(basis.get_chemical_formula(), "OO_up2")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis[0:3] = "N"
        self.assertEqual(basis.get_chemical_formula(), "N3")
        self.assertEqual(len(basis.species), 1)
        self.assertEqual(len(basis.get_species_symbols()), 1)
        basis[:] = "Ne"
        self.assertEqual(basis.get_chemical_formula(), "Ne3")
        self.assertEqual(len(basis.species), 1)
        self.assertEqual(len(basis.get_species_symbols()), 1)
        basis[-2:] = "H"
        self.assertEqual(basis.get_chemical_formula(), "H2Ne")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis[0:3] = "O"
        self.assertEqual(basis.get_chemical_formula(), "O3")
        self.assertEqual(len(basis.species), 1)
        self.assertEqual(len(basis.get_species_symbols()), 1)


def generate_fcc_lattice(a=4.2):
    positions = [[0, 0, 0]]
    cell = (np.ones((3, 3)) - np.eye(3)) * 0.5 * a
    return positions, cell


if __name__ == '__main__':
    unittest.main()
