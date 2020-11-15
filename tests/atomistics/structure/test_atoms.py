# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
import os
import warnings
from pyiron.atomistics.structure.atom import Atom
from pyiron.atomistics.structure.atoms import Atoms, CrystalStructure
from pyiron.atomistics.structure.generator import create_ase_bulk, create_surface, create_hkl_surface, create_structure
from pyiron.atomistics.structure.sparse_list import SparseList
from pyiron.atomistics.structure.periodic_table import PeriodicTable, ChemicalElement
from pyiron_base import FileHDFio, ProjectHDFio, Project
from ase.cell import Cell as ASECell
from ase.atoms import Atoms as ASEAtoms


class TestAtoms(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        if os.path.isfile(
            os.path.join(file_location, "../../static/atomistics/test_hdf")
        ):
            os.remove(
                os.path.join(file_location, "../../static/atomistics/test_hdf")
            )
        if os.path.isfile(
            os.path.join(file_location, "../../static/atomistics/test.h5")
        ):
            os.remove(
                os.path.join(file_location, "../../static/atomistics/test.h5")
            )

    @classmethod
    def setUpClass(cls):
        C = Atom("C").element
        cls.C3 = Atoms([C, C, C], positions=[[0, 0, 0], [0, 0, 2], [0, 2, 0]])
        cls.C2 = Atoms(2 * [Atom("C")])

    def setUp(self):
        # These atoms are reset before every test.
        self.CO2 = Atoms("CO2", positions=[[0, 0, 0], [0, 0, 1.5], [0, 1.5, 0]])

    def test__init__(self):
        pos, cell = generate_fcc_lattice()
        pse = PeriodicTable()
        el = pse.element("Al")
        basis = Atoms()
        ase_basis = ASEAtoms()
        self.assertIsInstance(ase_basis, ASEAtoms)
        self.assertIsInstance(ase_basis.info, dict)
        self.assertIsInstance(ase_basis.arrays, dict)
        self.assertIsInstance(ase_basis.pbc, (bool, list, np.ndarray))
        self.assertIsInstance(ase_basis._cellobj, ASECell)
        self.assertIsInstance(basis, Atoms)
        self.assertIsInstance(basis.info, dict)
        self.assertIsInstance(basis.arrays, dict)
        self.assertIsInstance(basis.units, dict)
        self.assertIsInstance(basis.pbc, (bool, list, np.ndarray))
        self.assertIsInstance(basis.indices, np.ndarray)
        self.assertEqual(len(basis.positions), 0)
        self.assertIsInstance(basis.species, list)
        self.assertIsInstance(basis.elements, np.ndarray)
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        self.assertIsInstance(basis, Atoms)
        self.assertEqual(basis.get_spacegroup()["Number"], 225)
        basis = Atoms(elements="Al", positions=pos, cell=cell)
        self.assertIsInstance(basis, Atoms)
        basis = Atoms(elements=["Al"], positions=pos, cell=cell)
        self.assertIsInstance(basis, Atoms)
        self.assertRaises(
            ValueError, Atoms, symbols="Pt", elements="Al", positions=pos, cell=cell
        )
        basis = Atoms(numbers=[13], positions=pos, cell=cell)
        self.assertEqual(basis.get_majority_species()["symbol"], "Al")
        basis = Atoms(species=[el], indices=[0], positions=pos, cell=cell)
        self.assertEqual(basis.get_majority_species()["symbol"], "Al")
        self.assertIsInstance(basis, Atoms)
        self.assertIsInstance(basis.info, dict)
        self.assertIsInstance(basis.arrays, dict)
        self.assertIsInstance(basis.units, dict)
        self.assertIsInstance(basis.pbc, (bool, list, np.ndarray))
        self.assertIsInstance(basis.indices, np.ndarray)
        self.assertIsInstance(basis.species, list)
        self.assertIsInstance(basis.cell, ASECell)
        self.assertIsInstance(basis.positions, np.ndarray)
        self.assertIsInstance(basis.get_scaled_positions(), np.ndarray)
        self.assertIsInstance(basis.elements, np.ndarray)

    def test_set_species(self):
        pos, cell = generate_fcc_lattice()
        pse = PeriodicTable()
        el = pse.element("Pt")
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        self.assertEqual(basis.get_chemical_formula(), "Al")
        basis.set_species([el])
        self.assertEqual(basis.get_chemical_formula(), "Pt")
        self.assertTrue(
            "Al" not in [sp.Abbreviation] for sp in basis._species_to_index_dict.keys()
        )
        self.assertTrue(
            "Pt" in [sp.Abbreviation] for sp in basis._species_to_index_dict.keys()
        )

    def test_new_array(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        basis.set_repeat([10, 10, 10])
        spins = np.ones(len(basis))
        basis.new_array(name="spins", a=spins)
        self.assertTrue(np.array_equal(basis.arrays["spins"], spins))

    def test_set_array(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        basis.set_repeat([10, 10, 10])
        spins = np.ones(len(basis), dtype=float)
        basis.set_array(name="spins", a=2 * spins, dtype=int)
        self.assertTrue(np.array_equal(basis.arrays["spins"], 2 * spins))

    def test_get_array(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        basis.set_repeat([10, 10, 10])
        spins = np.ones(len(basis), dtype=float)
        basis.set_array(name="spins", a=2 * spins, dtype=int)
        self.assertTrue(np.array_equal(basis.arrays["spins"], 2 * spins))
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
        self.CO2.set_pbc([True, True, False])
        self.assertTrue(np.array_equal(self.CO2.pbc, self.CO2.get_pbc()))
        self.assertTrue(np.array_equal([True, True, False], self.CO2.get_pbc()))
        self.CO2.set_pbc(False)
        self.assertTrue(np.array_equal([False, False, False], self.CO2.get_pbc()))
        self.assertTrue(np.array_equal(self.CO2.pbc, self.CO2.get_pbc()))

    def test_chemical_element(self):
        conv = self.CO2.convert_element("C")
        self.assertIsInstance(conv, ChemicalElement)
        self.assertIsInstance(self.CO2.convert_element(conv), ChemicalElement)
        self.assertIsInstance(self.CO2.convert_element(self.CO2[0]), ChemicalElement)
        with self.assertRaises(ValueError):
            self.assertIsInstance(self.CO2.convert_element(self.CO2), ChemicalElement)
        self.assertEqual(len(self.CO2.species), 2)

    def test_copy(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        basis_copy = basis.copy()
        basis_copy.positions[0, 0] += 5
        self.assertNotEqual(basis_copy.positions[0, 0], basis.positions[0, 0])
        basis_copy.cell[0, 0] += 5
        self.assertNotEqual(basis_copy.cell[0, 0], basis.cell[0, 0])
        basis_copy = basis.copy()
        self.assertEqual(basis, basis_copy)
        basis_copy[:] = "Pt"
        self.assertNotEqual(basis, basis_copy)

    def test_numbers_to_elements(self):
        num_list = [1, 12, 13, 6]
        self.assertTrue(
            np.array_equal(
                [el.Abbreviation for el in self.CO2.numbers_to_elements(num_list)],
                ["H", "Mg", "Al", "C"],
            )
        )

    def test_scaled_pos_xyz(self):
        basis = Atoms(symbols="AlAl", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3))
        pos_xyz = basis.pos_xyz()
        self.assertAlmostEqual(np.linalg.norm(pos_xyz[0] - np.array([0, 1])), 0)
        scaled_pos_xyz = basis.scaled_pos_xyz()
        self.assertAlmostEqual(
            np.linalg.norm(pos_xyz[0] - basis.cell[0, 0] * scaled_pos_xyz[0]), 0
        )

    def test_to_hdf(self):
        filename = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "../../static/atomistics/test_hdf",
        )
        abs_filename = os.path.abspath(filename)
        hdf_obj = FileHDFio(abs_filename)
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        basis.set_repeat([2, 2, 2])
        basis.to_hdf(hdf_obj, "test_structure")
        self.assertTrue(
            np.array_equal(hdf_obj["test_structure/positions"], basis.positions)
        )
        basis_new = Atoms().from_hdf(hdf_obj, "test_structure")
        self.assertEqual(basis, basis_new)

    def test_from_hdf(self):
        filename = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "../../static/atomistics/test_hdf",
        )
        abs_filename = os.path.abspath(filename)
        hdf_obj = FileHDFio(abs_filename)
        pos, cell = generate_fcc_lattice()
        basis_store = Atoms(symbols="Al", positions=pos, cell=cell)
        basis_store.set_repeat([2, 2, 2])
        basis_store.add_tag(selective_dynamics=[False, False, False])
        basis_store.selective_dynamics[7] = [True, True, True]
        basis_store.to_hdf(hdf_obj, "simple_structure")
        basis = Atoms().from_hdf(hdf_obj, group_name="simple_structure")
        self.assertEqual(len(basis), 8)
        self.assertEqual(basis.get_majority_species()["symbol"], "Al")
        self.assertEqual(basis.get_spacegroup()["Number"], 225)
        self.assertTrue(basis.selective_dynamics[7][0])
        self.assertFalse(basis.selective_dynamics[0][0])
        basis.add_tag(selective_dynamics=[False, False, False])
        basis.selective_dynamics[6] = [True, True, True]
        self.assertTrue(basis.selective_dynamics[6][0])
        self.assertFalse(basis.selective_dynamics[5][0])

    def test_to_object(self):
        filename = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "../../static/atomistics",
        )
        abs_filename = os.path.abspath(filename)
        hdf_obj = ProjectHDFio(
            project=Project(abs_filename),
            file_name="test.h5"
        )
        pos, cell = generate_fcc_lattice()
        basis_store = Atoms(symbols="Al", positions=pos, cell=cell)
        basis_store.set_repeat([2, 2, 2])
        basis_store.to_hdf(hdf_obj, "simple_structure")
        basis = hdf_obj["simple_structure"].to_object()
        self.assertEqual(len(basis), 8)
        self.assertEqual(basis.get_majority_species()["symbol"], "Al")
        self.assertEqual(basis.get_spacegroup()["Number"], 225)

    def test_create_Fe_bcc(self):
        self.pse = PeriodicTable()
        self.pse.add_element("Fe", "Fe_up", spin="up", pseudo_name="GGA")
        self.pse.add_element("Fe", "Fe_down", spin="down", pseudo_name="GGA")
        Fe_up = self.pse.element("Fe_up")
        Fe_down = self.pse.element("Fe_down")
        self.Fe_bcc = Atoms(
            [Fe_up, Fe_down],
            scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]],
            cell=np.identity(3),
        )
        self.Fe_bcc.add_tag("group")
        self.Fe_bcc.group[:] = 0

    def test_convert_formula(self):
        self.assertEqual(self.CO2.convert_formula("C"), ["C"])
        self.assertEqual(self.CO2.convert_formula("C3"), ["C", "C", "C"])
        self.assertEqual(self.CO2.convert_formula("CO2"), ["C", "O", "O"])
        self.assertEqual(self.CO2.convert_formula("CO2Fe"), ["C", "O", "O", "Fe"])
        self.assertEqual(
            self.CO2.convert_formula("CO2FeF21"), ["C", "O", "O", "Fe", "F", "F"]
        )

    def test__getitem__(self):
        self.assertEqual(self.CO2[0].symbol, "C")
        self.assertEqual(self.C3[2].position.tolist(), [0, 2, 0])
        self.assertTrue(
            (self.C3[1:].positions == np.array([[0, 0, 2], [0, 2, 0]])).all()
        )
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
        basis_O.positions += [0.0, 0.0, 0.5]
        basis = basis_Mg + basis_O
        basis.center_coordinates_in_unit_cell()
        basis.set_repeat([3, 3, 3])
        mg_indices = basis.select_index("Mg")
        o_indices = basis.select_index("O")
        basis_new = basis[mg_indices] + basis[o_indices]
        self.assertEqual(
            len(basis_new._tag_list), len(basis[mg_indices]) + len(basis[o_indices])
        )
        self.assertEqual(basis_new.get_spacegroup()["Number"], 225)
        self.assertEqual(basis[:-3], basis[0:len(basis)-3])
        self.assertEqual(basis.dimension, basis[mg_indices].dimension)
        self.assertTrue(np.array_equal(basis.pbc, basis[mg_indices].pbc))
        self.assertRaises(IndexError, basis_new.__getitem__, [True, True, False])
        self.assertEqual(basis_new, basis_new[[True] * len(basis_new)])
        bool_array = np.array([True] * len(basis_new))
        bool_array[[10, 20, 40]] = False
        self.assertEqual(len(basis_new[bool_array]), len(basis_new) - 3)
        bool_array = np.array([True] * len(basis))
        bool_array[mg_indices] = False
        self.assertEqual(len(basis[bool_array]), len(o_indices))
        self.assertEqual(len(basis[0:10]), 10)
        self.assertEqual(basis[0, 10], basis[[0, 10]])

    def test_positions(self):
        self.assertEqual(self.CO2[1:].positions[1:].tolist(), [[0.0, 1.5, 0.0]])
        self.CO2.positions[1][0] = 5.0
        self.assertEqual(self.CO2.positions[1].tolist(), [5.0, 0, 1.5])

    def test_set_positions(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell)
        basis.set_positions(np.array([[2.5, 2.5, 2.5]]))
        self.assertTrue(np.array_equal(basis.positions, [[2.5, 2.5, 2.5]]))

    def test_set_scaled_positions(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell, a=4.2)
        basis.set_scaled_positions(np.array([[0.5, 0.5, 0.5]]))
        self.assertTrue(np.array_equal(basis.get_scaled_positions(), [[0.5, 0.5, 0.5]]))
        self.assertTrue(
            np.array_equal(basis.positions, np.dot([[0.5, 0.5, 0.5]], basis.cell))
        )
        with warnings.catch_warnings(record=True):
            basis.scaled_positions = np.array([[0.5, 0.5, 0.5]])
            self.assertTrue(np.array_equal(basis.scaled_positions, [[0.5, 0.5, 0.5]]))

    def test_cell(self):
        CO = Atoms(
            "CO",
            positions=[[0, 0, 0], [0, 0, 2]],
            cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            pbc=[True, True, True],
        )
        self.assertTrue((CO.get_cell() == np.identity(3)).all())
        self.assertTrue((CO.cell == np.identity(3)).all())
        cell = CO.cell.copy()
        cell[2][2] = 10.0
        CO.set_cell(cell)
        self.assertEqual(CO.cell[2, 2], 10.0)
        self.assertAlmostEqual(CO.get_volume(), 10)
        self.assertAlmostEqual(CO.get_volume(per_atom=True), 0.5 * 10)
        CO.set_cell(-np.eye(3))
        with self.assertRaises(ValueError):
            CO.set_cell([2, 1])
        dx = 1.0
        r_o = [0, 0, 0]
        r_h1 = [dx, 0, 0]
        r_h2 = [0, dx, 0]
        water = Atoms(elements=['H', 'H', 'O'], positions=[r_h1, r_h2, r_o])
        self.assertEqual(water.center_coordinates_in_unit_cell(), water)
        water.set_cell(np.zeros((3, 3)))
        self.assertTrue(np.array_equal(water.cell, np.zeros((3, 3))))
        self.assertTrue(np.array_equal(water.get_scaled_positions(), water.positions))
        self.assertEqual(water.center_coordinates_in_unit_cell(), water)

    def test_add(self):
        COX = self.C2 + Atom("O", position=[0, 0, -2])
        COX += Atom("O", position=[0, 0, -4])
        COX += COX
        n_objects = len(set(COX.get_species_objects()))
        n_species = len(set(COX.get_chemical_elements()))
        self.assertEqual(n_objects, n_species)

    def test_pbc(self):
        CO = Atoms(
            "CO",
            positions=[[0, 0, 0], [0, 0, 2]],
            cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            pbc=[True, True, True],
        )
        self.assertTrue((CO.pbc == np.array([True, True, True])).all())
        CO.set_pbc((True, True, False))

    def test_get_masses_DOF(self):
        self.assertEqual(
            len(self.CO2.get_masses_dof()), len(self.CO2.positions.flatten())
        )

    def test_get_center_of_mass(self):
        basis = Atoms(
            elements="AlFe", positions=[3 * [0.5], 3 * [1.5]], cell=2 * np.eye(3)
        )
        mass = np.array(basis.get_masses())
        self.assertAlmostEqual(
            (mass[0] * 0.5 + mass[1] * 1.5) / mass.sum(), basis.get_center_of_mass()[0]
        )
        basis.set_repeat(2)
        self.assertAlmostEqual(
            (mass[0] * 0.5 + mass[1] * 1.5) / mass.sum() + 1,
            basis.get_center_of_mass()[0],
        )

    def test_rotate(self):
        unitcell = Atoms(
            elements="AlFe", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3), pbc=True)
        basis = unitcell.copy()
        basis.rotate(a=10.0, v=[0, 0, 0.1 * np.pi])
        self.assertAlmostEqual(np.arccos(basis.positions[1, :2].sum() / 2) * 180 / np.pi, 10.0)
        basis = unitcell.copy()
        basis.rotate(v=[0, 0, 1], a=0.1)
        self.assertAlmostEqual(np.arccos(basis.positions[1, :2].sum() / 2) * 180 / np.pi, 0.1)
        basis = unitcell.copy()
        center_of_mass = basis.get_center_of_mass()
        basis.rotate(v=[0, 0, 0.1 * np.pi], center="com")
        self.assertTrue(np.allclose(basis.get_center_of_mass(), center_of_mass))
        basis = unitcell.copy()
        center_of_positions = basis.positions.mean(axis=0)
        basis.rotate(v=[0, 0, 1], center="cop")
        self.assertTrue(np.allclose(center_of_positions, basis.positions.mean(axis=0)))
        basis = unitcell.copy()
        position = basis.positions[1]
        basis.rotate(v=[0, 0, 1], center="cou")
        self.assertTrue(np.allclose(position, basis.positions[1]))
        basis = unitcell.copy()
        basis.rotate(v=np.random.random(3), rotate_cell=True)
        self.assertAlmostEqual(basis.get_scaled_positions()[1, 0], 0.5)
        basis = unitcell.copy()
        basis.rotate(v=np.random.random(3), index_list=[0])
        self.assertTrue(
            np.allclose(unitcell.positions.flatten(), basis.positions.flatten())
        )

    def test_rotate_euler(self):
        unitcell = Atoms(
            elements="AlFe", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3)
        )
        basis = unitcell.copy()
        basis.rotate_euler(phi=0.1 * np.pi)
        self.assertAlmostEqual(np.arccos(basis.positions[1, :2].sum() / 2) / np.pi, 0.1)
        basis = unitcell.copy()
        center_of_mass = basis.get_center_of_mass()
        basis.rotate_euler(phi=0.1 * np.pi, center="com")
        self.assertTrue(np.allclose(basis.get_center_of_mass(), center_of_mass))
        basis = unitcell.copy()
        center_of_positions = basis.positions.mean(axis=0)
        basis.rotate_euler(phi=0.1 * np.pi, center="cop")
        self.assertTrue(np.allclose(center_of_positions, basis.positions.mean(axis=0)))
        basis = unitcell.copy()
        position = basis.positions[1]
        basis.rotate_euler(phi=0.1 * np.pi, center="cou")
        self.assertTrue(np.allclose(position, basis.positions[1]))

    def test_set_initial_magnetic_moments(self):
        pos, cell = generate_fcc_lattice()
        basis = Atoms(symbols="Al", positions=pos, cell=cell, a=4.2, pbc=True)
        basis *= 2
        basis.set_initial_magnetic_moments(magmoms=np.ones(len(basis)))
        basis = Atoms(symbols="Al", positions=pos, cell=cell, a=4.2, pbc=True)
        basis.set_initial_magnetic_moments(magmoms=np.ones((len(basis), 3)))
        basis = Atoms(symbols="Al", positions=pos, cell=cell, a=4.2, pbc=True)
        basis *= 2
        basis.set_initial_magnetic_moments(magmoms=np.ones(len(basis)))
        self.assertTrue(np.allclose(basis.arrays["initial_magmoms"], np.ones(len(basis))))
        # set new magnetic moments with different shape
        basis.set_initial_magnetic_moments(magmoms=np.ones((len(basis), 3)))
        self.assertTrue(np.allclose(basis.arrays["initial_magmoms"], np.ones((len(basis), 3))))
        with self.assertRaises(ValueError):
            basis.set_initial_magnetic_moments(magmoms=np.ones(4))

    def test_get_parent_basis(self):
        periodic_table = PeriodicTable()
        periodic_table.add_element(parent_element="O", new_element="O_up")
        O_up = periodic_table.element("O_up")

        O_basis = Atoms(
            [O_up], cell=10.0 * np.eye(3), scaled_positions=[[0.5, 0.5, 0.5]]
        )
        O_simple = Atoms(
            ["O"], cell=10.0 * np.eye(3), scaled_positions=[[0.5, 0.5, 0.5]]
        )
        O_parent = O_basis.get_parent_basis()
        self.assertNotEqual(O_basis, O_parent)
        self.assertEqual(O_simple, O_parent)
        self.assertEqual(O_parent[0].symbol, "O")
        periodic_table.add_element(parent_element="O", new_element="O_down")
        O_down = periodic_table.element("O_down")
        O_basis = Atoms(
            [O_up, O_down],
            cell=10.0 * np.eye(3),
            scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        O_simple = Atoms(
            ["O", "O"], cell=10.0 * np.eye(3), scaled_positions=[[0., 0., 0.], [0.5, 0.5, 0.5]]
        )
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
        b = a / 2.0
        fcc = Atoms(["Au"], cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)
        # print fcc
        # print "volume: ", fcc.get_volume()

    def test_set_absolute(self):
        a = 4.05  # Gold lattice constant
        b = a / 2.0
        positions = np.array([(0.5, 0.4, 0.0)])
        fcc = Atoms(
            symbols=["Au"],
            scaled_positions=positions,
            cell=[(0, b, b), (b, 0, b), (b, b, 0)],
            pbc=True,
        )
        # fcc.set_absolute()
        # print fcc.positions
        # fcc.set_relative()
        self.assertTrue(np.linalg.norm(fcc.get_scaled_positions() - positions) < 1e-10)

    def test_set_relative(self):
        lattice = CrystalStructure(
            element="Al", bravais_basis="fcc", lattice_constants=4
        )
        basis_relative = lattice.copy()
        cell = basis_relative.cell.copy()
        cell[0, 0] = 6
        basis_relative.set_cell(cell, True)
        basis_absolute = lattice.copy()
        cell = basis_absolute.cell.copy()
        cell[0, 0] = 6
        basis_absolute.set_cell(cell)
        self.assertAlmostEqual(
            basis_relative.positions[-1, 0] * 1.5, basis_absolute.positions[-1, 0]
        )
        basis = lattice.copy()
        self.assertAlmostEqual(
            basis.get_scaled_positions(wrap=False)[-1, 0],
            basis_relative.get_scaled_positions(wrap=False)[-1, 0],
        )
        cell = basis.cell.copy()
        cell[0, 0] = 6
        basis.set_cell(cell)
        self.assertAlmostEqual(basis.positions[-1, 0], basis_absolute.positions[-1, 0])
        basis = lattice.copy()
        basis_relative = lattice.copy()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            basis_relative.set_relative()
            self.assertEqual(len(w), 1)
            self.assertIsInstance(w[-1].message, DeprecationWarning)
        basis.positions[-1, 0] = 0.5
        basis_relative.positions[-1, 0] = 0.5
        self.assertAlmostEqual(basis.positions[-1, 0], basis_relative.positions[-1, 0])
        basis.set_cell(3 * np.ones(3))
        self.assertAlmostEqual(basis.get_volume(), 27)
        basis.set_cell(np.append(np.ones(3), 90 - np.random.random(3)).flatten())
        self.assertLess(basis.get_volume(), 1)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            basis_absolute.set_absolute()
            self.assertEqual(len(w), 1)
            self.assertIsInstance(w[-1].message, DeprecationWarning)

    def test_repeat(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        basis_O.set_scaled_positions(basis_O.get_scaled_positions() + [0.0, 0.0, 0.5])
        with self.assertRaises(ValueError):
            basis_O.set_repeat(5.0)
        with self.assertRaises(AssertionError):
            basis_O.set_repeat([2, 2])
        basis = basis_Mg + basis_O
        basis.center_coordinates_in_unit_cell()
        basis.add_tag(selective_dynamics=[True, True, True])
        basis.selective_dynamics[basis.select_index("O")] = [False, False, False]
        len_before = len(basis)
        sel_dyn_before = np.array(basis.selective_dynamics.list())
        self.assertTrue(
            np.alltrue(
                np.logical_not(
                    np.alltrue(sel_dyn_before[basis.select_index("O")], axis=1)
                )
            )
        )
        self.assertTrue(
            np.alltrue(np.alltrue(sel_dyn_before[basis.select_index("Mg")], axis=1))
        )
        basis.set_repeat([3, 3, 2])
        sel_dyn_after = np.array(basis.selective_dynamics.list())
        len_after = len(basis)
        self.assertEqual(basis.get_spacegroup()["Number"], 225)
        self.assertEqual(len_before * 18, len_after)
        self.assertEqual(len(sel_dyn_before) * 18, len(sel_dyn_after))
        self.assertTrue(
            np.alltrue(
                np.logical_not(
                    np.alltrue(sel_dyn_after[basis.select_index("O")], axis=1)
                )
            )
        )
        self.assertTrue(
            np.alltrue(np.alltrue(sel_dyn_after[basis.select_index("Mg")], axis=1))
        )
        basis = basis_Mg + basis_O
        basis.add_tag(spin=None)
        basis.spin[basis.select_index("Mg")] = 1
        basis.spin[basis.select_index("O")] = -1
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("Mg")].list(),
                1 * np.ones(len(basis.select_index("Mg"))),
            )
        )
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("O")].list(),
                -1 * np.ones(len(basis.select_index("O"))),
            )
        )
        basis.set_repeat(2)
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("Mg")].list(),
                1 * np.ones(len(basis.select_index("Mg"))),
            )
        )
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("O")].list(),
                -1 * np.ones(len(basis.select_index("O"))),
            )
        )
        basis = basis_Mg + basis_O
        basis.add_tag(spin=None)
        # Indices set as int
        Mg_indices = np.array(basis.select_index("Mg"), dtype=int).tolist()
        for ind in Mg_indices:
            basis.spin[ind] = 1
        O_indices = np.array(basis.select_index("O"), dtype=int).tolist()
        for ind in O_indices:
            basis.spin[ind] = -1
        basis.set_repeat(2)
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("Mg")].list(),
                1 * np.ones(len(basis.select_index("Mg"))),
            )
        )
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("O")].list(),
                -1 * np.ones(len(basis.select_index("O"))),
            )
        )
        # Indices set as numpy.int
        Mg_indices = np.array(basis.select_index("Mg"), dtype=np.int)
        for ind in Mg_indices:
            basis.spin[ind] = 1
        O_indices = np.array(basis.select_index("O"), dtype=np.int)
        for ind in O_indices:
            basis.spin[ind] = -1
        basis.set_repeat(2)
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("Mg")].list(),
                1 * np.ones(len(basis.select_index("Mg"))),
            )
        )
        self.assertTrue(
            np.array_equal(
                basis.spin[basis.select_index("O")].list(),
                -1 * np.ones(len(basis.select_index("O"))),
            )
        )
        self.assertEqual(8 * len(self.CO2), len(self.CO2.repeat(np.int64(2))))

    def test_get_distance(self):
        cell = 2.2 * np.identity(3)
        NaCl = Atoms("NaCl", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        self.assertAlmostEqual(NaCl.get_distance(0, 1), 2.2 * 0.5 * np.sqrt(3))
        self.assertAlmostEqual(NaCl.get_distance(0, [0, 0, 0.5]), 0.5)
        self.assertAlmostEqual(NaCl.get_distance([0, 0, 0], [0, 0, 0.5]), 0.5)

    def test_find_neighbors_by_vector(self):
        basis = Atoms(symbols=2*["Fe"],
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)],
                      cell=np.identity(3),
                      pbc=True)
        id_lst = basis.find_neighbors_by_vector([0, 0, 1],
                                                num_neighbors=14)
        self.assertEqual(len(np.unique(np.unique(id_lst, return_counts=True)[1])), 1)

    def test_get_neighborhood(self):
        basis = Atoms(
            "FeFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=np.identity(3), pbc=True
        )
        neigh = basis.get_neighborhood([0, 0, 0.1])
        self.assertEqual(neigh.distances[0], 0.1)

    def test_get_neighbors_update_vectors(self):
        structure = CrystalStructure(elements='Fe', lattice_constants=1, bravais_basis='bcc', pbc=True)
        neigh = structure.get_neighbors(num_neighbors=8)
        with self.assertRaises(AssertionError):
            neigh.update_vectors()
        structure = CrystalStructure(elements='Fe', lattice_constants=1, bravais_basis='bcc', pbc=True).repeat(2)
        neigh = structure.get_neighbors(num_neighbors=8)
        self.assertAlmostEqual(np.min(neigh.distances), np.sqrt(3)/2)
        structure.positions[0] += 0.01
        neigh.update_vectors()
        self.assertAlmostEqual(np.min(neigh.distances), np.sqrt(3)*0.49)

    def test_get_neighbors(self):
        basis = Atoms(symbols="FeFe", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3), pbc=True)
        neigh = basis.get_neighbors(num_neighbors=58)
        self.assertAlmostEqual(neigh.distances[0][0], np.sqrt(3))
        counts = np.unique(neigh.shells[0], return_counts=True)
        self.assertTrue(np.array_equal(counts[0], np.arange(5)+1))
        self.assertTrue(np.array_equal(counts[1], np.array([ 8,  6, 12, 24,  8])))
        basis = Atoms(symbols="FeFe", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3), pbc=True)
        basis.pbc = np.array([True, False, False])
        neigh = basis.get_neighbors(num_neighbors=10)
        self.assertAlmostEqual(neigh.distances[0][0], np.sqrt(3))
        self.assertAlmostEqual(neigh.distances[0][2], 2)
        self.assertAlmostEqual(neigh.distances[0][4], np.sqrt(11))
        self.assertAlmostEqual(neigh.distances[0][6], 4)
        self.assertAlmostEqual(neigh.distances[0][8], np.sqrt(27))
        basis.pbc = True
        basis.set_repeat(2)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            neigh = basis.get_neighbors(boundary_width_factor=0.1)
            self.assertGreaterEqual(len(w), 1)
        with self.assertRaises(ValueError):
            neigh = basis.get_neighbors(boundary_width_factor=0.001, num_neighbors=100)
        basis = Atoms(symbols="FeFe", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3))
        neigh = basis.get_neighbors(num_neighbors=1)

    def test_get_neighbors_by_distance(self):
        basis = Atoms(symbols="FeFeFe", positions=[3 * [0], 3 * [1], [0, 0, 1]], cell=2 * np.eye(3), pbc=True)
        neigh = basis.get_neighbors_by_distance(1.5)
        self.assertEqual(len(neigh.distances[0]), 2)
        self.assertEqual(len(neigh.distances[1]), 4)
        self.assertEqual(len(neigh.distances[2]), 6)
        self.assertEqual(neigh.distances[0][0], 1.)
        self.assertAlmostEqual(neigh.distances[1][0], np.sqrt(2))
        self.assertEqual(neigh.distances[2][0], 1.)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            basis.get_neighbors_by_distance(cutoff_radius=1.5, num_neighbors_estimate_buffer=0)
            self.assertGreaterEqual(len(w), 1)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            neigh = basis.get_neighbors_by_distance(cutoff_radius=1.5, num_neighbors=5)
            self.assertGreaterEqual(len(w), 1)
        self.assertEqual(len(neigh.distances[2]), 5)
        with self.assertRaises(ValueError):
            basis.get_neighbors_by_distance(num_neighbors_estimate_buffer=-1)
        # Check with large cell with few atoms
        dx = 0.7
        r_O = [0, 0, 0]
        r_H1 = [dx, dx, 0]
        r_H2 = [-dx, dx, 0]
        unit_cell = 10 * np.eye(3)
        water = Atoms(elements=['H', 'H', 'O'], positions=[r_H1, r_H2, r_O], cell=unit_cell, pbc=True)
        self.assertIsInstance(water.get_neighbors_by_distance(1.3).indices, list)
        water_new = water[[0, 1]]
        self.assertTrue(np.array_equal(water_new.get_neighbors_by_distance(1.3).indices, [np.array([]), np.array([])]))

    def test_get_number_of_neighbors_in_sphere(self):
        basis = Atoms(symbols="FeFeFe", positions=[3 * [0], 3 * [1], [0, 0, 1]], cell=2 * np.eye(3), pbc=True)
        num_neighbors_per_atom = basis.get_numbers_of_neighbors_in_sphere(cutoff_radius=2,
                                                                          num_neighbors_estimate_buffer=0)
        self.assertEqual(num_neighbors_per_atom[0], 10)
        self.assertEqual(num_neighbors_per_atom[1], 12)
        self.assertEqual(num_neighbors_per_atom[2], 6)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            num_neighbors_per_atom = basis.get_numbers_of_neighbors_in_sphere(cutoff_radius=1.5, num_neighbors=5)
            self.assertGreaterEqual(len(w), 1)
        self.assertEqual(num_neighbors_per_atom[2], 5)
        with self.assertRaises(ValueError):
            basis.get_numbers_of_neighbors_in_sphere(num_neighbors_estimate_buffer=-1)

    def test_get_shell_matrix(self):
        structure = CrystalStructure(elements='Fe', lattice_constants=2.83, bravais_basis='bcc')
        shell_mat_atoms = structure.get_shell_matrix(num_neighbors=8)
        neigh = structure.get_neighbors(num_neighbors=8)
        self.assertEqual(shell_mat_atoms[0].sum(), neigh.get_shell_matrix()[0].sum())

    def test_center_coordinates(self):
        cell = 2.2 * np.identity(3)
        NaCl = Atoms("NaCl", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        NaCl.set_repeat([3, 3, 3])
        NaCl.positions += [2.2, 2.2, 2.2]
        NaCl.center_coordinates_in_unit_cell(origin=-0.5)
        self.assertTrue(-0.5 <= np.min(NaCl.get_scaled_positions()))
        self.assertTrue(np.max(NaCl.get_scaled_positions() < 0.5))
        NaCl.center_coordinates_in_unit_cell(origin=0.0)
        self.assertTrue(0 <= np.min(NaCl.positions))
        self.assertTrue(np.max(NaCl.get_scaled_positions() < 1))

    @unittest.skip("skip ovito because it is not installed in the test environment")
    def test_analyse_ovito_cna_adaptive(self):
        basis = Atoms(
            "FeFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=np.identity(3)
        )
        basis.analyse_ovito_cna_adaptive()["CommonNeighborAnalysis.counts.BCC"] == 2

    @unittest.skip("skip ovito because it is not installed in the test environment")
    def test_analyse_ovito_centro_symmetry(self):
        basis = Atoms(
            "FeFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=np.identity(3)
        )
        self.assertTrue(
            all(basis.analyse_ovito_centro_symmetry() == np.array([0.75, 0.75]))
        )

    @unittest.skip("skip ovito because it is not installed in the test environment")
    def test_analyse_ovito_voronoi_volume(self):
        basis = Atoms(
            "FeFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=np.identity(3)
        )
        self.assertTrue(
            all(basis.analyse_ovito_centro_symmetry() == np.array([0.5, 0.5]))
        )

    @unittest.skip("skip nglview because it is not installed in the test environment")
    def test_plot3d(self):
        basis = Atoms(
            "FeFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=np.identity(3)
        )
        view = basis.plot3d()

    @staticmethod
    def test_plot3d_plotly():
        basis = Atoms(
            "FeFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=np.identity(3)
        )
        basis.plot3d(mode='plotly')

    def test_group_points_by_symmetry(self):
        basis = Atoms("FeFe", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3))
        self.assertEqual(len(basis.group_points_by_symmetry([3 * [0.5], 3 * [1.5]])), 1)
        self.assertEqual(len(basis.group_points_by_symmetry([3 * [0.5], 3 * [1.4]])), 2)

    def test_get_equivalent_voronoi_vertices(self):
        basis = Atoms("FeFe", positions=[3 * [0], 3 * [1]], cell=2 * np.eye(3), pbc=True)
        vert = basis.get_equivalent_voronoi_vertices()
        self.assertEqual(len(vert), 1)
        self.assertGreater(
            np.min(np.linalg.norm(vert[0] - basis.positions[0], axis=-1)), 0.5
        )
        self.assertGreater(
            np.min(np.linalg.norm(vert[0] - basis.positions[1], axis=-1)), 0.5
        )

    def test_find_mic(self):
        cell = 0.1*(np.random.random((3,3))-0.5)+np.eye(3)
        basis = Atoms("Fe", positions=[3*[0.5]], cell=cell, pbc=True)
        v = 2*np.random.random(3)-1
        r = np.linalg.inv(cell.T).dot(v)
        r -= np.rint(r)
        self.assertTrue(np.isclose(
            r[0]*cell[0]+r[1]*cell[1]+r[2]*cell[2],
            basis.find_mic(v, vectors=True)
        ).all())
        for v in [np.ones(3), np.ones((3,3)), np.ones((3,3,3))]:
            self.assertTrue(np.array_equal(basis.find_mic(v).shape, v.shape))

    def test_get_distances_array(self):
        basis = Atoms("FeFe", positions=[3*[0], 3*[0.9]], cell=np.identity(3), pbc=True)
        self.assertAlmostEqual(basis.get_distances_array(mic=False)[0, 1], 0.9*np.sqrt(3))
        self.assertTrue(np.allclose(basis.get_distances_array(p1=0.5*np.ones(3)),
                                    basis.get_distances_array(p2=0.5*np.ones(3))))
        self.assertTrue(np.allclose(basis.get_distances_array(vectors=True)[0, 1], -0.1*np.ones(3)))

    def test_repeat_points(self):
        basis = Atoms("Fe", positions=np.random.rand(3).reshape(-1, 3), cell=np.identity(3))
        basis.cell[0, 1] = 0.01
        with self.assertRaises(ValueError):
            basis.repeat_points([0, 0, 0], [2 ,2])
        with self.assertRaises(ValueError):
            basis.repeat_points([0, 0], 2)
        v = np.random.rand(3)
        w = basis.repeat_points(v, 3)
        v += np.array([1, 0.01, 0])
        self.assertAlmostEqual(np.linalg.norm(w-v, axis=-1).min(), 0)
        v = np.random.rand(6).reshape(-1, 3)
        self.assertEqual(basis.repeat_points(v, 2).shape, (8, 2, 3))

    def test_get_extended_positions(self):
        basis = Atoms("FeFe", positions=[[0.01, 0, 0], [0.5, 0.5, 0.5]], cell=np.identity(3), pbc=True)
        with self.assertRaises(ValueError):
            basis.get_extended_positions(-0.1)
        self.assertTrue(np.array_equal(basis.get_extended_positions(0), basis.positions))

    def test_get_equivalent_points(self):
        basis = Atoms("FeFe", positions=[[0.01, 0, 0], [0.5, 0.5, 0.5]], cell=np.identity(3))
        arr = basis.get_equivalent_points([0, 0, 0.5])
        self.assertAlmostEqual(np.linalg.norm(arr-np.array([0.51, 0.5, 0]), axis=-1).min(), 0)

    def test_cluster_analysis(self):
        basis = CrystalStructure("Al", bravais_basis="fcc", lattice_constants=4.2).repeat(10)
        key, counts = basis.cluster_analysis(id_list=[0,1], return_cluster_sizes=True)
        self.assertTrue(np.array_equal(key[1], [0,1]))
        self.assertEqual(counts[0], 2)
        key, counts = basis.cluster_analysis(id_list=[0,int(len(basis)/2)], return_cluster_sizes=True)
        self.assertTrue(np.array_equal(key[1], [0]))
        self.assertEqual(counts[0], 1)

    def test_get_bonds(self):
        basis = CrystalStructure("Al", bravais_basis="fcc", lattice_constants=4.2).repeat(5)
        bonds = basis.get_bonds()
        neigh = basis.get_neighbors()
        self.assertTrue(np.array_equal(np.sort(bonds[0]['Al'][0]),
                        np.sort(neigh.indices[0, neigh.shells[0]==1])))


    def test_get_symmetr(self):
        cell = 2.2 * np.identity(3)
        Al = Atoms("AlAl", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        with self.assertRaises(ValueError):
            Al.symmetrize_vectors(1)
        v = np.random.rand(6).reshape(-1, 3)
        self.assertAlmostEqual(np.linalg.norm(Al.symmetrize_vectors(v)), 0)
        Al.positions[0,0] += 0.01
        w = Al.symmetrize_vectors(v, force_update=True)
        self.assertAlmostEqual(np.absolute(w[:,0]).sum(), np.linalg.norm(w, axis=-1).sum())

    def test_get_symmetry(self):
        cell = 2.2 * np.identity(3)
        Al = Atoms("AlAl", positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell).repeat(2)
        self.assertEqual(len(set(Al.get_symmetry()["equivalent_atoms"])), 1)
        self.assertEqual(len(Al.get_symmetry()["translations"]), 96)
        self.assertEqual(
            len(Al.get_symmetry()["translations"]), len(Al.get_symmetry()["rotations"])
        )

    def test_get_voronoi_vertices(self):
        cell = 2.2 * np.identity(3)
        Al = Atoms("AlAl", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell, pbc=True)
        pos, box = Al._get_voronoi_vertices()
        self.assertEqual(len(pos), 14)

    def test_get_parent_symbols(self):
        self.assertTrue(np.array_equal(self.CO2.get_parent_symbols(), ["C", "O", "O"]))
        self.assertTrue(
            np.array_equal(
                self.CO2.get_parent_symbols(), self.CO2.get_chemical_symbols()
            )
        )
        cell = np.eye(3) * 10.0
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis = Atoms([o_up], scaled_positions=[[0.27, 0.27, 0.27]], cell=cell)
        self.assertTrue(np.array_equal(basis.get_parent_symbols(), ["O"]))
        self.assertFalse(
            np.array_equal(basis.get_parent_symbols(), basis.get_chemical_symbols())
        )

    def test_get_chemical_symbols(self):
        self.assertTrue(
            np.array_equal(self.CO2.get_chemical_symbols(), ["C", "O", "O"])
        )
        cell = np.eye(3) * 10.0
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis = Atoms([o_up], scaled_positions=[[0.27, 0.27, 0.27]], cell=cell)
        self.assertTrue(np.array_equal(basis.get_chemical_symbols(), ["O_up"]))

    def test_get_symmetry_dataset(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms("AlAl", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([2, 2, 2])
        self.assertEqual(Al_sc.get_symmetry_dataset()["number"], 229)

    def test_get_space_group(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms("AlAl", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        self.assertEqual(Al_sc.get_spacegroup()["InternationalTableSymbol"], "Im-3m")
        self.assertEqual(Al_sc.get_spacegroup()["Number"], 229)
        cell = 4.2 * (0.5 * np.ones((3, 3)) - 0.5 * np.eye(3))
        Al_fcc = Atoms("Al", scaled_positions=[(0, 0, 0)], cell=cell)
        self.assertEqual(Al_fcc.get_spacegroup()["InternationalTableSymbol"], "Fm-3m")
        self.assertEqual(Al_fcc.get_spacegroup()["Number"], 225)
        a = 3.18
        c = 1.623 * a
        cell = np.eye(3)
        cell[0, 0] = a
        cell[2, 2] = c
        cell[1, 0] = -a / 2.0
        cell[1, 1] = np.sqrt(3) * a / 2.0
        pos = np.array([[0.0, 0.0, 0.0], [1.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0]])
        Mg_hcp = Atoms("Mg2", scaled_positions=pos, cell=cell)
        self.assertEqual(Mg_hcp.get_spacegroup()["Number"], 194)
        cell = np.eye(3)
        cell[0, 0] = a
        cell[2, 2] = c
        cell[1, 1] = np.sqrt(3) * a
        pos = np.array(
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.5, 0.16666667, 0.5],
                [0.0, 0.66666667, 0.5],
            ]
        )
        Mg_hcp = Atoms("Mg4", scaled_positions=pos, cell=cell)
        self.assertEqual(Mg_hcp.get_spacegroup()["Number"], 194)

    def test_get_primitive_cell(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms("AlFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([2, 2, 2])
        primitive_cell = Al_sc.get_primitive_cell()
        self.assertEqual(primitive_cell.get_spacegroup()["Number"], 221)

    def test_get_ir_reciprocal_mesh(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms("AlAl", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        self.assertEqual(len(Al_sc.get_ir_reciprocal_mesh([3, 3, 3])[0]), 27)

    def test_get_number_species_atoms(self):
        self.assertEqual(list(self.CO2.get_number_species_atoms().values()), [1, 2])

    def test_get_chemical_formula(self):
        self.assertEqual(self.CO2.get_chemical_formula(), "CO2")

    def test_get_equivalent_atoms(self):
        cell = 2.2 * np.identity(3)
        Al_sc = Atoms("AlFe", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)], cell=cell)
        Al_sc.set_repeat([2, 2, 2])

    def test_center(self):
        old_pos = self.CO2.positions.copy()
        self.CO2.center(vacuum=5)
        new_array = old_pos + 5 * np.ones(3)
        self.assertTrue(np.array_equal(self.CO2.positions, new_array))

    def test_get_positions(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constants=4.2)
        self.assertTrue(np.array_equal(basis_Mg.positions, basis_Mg.get_positions()))

    def test_get_scaled_positions(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constants=4.2)
        basis_Mg.set_cell(basis_Mg.cell+0.1 * np.random.random((3, 3)))
        basis_Mg = basis_Mg.center_coordinates_in_unit_cell()
        self.assertTrue(
            np.allclose(
                np.dot(np.linalg.inv(basis_Mg.cell).T, basis_Mg.positions.T).T,
                basis_Mg.get_scaled_positions(),
            )
        )

    def test_apply_strain(self):
        basis_Fe = CrystalStructure("Fe", bravais_basis="bcc", lattice_constants=2.85)
        with self.assertRaises(ValueError):
            basis_Fe.apply_strain(-2)
        basis_new = basis_Fe.apply_strain(0.01, return_box=True)
        self.assertAlmostEqual(basis_new.cell[0,0], 2.85*1.01)
        self.assertAlmostEqual(basis_new.positions[1,0], 0.5*2.85*1.01)
        self.assertAlmostEqual(basis_Fe.cell[0, 0], 2.85)
        basis_Fe.apply_strain(0.01)
        self.assertAlmostEqual(basis_Fe.cell[0,0], 2.85*1.01)
        basis_Fe = CrystalStructure("Fe", bravais_basis="bcc", lattice_constants=2.85)
        basis_Fe.apply_strain(0.01*np.eye(3))
        self.assertAlmostEqual(basis_Fe.cell[0,0], 2.85*1.01)

    def test_get_spherical_coordinates(self):
        basis_Fe = CrystalStructure("Fe", bravais_basis="bcc", lattice_constants=2.85)
        x = basis_Fe.get_spherical_coordinates([0, 0, 0])
        self.assertAlmostEqual(x[0, 0], 0)
        x = basis_Fe.get_spherical_coordinates()
        self.assertAlmostEqual(x[1, 2], 0.25*np.pi)

    def test_occupy_lattice(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        basis_O.set_scaled_positions(basis_O.get_scaled_positions() + [0.0, 0.0, 0.5])
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

    def test_get_majority_species(self):
        basis = Atoms(
            symbols=4 * ["Fe"], positions=np.random.random((4, 3)), cell=np.eye(3)
        )
        self.assertEqual(basis.get_majority_species()["count"], 4)
        self.assertEqual(basis.get_majority_species()["symbol"], "Fe")
        basis = Atoms(
            symbols=["Fe", "Cu", "Ni", "Al"],
            positions=np.random.random((4, 3)),
            cell=np.eye(3),
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            basis.get_majority_species()
            self.assertEqual(len(w), 1)

    def test_select_index(self):
        basis = Atoms(
            symbols=["Fe", "Cu", "Ni", "Al"],
            positions=np.random.random((4, 3)),
            cell=np.eye(3),
        )
        self.assertTrue(np.array_equal(basis.select_index("Fe"), [0]))
        self.assertTrue(np.array_equal(basis.select_index("Ni"), [2]))
        self.assertTrue(np.array_equal(basis.select_index(["Cu", "Al"]), [1, 3]))
        Fe = basis.convert_element("Fe")
        Ni = basis.convert_element("Ni")
        self.assertTrue(np.array_equal(basis.select_index([Fe, Ni]), [0, 2]))
        pse = PeriodicTable()
        pse.add_element("Ni", "Ni_up", spin=1)
        ni_up = pse.element("Ni_up")
        basis = Atoms(
            symbols=["Fe", "Cu", ni_up, "Al"],
            positions=np.random.random((4, 3)),
            cell=np.eye(3),
        )
        self.assertTrue(np.array_equal(basis.select_index("Fe"), [0]))
        self.assertTrue(np.array_equal(basis.select_index(ni_up), [2]))
        self.assertTrue(np.array_equal(basis.select_index(["Cu", "Al"]), [1, 3]))
        Fe = basis.convert_element("Fe")
        Ni = basis.convert_element(ni_up)
        self.assertTrue(np.array_equal(basis.select_index([Fe, Ni]), [0, 2]))

    def test_parent_index(self):
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_O = CrystalStructure("O", bravais_basis="fcc", lattice_constant=4.2)
        basis_O.positions += [0.0, 0.0, 0.5]
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
        basis_3 = Atoms(
            ["H", "O", "N"],
            scaled_positions=[[0.35, 0.35, 0.35], [0.0, 0.0, 0.0], [0.0, 0.0, 0.1]],
            cell=cell,
        )

        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        basis_4 = Atoms(
            [o_up], scaled_positions=[[0.27, 0.27, 0.27]], cell=np.eye(3) * 20.0
        )
        b = basis_0 + basis_1
        self.assertEqual(b.get_chemical_formula(), "HO")
        b = basis_0 + basis_1 + basis_2
        self.assertEqual(b.get_chemical_formula(), "H2O")
        b += basis_2
        self.assertEqual(b.get_chemical_formula(), "H3O")
        b = basis_0 + basis_1 + basis_2 + basis_3
        self.assertEqual(b.get_chemical_formula(), "H3NO2")
        self.assertTrue(
            np.array_equal(
                b.get_scaled_positions()[b.select_index("N")], [[0.0, 0.0, 0.1]]
            )
        )
        self.assertTrue(
            np.allclose(
                b.get_scaled_positions()[b.select_index("H")],
                [[0.75, 0.75, 0.75], [0.25, 0.25, 0.25], [0.35, 0.35, 0.35]],
            )
        )
        self.assertTrue(
            np.allclose(
                b.get_scaled_positions()[b.select_index("O")],
                [[0.5, 0.5, 0.5], [0.0, 0.0, 0.0]],
            )
        )
        b.set_repeat([2, 2, 2])
        self.assertEqual(b.get_chemical_formula(), "H24N8O16")
        b += basis_4
        self.assertEqual(b.get_chemical_formula(), "H24N8O16O_up")
        self.assertTrue(
            np.allclose(
                b.get_scaled_positions()[b.select_index(o_up)], [[0.27, 0.27, 0.27]]
            )
        )
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
        basis_O.set_scaled_positions([0.0, 0.0, 0.5] + basis_O.get_scaled_positions())
        basis = basis_Mg + basis_O
        self.assertEqual(
            len(basis._tag_list), len(basis_Mg._tag_list) + len(basis_O._tag_list)
        )
        basis.center_coordinates_in_unit_cell()
        self.assertEqual(basis.get_spacegroup()["Number"], 225)
        # Adding an ASE instance to a pyiron instance
        ase_basis = ASEAtoms("O", scaled_positions=[[0, 0, 0]], cell=np.eye(3) * 10)
        pyiron_basis = Atoms("O", scaled_positions=[[0.5, 0.5, 0.5]], cell=np.eye(3) * 10)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            pyiron_basis += ase_basis
            self.assertEqual(len(pyiron_basis), 2)
            self.assertEqual(len(ase_basis), 1)
            self.assertIsInstance(pyiron_basis, Atoms)
        ase_basis += pyiron_basis
        self.assertEqual(len(ase_basis), 3)
        self.assertIsInstance(ase_basis, ASEAtoms)
        self.assertNotIsInstance(ase_basis, Atoms)
        self.assertEqual(len(w), 1)
        pyiron_basis += ase_basis[0]
        self.assertEqual(len(pyiron_basis), 3)
        pyiron_basis = Atoms("O", scaled_positions=[[0.5, 0.5, 0.5]], cell=np.eye(3) * 10, pbc=True)
        larger_cell = pyiron_basis.repeat(2)
        larger_cell.positions += 2.5
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            larger_cell += pyiron_basis
            self.assertEqual(len(w), 1)
        basis_1 = Atoms("O", scaled_positions=[[0.5, 0.5, 0.5]], cell=np.eye(3) * 10)
        basis_2 = Atoms("O", scaled_positions=[[0., 0.5, 0.5]], cell=np.eye(3) * 10, pbc=True)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            basis_1 += basis_2
            self.assertEqual(len(w), 1)
        a_0 = 2.86
        structure = create_structure('Fe', 'bcc', a_0)
        carbon = Atoms(symbols=['C'], positions=[[0, 0, 0.5 * a_0]])
        structure += carbon
        self.assertEqual(carbon.indices[0], 0)

    def test_append(self):
        a_0 = 2.86
        structure = create_structure('Fe', 'bcc', a_0)
        carbon = Atoms(symbols=['C'], positions=[[0, 0, 0.5 * a_0]], pbc=True)
        with warnings.catch_warnings(record=True) as w:
            structure.append(carbon)
            self.assertEqual(len(w), 0)
            structure = create_structure('Fe', 'bcc', a_0)
            carbon.cell = np.random.rand(3)
            structure.append(carbon)
            self.assertEqual(len(w), 1)

    def test__delitem__(self):
        cell = np.eye(3) * 10.0
        basis_0 = Atoms(["O"], scaled_positions=[[0.5, 0.5, 0.5]], cell=cell)
        basis_1 = Atoms(["H"], scaled_positions=[[0.75, 0.75, 0.75]], cell=cell)
        basis_2 = Atoms(["H"], scaled_positions=[[0.25, 0.25, 0.25]], cell=cell)
        basis_3 = Atoms(
            ["H", "O", "N"],
            scaled_positions=[[0.35, 0.35, 0.35], [0.0, 0.0, 0.0], [0.0, 0.0, 0.1]],
            cell=cell,
        )

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
        basis[0] = "H"
        basis[1] = "H"
        self.assertEqual(basis.get_chemical_formula(), "H2O")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis = self.CO2.copy()
        basis[0] = "H"
        basis[np.int64(0)] = "H"
        self.assertEqual(basis.get_chemical_formula(), "HO2")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis[0] = "O"
        self.assertEqual(basis.get_chemical_formula(), "O3")
        self.assertEqual(len(basis.species), 1)
        self.assertEqual(len(basis.get_species_symbols()), 1)
        basis = self.CO2.copy()
        basis[[2]] = "N"
        self.assertEqual(basis.get_chemical_formula(), "CNO")
        self.assertEqual(len(basis.species), 3)
        self.assertEqual(len(basis.get_species_symbols()), 3)
        basis = self.CO2.copy()
        basis[[0]] = "H"
        basis[np.array([0])] = "H"
        self.assertEqual(basis.get_chemical_formula(), "HO2")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)

        basis = self.CO2.copy()
        basis[[0]] = "N"
        self.assertEqual(basis.get_chemical_formula(), "NO2")
        self.assertEqual(len(basis.species), 2)
        self.assertEqual(len(basis.get_species_symbols()), 2)
        basis[[0]] = "O"
        self.assertEqual(basis.get_chemical_formula(), "O3")
        self.assertEqual(len(basis.species), 1)
        self.assertEqual(len(basis.get_species_symbols()), 1)
        basis[[0, 2]] = "H"
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
        lat_0 = CrystalStructure(
            "Al", bravais_basis="fcc", lattice_constant=4.0
        ).repeat(3)
        # lat_0.set_SQS(['Al', 'Mg'], x=1/4)  # simple access to SQS
        lat_0[:] = "V"
        self.assertEqual(lat_0.get_chemical_formula(), "V108")
        lat_0[[1, 3, 5]] = "Mg"  # direct occupation
        self.assertEqual(lat_0.get_chemical_formula(), "Mg3V105")
        # lat_0[[0]] = 'V'                     # vacancy (note: do not delete atom)
        lat_1 = lat_0.copy()
        lat_1.set_scaled_positions(1 / 4 + lat_1.get_scaled_positions())
        lat_1[:] = "V"
        self.assertEqual(lat_1.get_chemical_formula(), "V108")
        lat_1[[1, 4, 9]] = "H"
        lat_1[[2, 5, 8]] = "C"
        self.assertEqual(lat_1.get_chemical_formula(), "C3H3V102")
        lat_1.set_scaled_positions(1 / 4 + lat_1.get_scaled_positions())
        lat_1[:] = "V"  # vacancies
        self.assertEqual(lat_1.get_chemical_formula(), "V108")
        basis_Mg = CrystalStructure("Mg", bravais_basis="fcc", lattice_constant=4.2)
        basis_Mg.set_repeat(3)
        basis_Mg[:-3] = "Al"
        self.assertEqual(basis_Mg.get_chemical_formula(), 'Al105Mg3')
        basis_Mg[4:-len(basis_Mg)+7] = "C"
        self.assertEqual(basis_Mg.get_chemical_formula(), 'Al102C3Mg3')
        basis_Mg[4:] = "C"
        self.assertEqual(basis_Mg.get_chemical_formula(), 'Al4C104')
        basis_Mg[:] = "Mg"
        self.assertEqual(basis_Mg.get_chemical_formula(), 'Mg108')
        basis_Mg[::2] = "Al"
        self.assertEqual(basis_Mg.get_chemical_formula(), 'Al54Mg54')
        struct = CrystalStructure("Al", bravais_basis="fcc", lattice_constant=4.2, bravais_lattice="cubic")
        struct[0] = 'Mg'
        self.assertEqual(struct.get_chemical_formula(), 'Al3Mg')
        struct[1] = 'Cu'
        self.assertEqual(struct.get_chemical_formula(), 'Al2CuMg')
        struct[0] = 'Cu'
        struct[1] = 'Cu'
        struct[2] = 'Cu'
        struct[3] = 'Cu'
        self.assertEqual(struct.get_chemical_formula(), 'Cu4')
        struct = CrystalStructure("Al", bravais_basis="fcc", lattice_constant=4.2, bravais_lattice="cubic")
        struct[0] = 'Mg'
        self.assertEqual(struct.get_chemical_formula(), 'Al3Mg')
        struct[1] = 'Cu'
        self.assertEqual(struct.get_chemical_formula(), 'Al2CuMg')
        struct[0:] = 'N'
        self.assertEqual(struct.get_chemical_formula(), 'N4')
        struct = CrystalStructure("Al", bravais_basis="fcc", lattice_constant=4.2, bravais_lattice="cubic")
        struct[0] = 'Mg'
        self.assertEqual(struct.get_chemical_formula(), 'Al3Mg')
        struct[1] = 'Cu'
        self.assertEqual(struct.get_chemical_formula(), 'Al2CuMg')
        struct[0:] = 'Cu'
        self.assertEqual(struct.get_chemical_formula(), 'Cu4')
        struct = CrystalStructure("Al", bravais_basis="fcc", lattice_constant=4.2, bravais_lattice="cubic")
        struct[0] = 'Mg'
        self.assertEqual(struct.get_chemical_formula(), 'Al3Mg')
        struct[1] = 'Cu'
        self.assertEqual(struct.get_chemical_formula(), 'Al2CuMg')
        struct[0:] = 'Mg'
        self.assertEqual(struct.get_chemical_formula(), 'Mg4')

    def test_static_functions(self):
        Al_bulk = create_ase_bulk("Al")
        self.assertIsInstance(Al_bulk, Atoms)
        self.assertTrue(all(Al_bulk.pbc))
        surface = create_surface("Al", "fcc111", size=(4, 4, 4), vacuum=10)
        self.assertTrue(all(surface.pbc))
        surface = create_surface("Al", "fcc111", size=(4, 4, 4), vacuum=10, pbc=[True, True, False])
        self.assertTrue(all(surface.pbc[0:2]))
        self.assertFalse(surface.pbc[2])
        self.assertIsInstance(surface, Atoms)
        hkl_surface = create_hkl_surface(Al_bulk, [10, 8, 7], layers=20, vacuum=10)
        self.assertIsInstance(hkl_surface, Atoms)
        self.assertTrue(all(hkl_surface.pbc))
        hkl_surface_center = create_hkl_surface(
            Al_bulk, [10, 8, 7], layers=20, vacuum=10, center=True
        )
        mean_z = np.mean([p[2] for p in hkl_surface_center.positions])
        self.assertAlmostEqual(mean_z, hkl_surface_center.cell[2][2]/2)

    def test_non_periodic(self):
        structure = CrystalStructure("Fe", bravais_basis="bcc", lattice_constant=4.2)
        pos = structure.repeat([1, 1, 2]).positions.copy()
        structure = CrystalStructure("Fe", bravais_basis="bcc", lattice_constant=4.2)
        structure.pbc = [False, False, True]
        pos_new = structure.repeat([1, 1, 2]).positions.copy()
        self.assertTrue(np.allclose(pos, pos_new))
        c3 = Atoms("C3", positions=[[0, 0, 0], [0, 0, 2], [0, 2, 0]])
        c3.get_scaled_positions()
        c3 = Atoms("C3", positions=[[0, 0, 0], [0, 0, 2], [0, 2, 0]], cell=np.eye(3)*10)
        c3.get_scaled_positions()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            c3.get_scaled_positions()
            self.assertEqual(len(w), 0)

    def test_get_wrapped_coordinates(self):
        structure = CrystalStructure("Fe", bravais_basis="bcc", lattice_constant=4.2, pbc=True)
        position = structure.get_wrapped_coordinates(structure.cell*1.1)
        self.assertAlmostEqual(
            np.linalg.norm(position-structure.cell*0.1), 0
        )


def generate_fcc_lattice(a=4.2):
    positions = [[0, 0, 0]]
    cell = (np.ones((3, 3)) - np.eye(3)) * 0.5 * a
    return positions, cell


if __name__ == "__main__":
    unittest.main()
