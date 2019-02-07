import unittest
import os
import posixpath
import numpy as np

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.vasp.structure import read_atoms, write_poscar, vasp_sorter, atoms_from_string
from pyiron.atomistics.structure.sparse_list import SparseList


class TestVaspStructure(unittest.TestCase):

    """
    Testing routines in the vasp/structure module.
    """

    def setUp(self):
        self.file_location = os.path.dirname(os.path.abspath(__file__))
        poscar_directory = os.path.join(self.file_location, "../static/vasp_test_files/poscar_samples")
        file_list = os.listdir(poscar_directory)
        self.file_list = [posixpath.join(poscar_directory, f) for f in file_list]
        atom_numbers = np.random.randint(low=1, high=99, size=(1, 3)).flatten()
        cell = 10.0 * np.eye(3)
        pos = 0.5 * np.ones((3, 3)) - 0.5 * np.eye(3)
        self.structure = Atoms(numbers=atom_numbers, cell=cell, positions=pos)
        self.assertIsInstance(self.structure, Atoms)
        self.structure.repeat([2, 2, 2])
        self.element_list = self.structure.get_chemical_elements()

    def test_atoms_from_string(self):
        for poscar_file in self.file_list:
            with open(poscar_file, "r") as f:
                lines = f.readlines()
                atoms = atoms_from_string(string=lines)
                self.assertIsInstance(atoms, Atoms)

    def test_read_atoms(self):
        for f in self.file_list:
            atoms = read_atoms(filename=f)
            self.assertIsInstance(atoms, Atoms)
            if f.split("/")[-1] == "POSCAR_1":
                self.assertEqual(len(atoms), 744)
                self.assertEqual(len(atoms.select_index("H")), 432)
                self.assertEqual(len(atoms.select_index("O")), 216)
                self.assertEqual(len(atoms.select_index("Mg")), 96)
            if f.split("/")[-1] == "POSCAR_scaled":
                self.assertEqual(len(atoms), 256)
                self.assertEqual(len(atoms.select_index("Cu")), 256)
                cell = np.eye(3) * 4. * 3.63
                self.assertTrue(np.array_equal(atoms.cell, cell))
                self.assertEqual(atoms.get_spacegroup()["Number"], 225)
            if f.split("/")[-1] == "POSCAR_volume_scaled":
                self.assertEqual(len(atoms), 256)
                self.assertEqual(len(atoms.select_index("Cu")), 256)
                cell = np.eye(3) * 4. * 3.63
                self.assertTrue(np.array_equal(atoms.cell, cell))
                self.assertEqual(atoms.get_spacegroup()["Number"], 225)
            if f.split("/")[-1] == "POSCAR_random":
                self.assertEqual(len(atoms), 33)
                self.assertEqual(len(atoms.selective_dynamics), 33)
                self.assertEqual(len(atoms.select_index("Zn")), 1)
                self.assertIsInstance(atoms.selective_dynamics, SparseList)
                truth_array = np.empty_like(atoms.positions, dtype=bool)
                truth_array[:] = [True, True, True]
                truth_array[0] = [False, False, False]
                truth_array[-4:] = [False, False, False]
                self.assertTrue(np.array_equal(atoms.selective_dynamics.list(), truth_array))

    def test_write_poscar(self):
        write_poscar(structure=self.structure, filename=posixpath.join(self.file_location, "POSCAR_test"))
        test_atoms = read_atoms(posixpath.join(self.file_location, "POSCAR_test"))
        self.assertEqual(self.structure.get_chemical_formula(), test_atoms.get_chemical_formula())
        struct = self.structure.copy()
        struct.add_tag(selective_dynamics=[True, True, True])
        write_poscar(structure=struct, filename=posixpath.join(self.file_location, "POSCAR_test"))
        test_atoms = read_atoms(posixpath.join(self.file_location, "POSCAR_test"))
        truth_array = np.empty_like(struct.positions, dtype=bool)
        truth_array[:] = [True, True, True]
        self.assertTrue(np.array_equal(np.array(test_atoms.selective_dynamics.list()), truth_array))
        os.remove(posixpath.join(self.file_location, "POSCAR_test"))

    def test_vasp_sorter(self):
        write_poscar(structure=self.structure, filename=posixpath.join(self.file_location, "POSCAR_test"))
        test_atoms = read_atoms(posixpath.join(self.file_location, "POSCAR_test"))
        vasp_order = vasp_sorter(self.structure)
        self.assertEqual(len(self.structure), len(test_atoms))
        self.assertEqual(self.structure[vasp_order], test_atoms)
        os.remove(posixpath.join(self.file_location, "POSCAR_test"))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
