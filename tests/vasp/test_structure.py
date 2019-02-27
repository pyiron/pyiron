import unittest
import os
import posixpath
import numpy as np

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.vasp.structure import read_atoms, write_poscar, vasp_sorter, atoms_from_string
from pyiron.atomistics.structure.sparse_list import SparseList
import warnings


class TestVaspStructure(unittest.TestCase):

    """
    Testing routines in the vasp/structure module.
    """

    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        poscar_directory = os.path.join(cls.file_location, "../static/vasp_test_files/poscar_samples")
        file_list = os.listdir(poscar_directory)
        cls.file_list = [posixpath.join(poscar_directory, f) for f in file_list]
        atom_numbers = np.random.randint(low=1, high=99, size=(1, 3)).flatten()
        cell = 10.0 * np.eye(3)
        pos = 0.5 * np.ones((3, 3)) - 0.5 * np.eye(3)
        cls.structure = Atoms(numbers=atom_numbers, cell=cell, positions=pos)
        cls.assertIsInstance(cls.structure, Atoms)
        cls.structure.repeat([2, 2, 2])
        cls.element_list = cls.structure.get_chemical_elements()

    def test_atoms_from_string(self):
        for poscar_file in self.file_list:
            with open(poscar_file, "r") as f:
                lines = f.readlines()
                if poscar_file.split("/")[-1] == "POSCAR_spoilt":
                    self.assertRaises(AssertionError, atoms_from_string, string=lines)
                else:
                    atoms = atoms_from_string(string=lines)
                    self.assertIsInstance(atoms, Atoms)

    def test_read_atoms(self):
        for f in self.file_list:
            if f.split("/")[-1] == "POSCAR_velocity":
                atoms, velocities = read_atoms(filename=f, return_velocities=True)
                self.assertEqual(len(atoms), 19)
                self.assertEqual(np.shape(velocities), (19, 3))
                self.assertEqual(len(atoms.selective_dynamics), 19)
                self.assertEqual(len(atoms.select_index("Mg")), 10)
                self.assertIsInstance(atoms.selective_dynamics, SparseList)
                neon_indices = atoms.select_index("Ne")
                hydrogen_indices = atoms.select_index("H")
                oxygen_indices = atoms.select_index("O")
                truth_array = np.empty_like(atoms.positions[neon_indices], dtype=bool)
                truth_array[:, :] = True
                sel_dyn = np.array(atoms.selective_dynamics.list())
                self.assertTrue(np.array_equal(sel_dyn[neon_indices], np.logical_not(truth_array)))
                truth_array = np.empty_like(atoms.positions[oxygen_indices], dtype=bool)
                truth_array[:, :] = True
                self.assertTrue(np.array_equal(sel_dyn[oxygen_indices], truth_array))
                truth_array = np.empty_like(atoms.positions[hydrogen_indices], dtype=bool)
                truth_array[:, :] = True
                self.assertTrue(np.array_equal(sel_dyn[hydrogen_indices], truth_array))
                velocities_neon = np.zeros_like(np.array(velocities)[neon_indices])
                self.assertTrue(np.array_equal(np.array(velocities)[neon_indices], velocities_neon))

            if f.split("/")[-1] == "POSCAR_no_species":
                atoms = read_atoms(filename=f)
                self.assertEqual(len(atoms), 33)
                self.assertEqual(len(atoms.selective_dynamics), 33)

            elif f.split("/")[-1] != "POSCAR_spoilt":
                atoms = read_atoms(filename=f)
                self.assertIsInstance(atoms, Atoms)
                if f.split("/")[-1] == "POSCAR_1":
                    self.assertEqual(len(atoms), 744)
                    self.assertEqual(len(atoms.select_index("H")), 432)
                    self.assertEqual(len(atoms.select_index("O")), 216)
                    self.assertEqual(len(atoms.select_index("Mg")), 96)
                    with warnings.catch_warnings(record=True) as w:
                        warnings.simplefilter("always")
                        atoms_new, velocities = read_atoms(filename=f, return_velocities=True)
                        self.assertEqual(w[-1].category, UserWarning)
                        warning_string = "The velocities are either not available or they are incomplete/corrupted. " \
                                         "Returning empty list instead"
                        self.assertEqual(str(w[-1].message), warning_string)
                    self.assertEqual(atoms_new, atoms)
                    self.assertEqual(velocities, list())

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
        struct = self.structure.copy()
        struct.add_tag(selective_dynamics=[True, True, True])
        write_poscar(structure=struct, filename=posixpath.join(self.file_location, "POSCAR_test"), cartesian=False)
        test_atoms_new = read_atoms(posixpath.join(self.file_location, "POSCAR_test"))
        self.assertEqual(test_atoms, test_atoms_new)
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
