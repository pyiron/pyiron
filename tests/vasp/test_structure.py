import unittest
from pyiron_base.core.settings.generic import Settings
import os

s = Settings(config={'file': 'structure.db',
                     'top_level_dirs': os.path.abspath(os.getcwd()),
                     'resource_paths': os.path.join(os.path.abspath(os.getcwd()), '../static')})

import posixpath
import numpy as np

from pyiron_atomistics.structure.atoms import Atoms
from pyiron_vasp.structure import read_atoms, write_poscar, vasp_sorter, \
    atoms_from_string

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class TestVaspStructure(unittest.TestCase):

    """
    Testing routines in the vasp/structure module.
    """

    def setUp(self):
        poscar_directory = "../static/vasp_test_files/poscar_samples"
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
                self.assertFalse(np.array_equal(atoms.selective_dynamics[0], [True, True, True]))
                self.assertTrue(np.array_equal(atoms.selective_dynamics[0], [False, False, False]))
                self.assertTrue(np.array_equal(atoms.selective_dynamics[-5], [True, True, True]))
                self.assertTrue(np.array_equal(atoms.selective_dynamics[-4], [False, False, False]))

    def test_write_poscar(self):
        write_poscar(structure=self.structure, filename="POSCAR_test")
        test_atoms = read_atoms("POSCAR_test")
        self.assertEqual(self.structure.get_chemical_formula(), test_atoms.get_chemical_formula())
        struct = self.structure.copy()
        struct.add_tag(selective_dynamics=[True, True, True])
        write_poscar(structure=struct, filename="POSCAR_test")
        test_atoms = read_atoms("POSCAR_test")
        truth_array = np.empty_like(struct.positions, dtype=bool)
        truth_array[:] = [True, True, True]
        self.assertTrue(np.array_equal(np.array(test_atoms.selective_dynamics), truth_array))

    def test_vasp_sorter(self):
        write_poscar(structure=self.structure, filename="POSCAR_test")
        test_atoms = read_atoms("POSCAR_test")
        vasp_order = vasp_sorter(self.structure)
        self.assertEqual(len(self.structure), len(test_atoms))
        self.assertEqual(self.structure[vasp_order], test_atoms)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
