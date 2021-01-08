# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import posixpath
import unittest
import os
import numpy as np
import scipy.constants
from pyiron_atomistic.atomistics.structure.atoms import Atoms
from pyiron_atomistic.sphinx.structure import read_atoms

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Feb 4, 2018"

BOHR_TO_ANGSTROM = (
    scipy.constants.physical_constants["Bohr radius"][0] / scipy.constants.angstrom
)


class TestSphinxStructure(unittest.TestCase):

    """
    Testing routines in the sphinx/structure module.
    """

    def setUp(self):
        self.file_location = os.path.dirname(os.path.abspath(__file__))
        structure_directory = "../static/sphinx/sphinx_test_files"
        file_list = ["structure_1.sx", "structure_2.sx"]
        self.file_list = [
            posixpath.join(self.file_location, structure_directory, f)
            for f in file_list
        ]
        atom_numbers = np.random.randint(low=1, high=99, size=(1, 3)).flatten()
        cell = 10.0 * np.eye(3)
        pos = 0.5 * np.ones((3, 3)) - 0.5 * np.eye(3)
        self.structure = Atoms(numbers=atom_numbers, cell=cell, positions=pos)
        self.assertIsInstance(self.structure, Atoms)
        self.structure.repeat([2, 2, 2])
        self.element_list = self.structure.get_chemical_elements()

    def test_read_atoms(self):
        for i, f in enumerate(self.file_list):
            atoms = read_atoms(filename=f)
            self.assertIsInstance(atoms, Atoms)
            if i == 0:
                self.assertEqual(atoms.get_chemical_formula(), "Mg5")
                self.assertTrue(
                    np.allclose(
                        atoms.cell / BOHR_TO_ANGSTROM,
                        [
                            [18.0936435257, 0.0, 0.0],
                            [-12.0624290171, 20.8927399203, 0.0],
                            [0.0, 0.0, 39.1932378013],
                        ],
                    )
                )
            if i == 1:
                self.assertEqual(atoms.get_chemical_formula(), "C2Mg5")
                self.assertTrue(
                    np.allclose(
                        atoms.cell / BOHR_TO_ANGSTROM,
                        [
                            [18.09364353, 0.00000000, 0.00000000],
                            [-6.03121451, 20.89273992, -0.00000000],
                            [0.00000000, 0.00000000, 39.19323780],
                        ],
                    )
                )

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
