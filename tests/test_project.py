# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron import Project
from pyiron.atomistics.structure.atoms import Atoms


class TestProject(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        # print(cls.execution_path)
        cls.project = Project(os.path.join(cls.execution_path, "test_project"))

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, "test_project"))
        project.remove(enable=True)

    def test_structure_creation(self):
        self.assertIsInstance(self.project.create_ase_bulk("Al"), Atoms)
        surface = self.project.create_surface("Al", "fcc111", size=(4, 4, 4), vacuum=10)
        self.assertTrue(all(surface.pbc))
        self.assertIsInstance(surface, Atoms)
        surface = self.project.create_surface("Al", "fcc111", size=(4, 4, 4), vacuum=10, pbc=[True, True, False])
        self.assertFalse(all(surface.pbc))


if __name__ == "__main__":
    unittest.main()