# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import mock
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
        project.remove_jobs_silently(recursive=True)
        project.remove(enable=True)

    def test_structure_creation(self):
        self.assertIsInstance(self.project.create_ase_bulk("Al"), Atoms)
        surface = self.project.create_surface("Al", "fcc111", size=(4, 4, 4), vacuum=10)
        self.assertTrue(all(surface.pbc))
        self.assertIsInstance(surface, Atoms)
        surface = self.project.create_surface("Al", "fcc111", size=(4, 4, 4), vacuum=10, pbc=[True, True, False])
        self.assertFalse(all(surface.pbc))
        self.assertIsInstance(self.project.create_structure("Al", "fcc", 4.05), Atoms)

    def test_remove_jobs(self):
        sample_job = self.project.create_job("ScriptJob", "Sample")
        sample_job.save()
        with mock.patch('builtins.input', return_value="n"):
            self.project.remove_jobs(recursive=True)
        self.assertEqual(len(self.project.list_nodes()), 1)
        with mock.patch('builtins.input', return_value="y"):
            self.project.remove_jobs(recursive=True)
        self.assertEqual(len(self.project.list_nodes()), 0)

if __name__ == "__main__":
    unittest.main()
