# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
import os
from pyiron.project import Project
from pyiron.vasp.vasp import Vasp


class TestVaspImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "vasp_import_testing"))

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.file_location, "vasp_import_testing"))
        project.remove_jobs_silently(recursive=True)
        project.remove(enable=True)

    def test_import(self):
        folder_path = os.path.join(
                self.file_location, "../static/vasp_test_files/full_job_sample"
            )
        self.project.import_from_path(path=folder_path, recursive=False)
        ham = self.project.load("full_job_sample")
        self.assertTrue(isinstance(ham, Vasp))
        self.assertEqual(ham.get_nelect(), 16)
        self.assertTrue(
            np.array_equal(ham.structure.get_initial_magnetic_moments(), [-1, -1])
        )
        self.assertRaises(IOError, ham.get_final_structure_from_file)
        self.assertIsInstance(ham.output.unwrapped_positions, np.ndarray)


if __name__ == "__main__":
    unittest.main()
