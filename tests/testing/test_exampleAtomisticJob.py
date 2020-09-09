# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
import numpy as np
from pyiron_base.project.generic import Project
from pyiron.atomistics.structure.atoms import Atoms


class TestExampleJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.count = 12
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(
            os.path.join(cls.file_location, "random_testing_atomistic")
        )
        cls.job = cls.project.create_job(
            "AtomisticExampleJob", "job_test_atomistic_run"
        )
        cls.job.input["count"] = cls.count
        cls.job.structure = Atoms(
            positions=[[0, 0, 0], [1, 1, 1]], elements=["Fe", "Fe"], cell=2 * np.eye(3)
        )
        cls.job.interactive_open()
        cls.job.run()

    @classmethod
    def tearDownClass(cls):
        cls.job.interactive_close()
        project = Project(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "random_testing_atomistic"
            )
        )
        job = project.load(project.get_job_ids()[0])
        job.remove()
        project.remove(enable=True)

    def test_get_structure(self):
        struct = self.job.get_structure()
        self.assertIsNotNone(struct)
        self.assertEqual(struct.positions.shape, (2, 3))
        self.assertEqual(struct.cell.shape, (3, 3))


if __name__ == "__main__":
    unittest.main()
