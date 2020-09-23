# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
import numpy as np
from pyiron_base import Project
from pyiron.atomistics.structure.atoms import Atoms


class TestHessianJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(
            os.path.join(cls.file_location, "hessian_class")
        )
        cls.project.remove_jobs_silently(recursive=True)
        cls.job = cls.project.create_job(
            "HessianJob", "job_test_hessian"
        )
        structure = Atoms(
            positions=[[0, 0, 0], [1, 1, 1]], elements=["Fe", "Fe"], cell=2 * np.eye(3)
        )
        cls.job.set_reference_structure(structure)
        cls.job.structure.apply_strain(0.01)
        cls.job.structure.positions[0, 0] = 0.1
        cls.job.structure.center_coordinates_in_unit_cell()
        cls.job.set_force_constants(force_constants=1)
        cls.job.set_elastic_moduli(bulk_modulus=1, shear_modulus=1)
        cls.job.server.run_mode.interactive = True
        cls.job.run()
        cls.job.structure.positions[0, 1] -= 0.1
        cls.job.run()

    @classmethod
    def tearDownClass(cls):
        cls.job.interactive_close()
        cls.job.remove()
        cls.project.remove(enable=True)

    def test_forces(self):
        self.assertAlmostEqual(self.job.output.forces[0, 0, 0], -0.1)
        self.assertAlmostEqual(self.job.output.forces[0, 0, 1], 0)
        self.assertAlmostEqual(self.job.output.forces[1, 0, 1], 0.1)
        self.assertAlmostEqual(self.job.output.pressures[0, 0, 0], -0.03)


if __name__ == "__main__":
    unittest.main()
