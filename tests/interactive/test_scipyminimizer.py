# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron.project import Project
import numpy as np


class TestSxExtOptInteractive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "../static/minim"))
        cls.basis = cls.project.create_structure("Fe", "bcc", 2.8)
        job = cls.project.create_job( 'EinsteinExampleJob', "job_single")
        job.server.run_mode.interactive = True
        job.structure = cls.basis
        cls.minim = cls.project.create_job("ScipyMinimizer", "job_scipy")
        cls.minim.ref_job = job

    @classmethod
    def tearDownClass(cls):
        cls.project.remove_jobs_silently(recursive=True)
        cls.project.remove(enable=True, enforce=True)

    def test_run(self):
        self.minim.run()
        self.assertAlmostEqual(np.linalg.norm(self.minim.ref_job['output/generic/forces'][-1], axis=-1).max(), 0)

    def test_minimizer(self):
        self.minim.minimizer = 'CG'
        self.assertEqual(self.minim.minimizer, 'CG')

if __name__ == "__main__":
    unittest.main()
