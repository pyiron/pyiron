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

    @classmethod
    def tearDownClass(cls):
        cls.project.remove_jobs_silently(recursive=True)
        cls.project.remove(enable=True, enforce=True)

    def test_run(self):
        basis = self.project.create_structure("Fe", "bcc", 2.8)
        job = self.project.create_job( 'HessianJob', "job_single")
        job.server.run_mode.interactive = True
        job.set_reference_structure(basis)
        job.set_force_constants(1)
        job.structure.positions[0,0] += 0.01
        minim = job.create_job("ScipyMinimizer", "job_scipy")
        force_tolerance = 1.0e-3
        minim.input.ionic_force_tolerance = force_tolerance
        minim.run()
        self.assertLess(np.linalg.norm(minim.ref_job['output/generic/forces'][-1], axis=-1).max(), force_tolerance)

    def test_run_pressure(self):
        basis = self.project.create_structure("Al", "fcc", 4)
        job = self.project.create_job( 'HessianJob', "job_pressure")
        job.server.run_mode.interactive = True
        job.set_reference_structure(basis)
        job.set_force_constants(1)
        job.set_elastic_moduli(1,1)
        job.structure.set_cell(job.structure.cell*1.01, scale_atoms=True)
        minim = job.create_job("ScipyMinimizer", "job_scipy_pressure")
        minim.calc_minimize(pressure=0, volume_only=True)
        pressure_tolerance = 1.0e-3
        minim.input.ionic_force_tolerance = pressure_tolerance
        minim.run()
        self.assertLess(np.absolute(minim.ref_job['output/generic/pressures'][-1]).max(), pressure_tolerance)

    def test_calc_minimize(self):
        minim = self.project.create_job('ScipyMinimizer', 'calc_minimize')
        with self.assertRaises(ValueError):
            minim.calc_minimize(volume_only=True, pressure=None)
        minim.calc_minimize(pressure=0)
        self.assertTrue(np.array_equal(minim.input.pressure, np.zeros(1)))
        minim.calc_minimize(pressure=1)
        self.assertTrue(np.array_equal(minim.input.pressure, np.ones(1)))
        minim.calc_minimize(pressure=[1, None, 0])
        self.assertTrue(np.array_equal(minim.input.pressure, np.array([1, None, 0, None, None, None])))

if __name__ == "__main__":
    unittest.main()
