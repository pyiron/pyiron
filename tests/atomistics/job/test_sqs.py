# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron.project import Project
from pyiron.atomistics.structure.atoms import CrystalStructure


class TestSQS(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "test_sqs"))
        cls.project.remove_jobs_silently(recursive=True)

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project.remove(enable=True, enforce=True)

    def test_run(self):
        try:  # Only test if the machine has access to sqsgenerator -- at time of writing Windows doesn't
            job = self.project.create_job(
                'SQSJob', "job_test"
            )
            structure = CrystalStructure("Al", bravais_basis="fcc", lattice_constant=4)
            job.structure = structure
            with self.assertRaises(ValueError):
                job.validate_ready_to_run()
            job.input.mole_fractions = dict()
            structure[0] = 'Sc'
            job.structure = structure
            job.validate_ready_to_run()
            job.input.mole_fractions = dict()
            mole_fractions = {'Al': 0.5, 'Sc': 0.5}
            job.input.mole_fractions = mole_fractions.copy()
            job.structure = structure
            job.validate_ready_to_run()
            self.assertAlmostEqual(job.input.mole_fractions['Al'], mole_fractions['Al'])
        except ImportError:
            pass


if __name__ == "__main__":
    unittest.main()
