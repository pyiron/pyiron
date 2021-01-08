# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron_atomistic.project import Project


class TestSxExtOptInteractive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "../static/sxextopt"))
        cls.basis = cls.project.create_structure("Fe", "bcc", 2.8)
        job = cls.project.create_job(
            cls.project.job_type.AtomisticExampleJob, "job_single"
        )
        job.server.run_mode.interactive = True
        job.structure = cls.basis
        cls.sxextoptint = cls.project.create_job("SxExtOptInteractive", "job_sxextopt")
        cls.sxextoptint.ref_job = job

    def test_input(self):
        self.assertEqual(self.sxextoptint.input["ionic_steps"], 1000)

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "../static/sxextopt"))
        cls.project.remove_jobs_silently(recursive=True)


if __name__ == "__main__":
    unittest.main()
