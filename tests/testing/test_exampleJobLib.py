# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron_base import Project


class TestExampleJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.count_run_one = 12
        cls.count_run_two = 12
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "random_testing_lib"))
        cls.ham = cls.project.create_job("ExampleJob", "job_test_run")
        cls.ham.input["count"] = cls.count_run_one
        cls.ham.server.run_mode.interactive = True
        cls.ham.run()
        cls.ham.input["count"] = cls.count_run_two
        cls.ham.run()
        cls.ham.interactive_close()

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, "random_testing_lib"))
        ham = project.load(project.get_job_ids()[0])
        ham.remove()
        project.remove(enable=True)

    def test_output(self):
        count = self.ham.get("output/generic/count")
        energy = self.ham.get("output/generic/energy")
        self.assertEqual(self.count_run_one, count[0])
        self.assertEqual(self.count_run_two, count[1])
        self.assertEqual(count[0], len(energy[0]))
        self.assertEqual(count[0], len(energy[1]))


if __name__ == "__main__":
    unittest.main()
