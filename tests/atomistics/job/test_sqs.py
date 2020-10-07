# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron.project import Project

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
            from sqsgenerator.core.sqs import ParallelSqsIterator

            self.project.create_job(
                'SQSJob', "job_test"
            )
        except ImportError:
            pass

if __name__ == "__main__":
    unittest.main()
