# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron.base.project.generic import Project
from pyiron.base.job.wrapper import job_wrapper_function


class TestJobWrapperFunction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "wrapper_testing"))
        cls.job = cls.project.create_job("ExampleJob", "job_test_run")
        cls.job.input["count"] = 12
        cls.job.server.run_mode.manual = True
        cls.job.run()

    @classmethod
    def tearDownClass(cls):
        project = Project(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "wrapper_testing")
        )
        job = project.load(project.get_job_ids()[0])
        job.remove()
        project.remove(enable=True)

    def test_job_wrapper_function(self):
        job_wrapper_function(
            working_directory=os.path.dirname(self.job.project_hdf5.file_name),
            job_id=None,
            file_path=self.job.project_hdf5.file_name + self.job.project_hdf5.h5_path,
            debug=False
        )
        self.assertEqual(self.job['status'], 'finished')
