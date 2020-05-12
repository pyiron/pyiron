# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron.base.project.generic import Project


class TestGenericJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "jobs_testing"))

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, "jobs_testing"))
        project.remove(enforce=True)

    # def test_generic_jobs(self):
    #     ham = self.project.create_job("ExampleJob", "job_single")
    #     job_ser = self.project.create_job("GenericMaster", "job_list")
    #     job_ser.append(ham)
    #     job_ser.to_hdf()
    #     job_ser_reload = self.project.create_job("GenericMaster", "job_list")
    #     job_ser_reload.from_hdf()
    #     self.assertTrue(job_ser_reload['job_single/input/input_inp'])
    #     job_ser.remove()
    #     ham.remove()
    #
    # def test_generic_jobs_ex(self):
    #     ham = self.project.create_job("ExampleJob", "job_single_ex")
    #     ham.to_hdf()
    #     job_ser = self.project.create_job("GenericMaster", "job_list_ex")
    #     job_ser.append(ham)
    #     job_ser.to_hdf()
    #     self.assertTrue(job_ser['job_single_ex/input/input_inp'])
    #     job_ser_reload = self.project.create_job("GenericMaster", "job_list_ex")
    #     job_ser_reload.from_hdf()
    #     self.assertTrue(job_ser_reload['job_single_ex/input/input_inp'])
    #     job_ser.remove()
    #     ham.remove()


if __name__ == "__main__":
    unittest.main()
