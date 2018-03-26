import unittest
from pyiron_base.core.settings.config.testing import ConfigTesting
from pyiron_base.core.settings.generic import Settings
import os

config = ConfigTesting(sql_lite_database='./testing_jobs.db',
                       path_potentials='../../../static/potentials/',
                       path_project=str(os.getcwd()))
s = Settings(config=config)

from pyiron_base.project import Project


class TestGenericJob(unittest.TestCase):
    def setUp(self):
        self.project = Project('jobs_testing')

    @classmethod
    def tearDownClass(cls):
        project = Project('jobs_testing')
        project.remove(enforce=True)
        s.close_connection()
        os.remove('testing_jobs.db')

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


if __name__ == '__main__':
    unittest.main()
