import os
import unittest
from pyiron.base.project import Project


class TestCopyTo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'testing_copyto'))

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, 'testing_copyto'))
        sub_project = project.open('sub_project_ex')
        sub_project.remove(enable=True)
        sub_project = project.open('sub_project')
        sub_project.remove(enable=True)

    # def test_copy_to_job(self):
    #     job_ser = self.project.create_job("SerialMaster", "sequence_single")
    #     ham = self.project.create_job('ScriptJob', "job_single")
    #     ham.copy_to(job_ser)
    #     self.assertTrue(job_ser['job_single/input/input_inp'])
    #     job_ser.remove()

    def test_copy_to_project(self):
        sub_project = self.project.copy()
        sub_project = sub_project.open('sub_project')
        ham = self.project.create_job('ScriptJob', "job_single_pr")
        ham.copy_to(sub_project)
        os.remove(os.path.join(self.file_location, 'testing_copyto/sub_project/job_single_pr.h5'))

    # def test_copy_to_job_ex(self):
    #     job_ser = self.project.create_job("SerialMaster", "sequence_single_ex")
    #     ham = self.project.create_job('ScriptJob', "job_single_ex")
    #     ham.to_hdf()
    #     ham.copy_to(job_ser)
    #     self.assertTrue(job_ser['job_single_ex/input/input_inp'])
    #     ham.remove()
    #     job_ser.remove()

    def test_copy_to_project_ex(self):
        sub_project = self.project.copy()
        print(sub_project.project_path)
        sub_project = sub_project.open('sub_project_ex')
        ham = self.project.create_job('ScriptJob', "job_single_pr_ex")
        ham.to_hdf()
        ham.copy_to(sub_project)
        ham.remove()
        os.remove(os.path.join(self.file_location, 'testing_copyto/sub_project_ex/job_single_pr_ex.h5'))

if __name__ == '__main__':
    unittest.main()
