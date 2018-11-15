import unittest
import os
from pyiron.project import Project
from pyiron.vasp.vasp import Vasp


class TestVaspImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'vasp_import_testing'))

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.file_location, 'vasp_import_testing'))
        project.remove_jobs(recursive=True)
        project.remove(enable=True)

    def test_import(self):
        self.project.import_from_path(path=os.path.join(self.file_location,
                                                        '../static/vasp_test_files/full_job_sample'),
                                      recursive=False)
        ham = self.project.load(1)
        self.assertTrue(isinstance(ham, Vasp))


if __name__ == '__main__':
    unittest.main()
