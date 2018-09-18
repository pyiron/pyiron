import unittest
from pyiron_base.core.settings.generic import Settings
import os

file_location = os.path.dirname(os.path.abspath(__file__))
s = Settings(config={'sql_file': os.path.join(file_location, 'import.db'),
                     'resource_paths': os.path.join(file_location, '../static')})

from pyiron.project import Project
from pyiron_vasp.vasp import Vasp


class TestVaspImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Project(os.path.join(file_location, 'vasp_import_testing'))

    @classmethod
    def tearDownClass(cls):
        project = Project(os.path.join(file_location, 'vasp_import_testing'))
        project.remove_jobs(recursive=True)
        project.remove(enable=True)
        s.close_connection()
        os.remove('import.db')

    def test_import(self):
        self.project.import_from_path(path=os.path.join(file_location, '../static/vasp_test_files/full_job_sample'),
                                      recursive=False)
        ham = self.project.load(1)
        self.assertTrue(isinstance(ham, Vasp))


if __name__ == '__main__':
    unittest.main()
