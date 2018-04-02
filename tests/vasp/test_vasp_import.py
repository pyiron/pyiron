import unittest
from pyiron_base.core.settings.generic import Settings
import os

s = Settings(config={'file': 'import.db',
                     'top_level_dirs': os.path.normpath(os.path.abspath(os.path.join(os.getcwd(), '..'))),
                     'resource_paths': os.path.join(os.path.abspath(os.getcwd()), '../static')})

from pyiron.project import Project
from pyiron_vasp.vasp import Vasp


class TestVaspImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Project('vasp_import_testing')

    @classmethod
    def tearDownClass(cls):
        project = Project('vasp_import_testing')
        project.remove_jobs(recursive=True)
        project.remove()
        s.close_connection()
        os.remove('import.db')

    def test_import(self):
        self.project.import_from_path(path='../../static/vasp_test_files/full_job_sample', recursive=False)
        ham = self.project.load(1)
        self.assertTrue(isinstance(ham, Vasp))


if __name__ == '__main__':
    unittest.main()
