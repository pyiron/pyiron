import os
import unittest

from pyiron_base.core.settings.config.testing import ConfigTesting
from pyiron_base.core.settings.generic import Settings

config = ConfigTesting(sql_lite_database='./testing_vasp_import.db', path_project=str(os.getcwd()),
                       path_potentials='../../static/potentials/')
s = Settings(config=config)

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
        os.remove('testing_vasp_import.db')

    def test_import(self):
        self.project.import_from_path(path='../vasp_test_files/full_job_sample', recursive=False)
        ham = self.project.load(1)
        self.assertTrue(isinstance(ham, Vasp))


if __name__ == '__main__':
    unittest.main()
