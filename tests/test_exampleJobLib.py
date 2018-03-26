import os
from pyiron_base.core.settings.config.testing import ConfigTesting
from pyiron_base.core.settings.generic import Settings
import unittest

config = ConfigTesting(sql_lite_database='./testing_random_lib.db', path_project=str(os.getcwd()),
                       path_potentials='../../../../static/potentials/')
s = Settings(config=config)

from pyiron.project import Project


class TestExampleJob(unittest.TestCase):
    def setUp(self):
        self.count = 12
        self.project = Project('random_testing_lib')
        self.ham = self.project.create_job("ExampleJob", "job_test_run")
        self.ham.input['count'] = self.count
        self.ham.library_activated = True
        self.ham.run()

    @classmethod
    def tearDownClass(cls):
        project = Project('random_testing_lib')
        ham = project.load(1)
        ham.remove()
        project.remove()
        s.close_connection()
        os.remove('testing_random_lib.db')

    def test_output(self):
        count = self.ham.get("output/generic/count")
        energy_length = len(self.ham.get("output/generic/energy"))
        self.assertEqual(self.count, count)
        self.assertEqual(count, energy_length)


if __name__ == '__main__':
    unittest.main()
