import os
import unittest
from pyiron_base.core.settings.generic import Settings, convert_path

s = Settings(config={'sql_file': 'library.db',
                     'project_paths': convert_path(os.getcwd()),
                     'resource_paths': convert_path(os.getcwd())})

from pyiron.project import Project


class TestExampleJob(unittest.TestCase):
    def setUp(self):
        self.count_run_one = 12
        self.count_run_two = 12
        self.project = Project('random_testing_lib')
        self.ham = self.project.create_job("ExampleJob", "job_test_run")
        self.ham.input['count'] = self.count_run_one
        self.ham.server.run_mode.interactive = True
        self.ham.run()
        self.ham.input['count'] = self.count_run_two
        self.ham.run()
        self.ham.interactive_close()

    @classmethod
    def tearDownClass(cls):
        project = Project('random_testing_lib')
        ham = project.load(1)
        ham.remove()
        project.remove(enable=True)
        s.close_connection()
        os.remove('library.db')

    def test_output(self):
        count = self.ham.get("output/generic/count")
        energy = self.ham.get("output/generic/energy")
        self.assertEqual(self.count_run_one, count[0])
        self.assertEqual(self.count_run_two, count[1])
        self.assertEqual(count[0], len(energy[0]))
        self.assertEqual(count[0], len(energy[1]))


if __name__ == '__main__':
    unittest.main()
