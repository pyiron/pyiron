import os
import unittest
from pyiron.project import Project


class TestExampleJob(unittest.TestCase):
    def setUp(self):
        self.count_run_one = 12
        self.count_run_two = 12
        self.file_location = os.path.dirname(os.path.abspath(__file__))
        self.project = Project(os.path.join(self.file_location, 'random_testing_lib'))
        self.ham = self.project.create_job("ExampleJob", "job_test_run")
        self.ham.input['count'] = self.count_run_one
        self.ham.server.run_mode.interactive = True
        self.ham.run()
        self.ham.input['count'] = self.count_run_two
        self.ham.run()
        self.ham.interactive_close()

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, 'random_testing_lib'))
        ham = project.load(project.get_job_ids()[0])
        ham.remove()
        project.remove(enable=True)

    def test_output(self):
        count = self.ham.get("output/generic/count")
        energy = self.ham.get("output/generic/energy")
        self.assertEqual(self.count_run_one, count[0])
        self.assertEqual(self.count_run_two, count[1])
        self.assertEqual(count[0], len(energy[0]))
        self.assertEqual(count[0], len(energy[1]))


if __name__ == '__main__':
    unittest.main()
