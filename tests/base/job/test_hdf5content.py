import unittest
import os
from pyiron.base.project.generic import Project


class DatabasePropertyIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'hdf5_content'))
        cls.ham = cls.project.create_job("ExampleJob", "job_test_run")
        cls.ham.run()

    @classmethod
    def tearDownClass(cls):
        project = Project(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hdf5_content'))
        ham = project.load(project.get_job_ids()[0])
        ham.remove()
        project.remove(enable=True)

    def test_inspect_job(self):
        job_inspect = self.project.inspect(self.ham.job_name)
        self.assertIsNotNone(job_inspect)
        self.assertEqual(job_inspect.content.input.input_inp.data_dict, job_inspect['input/input_inp/data_dict'])
        self.assertTrue(job_inspect.content.output.generic.energy_tot[-1], job_inspect['output/generic/energy_tot'][-1])
        self.assertEqual(job_inspect.content.output.generic.volume[-1], job_inspect['output/generic/volume'][-1])
        self.assertEqual(sorted(dir(job_inspect.content.output.generic)),
                         sorted(job_inspect['output/generic'].list_nodes()))
        self.assertEqual(job_inspect.content.output.__repr__(), job_inspect['output'].__repr__())


if __name__ == '__main__':
    unittest.main()
