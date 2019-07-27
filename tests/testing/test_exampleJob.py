import os
import unittest
import numpy as np
from pyiron.base.project.generic import Project


class TestExampleJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.count = 12
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'random_testing'))
        cls.ham = cls.project.create_job("ExampleJob", "job_test_run")
        cls.ham.input['count'] = cls.count
        cls.ham.run()

    @classmethod
    def tearDownClass(cls):
        project = Project(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'random_testing'))
        ham = project.load(project.get_job_ids()[0])
        ham.remove()
        project.remove(enable=True)

    def test_input(self):
        with open(os.path.join(self.file_location, 'random_testing/job_test_run_hdf5/job_test_run/input.inp')) as input_file:
            lines = input_file.readlines()
        input_lst = ["alat 3.2 #lattice constant (would be in a more realistic example in the structure file)\n",
                     "alpha 0.1 #noise amplitude\n",
                     "a_0 3 #equilibrium lattice constant\n",
                     "a_1 0\n",
                     "a_2 1.0 #2nd order in energy (corresponds to bulk modulus)\n",
                     "a_3 0.0 #3rd order\n",
                     "a_4 0.0 #4th order\n",
                     "count " + str(self.count) + " #number of calls (dummy)\n",
                     "write_restart True\n",
                     "read_restart False\n"]
        self.assertEqual(input_lst[0], lines[0])
        self.assertEqual(input_lst[1], lines[1])
        self.assertEqual(input_lst[2], lines[2])
        self.assertEqual(input_lst[3], lines[3])
        self.assertIn(lines[4], [input_lst[4], "a_2 1. #2nd order in energy (corresponds to bulk modulus)\n"])
        self.assertIn(lines[5], [input_lst[5], "a_3 0. #3rd order\n"])
        self.assertIn(lines[6], [input_lst[6], "a_4 0. #4th order\n"])
        self.assertEqual(input_lst[7], lines[7])
        self.assertEqual(input_lst[8], lines[8])
        self.assertEqual(input_lst[9], lines[9])

    def test_restart_file(self):
        with open(os.path.join(self.file_location, 'random_testing/job_test_run_hdf5/job_test_run/restart.out')) as restart_file:
            lines = restart_file.readlines()
        restart_lst = ["count " + str(self.count) + " \n"]
        self.assertEqual(restart_lst, lines)

    def test_output(self):
        energy_lst = self.ham.get("output/generic/energy_tot")
        self.assertEqual(len(energy_lst), self.count)
        with open(os.path.join(self.file_location, 'random_testing/job_test_run_hdf5/job_test_run/output.log')) as output_file:
            lines = output_file.readlines()
        output_lst = ["exampleExecutable logFile \n",
                      "alat 3.2 \n",
                      "count " + str(self.count) + " \n"]
        self.assertEqual(output_lst, lines[0:3])
        energy_output_lst = np.array([float(line.split('  ')[1]) for line in lines[3:]])
        for e_1, e_2 in zip(energy_lst, energy_output_lst):
            self.assertTrue(1 > e_1 > 0)
            self.assertEqual(e_1, e_2)


if __name__ == '__main__':
    unittest.main()
