"""
In contrast to many other tests, the non modal tests require a central database and can not be executed with the testing
configuration. Therefore these tests will use your default configuration.
"""
import unittest
import numpy as np
from pyiron.project import Project


class TestExampleJob(unittest.TestCase):
    def setUp(self):
        self.count = 12
        self.project = Project('random_testing_non_modal')
        self.project.remove_jobs(recursive=True)
        self.project.set_logging_level('INFO')

    @classmethod
    def tearDownClass(cls):
        # print('tear down')
        project = Project('random_testing_non_modal')
        project.remove_jobs()
        project.remove(enable=True)

    def test_non_modal_run(self):
        ham_non_modal = self.project.create_job("ExampleJob", "job_non_modal")
        ham_non_modal.input['count'] = self.count
        ham_non_modal.server.run_mode.non_modal = True
        ham_non_modal.run()
        self.assertFalse(ham_non_modal.status.finished)
        self.project.wait_for_job(ham_non_modal, interval_in_s=5, max_iterations=3)
        self.assertTrue(ham_non_modal.status.finished)

        lines = self.project['job_non_modal/input.inp']
        input_lst = ["alat 3.2 #lattice constant (would be in a more realistic example in the structure file)\n",
                     "alpha 0.1 #noise amplitude\n",
                     "a_0 3 #equilibrium lattice constant\n",
                     "a_1 0\n",
                     "a_2 1. #2nd order in energy (corresponds to bulk modulus)\n",
                     "a_3 0. #3rd order\n",
                     "a_4 0. #4th order\n",
                     "count " + str(self.count) + " #number of calls (dummy)\n",
                     "write_restart True\n",
                     "read_restart False\n"]
        self.assertEqual(input_lst, lines)

        lines = self.project['job_non_modal/restart.out']
        restart_lst = ["count " + str(self.count) + " \n"]
        self.assertEqual(restart_lst, lines)
        energy_lst = ham_non_modal.get("output/generic/energy_tot")
        self.assertEqual(len(energy_lst), self.count)

        lines = self.project['job_non_modal/output.log']
        output_lst = ["exampleExecutable logFile \n",
                      "alat 3.2 \n",
                      "count " + str(self.count) + " \n"]
        self.assertEqual(output_lst, lines[0:3])
        energy_output_lst = np.array([float(line.split('  ')[1]) for line in lines[3:]])
        for e_1, e_2 in zip(energy_lst, energy_output_lst):
            self.assertTrue(1 > e_1 > 0)
            self.assertEqual(e_1, e_2)

        ham_non_modal.remove()


if __name__ == '__main__':
    unittest.main()
