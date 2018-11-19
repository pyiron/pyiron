import os
from pyiron.base.objects.server.scheduler.cmmc import Cmmc
import unittest


class TestCmmcQueueSystem(unittest.TestCase):
    def setUp(self):
        dir_path = '.'
        self.cmmc_hydra_small = Cmmc(working_directory=dir_path)
        self.cmmc_hydra = Cmmc(working_directory=dir_path, scheduler_name='impi_hydra')
        self.cmmc_hydra_switch = Cmmc(working_directory=dir_path, scheduler_name='impi_hydra')
        self.cmmc_hydra_cores = Cmmc(working_directory=dir_path, scheduler_name='impi_hydra', cores=100)
        self.cmmc_hydra_runtime = Cmmc(working_directory=dir_path, scheduler_name='impi_hydra', runtime=40000)

    def test_active_que(self):
        self.assertEqual(self.cmmc_hydra_switch.active_scheduler.__name__, 'impi_hydra')
        self.cmmc_hydra_switch.active_scheduler = 'impi_hydra_small'
        self.assertEqual(self.cmmc_hydra_switch.active_scheduler.__name__, 'impi_hydra_small')
        self.cmmc_hydra_switch.active_scheduler = 'impi_hydra'

    def test_cores(self):
        self.assertEqual(self.cmmc_hydra.cores, self.cmmc_hydra.active_scheduler.minimum_number_of_cores)
        self.assertEqual(self.cmmc_hydra_cores.cores, 100)

    def test_run_time(self):
        self.assertEqual(self.cmmc_hydra.run_time, 20000)
        self.assertEqual(self.cmmc_hydra_runtime.run_time, 40000)
        try:
            self.cmmc_hydra.run_time = 259201
        except ValueError:
            pass
        self.assertNotEqual(self.cmmc_hydra.run_time, 259201)

    def test_working_directory(self):
        self.assertTrue(isinstance(self.cmmc_hydra_small.working_directory, str))
        self.assertEqual(self.cmmc_hydra_small.working_directory, '.')

    def test_que_script_name(self):
        self.assertTrue(isinstance(self.cmmc_hydra_small.script_name, str))
        self.assertEqual(self.cmmc_hydra_small.script_name, 'run_queue.sh')
        self.cmmc_hydra_small.script_name = 'run_queue_2.sh'
        self.assertTrue(isinstance(self.cmmc_hydra_small.script_name, str))

    def test_que_wrapper(self):
        self.assertTrue(isinstance(self.cmmc_hydra_small.wrapper, list))
        self.assertEqual(self.cmmc_hydra_small.wrapper, ['#!/bin/bash', '', 'python run_job.py'])
        self.cmmc_hydra_small.wrapper = '#!/bin/bash\n\npython run_job.py'
        self.assertTrue(isinstance(self.cmmc_hydra_small.wrapper, list))
        self.assertEqual(self.cmmc_hydra_small.wrapper, ['#!/bin/bash', '', 'python run_job.py'])

    def test_list_que_options(self):
        self.cmmc_hydra_small.script_name = 'run_queue_2.sh'
        self.assertEqual(self.cmmc_hydra_small.list_scheduler_options(),
                         ['qsub', '-terse', '-l', 'h_rt=20000', '-o', 'time.out', '-e', 'error.out',
                          '-wd', '.', '-V', '-S', '/bin/bash',
                          '-pe', 'impi_hydra_small', '1',
                          './run_queue_2.sh'])
        self.assertEqual(self.cmmc_hydra.list_scheduler_options(),
                         ['qsub', '-terse', '-l', 'h_rt=20000', '-o', 'time.out', '-e', 'error.out',
                          '-wd', '.', '-V', '-S', '/bin/bash',
                          '-pe', 'impi_hydra', '20',
                          './run_queue.sh'])
        self.cmmc_hydra_cores.cores = 200
        self.assertEqual(self.cmmc_hydra_cores.list_scheduler_options(),
                         ['qsub', '-terse', '-l', 'h_rt=20000', '-o', 'time.out', '-e', 'error.out',
                          '-wd', '.', '-V', '-S', '/bin/bash',
                          '-pe', 'impi_hydra', '200',
                          './run_queue.sh'])
        self.cmmc_hydra_runtime.run_time = 80000
        self.assertEqual(self.cmmc_hydra_runtime.list_scheduler_options(),
                         ['qsub', '-terse', '-l', 'h_rt=80000', '-o', 'time.out', '-e', 'error.out',
                          '-wd', '.', '-V', '-S', '/bin/bash',
                          '-pe', 'impi_hydra', '20',
                          './run_queue.sh'])

    def test_write_que_wrapper(self):
        dir_path = '.'
        self.cmmc_hydra.write_wrapper()
        with open(dir_path + '/run_queue.sh', 'r') as run_warpper:
            lines_lst = run_warpper.read()
            lines = lines_lst.split('\n')[:-1]
            self.assertEqual(lines, self.cmmc_hydra.wrapper)
        os.remove(dir_path + '/run_queue.sh')

    def test_available_ques_dict(self):
        self.assertEqual(str(self.cmmc_hydra_small.available_schedulers_dict()), str(self.cmmc_hydra.available_schedulers_dict()))


if __name__ == '__main__':
    unittest.main()