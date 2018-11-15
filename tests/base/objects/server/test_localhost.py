import os
from base.objects.server.scheduler.localhost import Localhost
import unittest


class TestLocalhost(unittest.TestCase):
    def setUp(self):
        dir_path = str(os.path.dirname(os.path.realpath(__file__)))
        self.localhost = Localhost(working_directory=dir_path)

    def test_cores(self):
        self.assertEqual(self.localhost.cores, 1)
        try:
            self.localhost.cores = 10
        except ValueError:
            pass
        self.assertNotEqual(self.localhost.cores, 10)
        self.assertEqual(self.localhost.cores, 1)

    def test_support_wait_for_prev_job(self):
        self.assertFalse(self.localhost.support_wait_for_prev_job)

    def test_support_run_time_limit(self):
        self.assertFalse(self.localhost.support_run_time_limit)

    def test_support_cores_limit(self):
        self.assertFalse(self.localhost.support_cores_limit)

    def test_support_working_directory(self):
        self.assertFalse(self.localhost.support_working_directory)

    def test_active_scheduler(self):
        self.assertEqual(self.localhost.active_scheduler.__name__, 'default')
        self.localhost.active_scheduler = 'default'
        self.assertEqual(self.localhost.active_scheduler.__name__, 'default')


if __name__ == '__main__':
    unittest.main()