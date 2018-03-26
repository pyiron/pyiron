import unittest
from pyiron_base.objects.server.queue import Queue


class TestQue(unittest.TestCase):
    def setUp(self):
        self.que1 = Queue(name='1_1_5_1', mini_cores=1, maxi_cores=5, divisor_list=1, run_time_limit=1)
        self.que2 = Queue(name='2_2_2_1', mini_cores=2, maxi_cores=2, divisor_list=[2], run_time_limit=1)

    def test_minimum_number_of_cores(self):
        self.assertEqual(self.que1.minimum_number_of_cores, 1)
        self.assertEqual(self.que2.minimum_number_of_cores, 2)

    def test_maximum_number_of_cores(self):
        self.assertEqual(self.que1.maximum_number_of_cores, 5)
        self.assertEqual(self.que2.maximum_number_of_cores, 2)

    def test_divisors_for_que(self):
        self.assertEqual(self.que1.divisors_for_que, [1])
        self.assertEqual(self.que2.divisors_for_que, [2])

    def test_cores(self):
        try:
            self.que1.cores = 1
        except ValueError():
            pass
        self.assertEqual(self.que1.cores, 1)

    def test_run_time(self):
        try:
            self.que1.run_time = 1
        except ValueError():
            pass
        self.assertEqual(self.que1.run_time, 1)
        try:
            self.que2.run_time = 2
        except ValueError:
            pass
        self.assertFalse(self.que2.run_time)


if __name__ == '__main__':
    unittest.main()