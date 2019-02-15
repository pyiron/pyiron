import unittest
import os
from pyiron.vasp.base import Input, Output
from pyiron.base.project.generic import Project

__author__ = "surendralal"


class TestVasp(unittest.TestCase):
    """
    Tests the pyiron.objects.hamilton.dft.vasp.Vasp class
    """

    def setUp(self):
        self.file_location = os.path.dirname(os.path.abspath(__file__))
        self.project = Project(os.path.join(self.file_location, 'test_vasp'))
        self.job = self.project.create_job("Vasp", "trial")

    @classmethod
    def tearDownClass(cls):
        pass

    def test_init(self):
        self.assertEqual(self.job.__name__, "Vasp")

    def tearDown(self):
        pass


class TestInput(unittest.TestCase):

    def setUp(self):
        pass


if __name__ == '__main__':
    unittest.main()