import unittest

import numpy as np

from pyiron.base.objects.volumetric.generic import VolumetricData

"""
@author: surendralal

Unittests for the pyiron.objects.volumetric module
"""


class TestVolumetricData(unittest.TestCase):

    def setUp(self):
        self.data_dict = dict()
        self.data_dict["cubic"] = np.ones((100, 100, 100))
        self.data_dict["non_cubic"] = np.zeros((200, 50, 100))

    def test_total_data_assertion(self):
        vd = VolumetricData()
        with self.assertRaises(ValueError):
            vd.total_data = list()
        with self.assertRaises(TypeError):
            vd.total_data = 1947.0
        for val in self.data_dict.values():
            vd.total_data = val
            self.assertIsInstance(vd.total_data, np.ndarray)

    def test_get_average_along_axis(self):
        vd = VolumetricData()
        for key, val in self.data_dict.items():
            vd.total_data = val
            for i in range(3):
                nz = vd.total_data.shape[i]
                if key == "cubic":
                    answer = np.ones(nz)
                    self.assertTrue(all(np.equal(answer, vd.get_average_along_axis(ind=i))))
                if key == "non-cubic":
                    answer = np.zeros(nz)
                    self.assertTrue(all(np.equal(answer, vd.get_average_along_axis(ind=i))))

if __name__ == '__main__':
    unittest.main()
