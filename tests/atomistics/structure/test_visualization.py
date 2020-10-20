# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
from pyiron.atomistics.structure.visualization import _get_flattened_orientation


class TestAtoms(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        pass

    @classmethod
    def setUpClass(cls):
        pass

    def test_get_flattened_orientation(self):
        R = np.random.random(9).reshape(-1, 3)
        R = np.array(_get_flattened_orientation(R, 1)).reshape(4, 4)
        self.assertAlmostEqual(np.linalg.det(R), 1)


if __name__ == "__main__":
    unittest.main()
