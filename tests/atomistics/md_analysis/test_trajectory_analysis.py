import unittest

import numpy as np

from pyiron.atomistics.md_analysis.trajectory_analysis import unwrap_coordinates

__author__ = "surendralal"


class TestTrajectoryAnalysis(unittest.TestCase):

    """
    Testing the TrajectoryAnalysis() module.
    """

    @classmethod
    def setUpClass(cls):
        cls.positions = np.array([[[0.2, 0.2, 0.99], [0., 0., 0.]], [[0.2, 0.2, 0.03], [0., 0., 0.]],
                                  [[0.2, 0.2, 0.02], [0., 0., 0.]], [[0.2, 0.2, 0.98], [0., 0., 0.]],
                                  [[0.2, 0.99, 0.99], [0., 0., 0.]], [[0.2, 0.01, 0.02], [0., 0., 0.]]])

    def test_unwrap_coordinates(self):
        unwrapped_positions = unwrap_coordinates(positions=self.positions, is_relative=True)
        self.assertTrue(np.array_equal(np.shape(self.positions), np.shape(unwrapped_positions)))
        expected_positions = np.array([[[0.2, 0.2, -0.01], [0., 0., 0.]], [[0.2, 0.2, 0.03], [0., 0., 0.]],
                  [[0.2, 0.2, 0.02], [0., 0., 0.]], [[0.2, 0.2, -0.02], [0., 0., 0.]],
                  [[0.2, -0.01, -0.01], [0., 0., 0.]], [[0.2, 0.01, 0.02], [0., 0., 0.]]])
        # print(unwrapped_positions[:, 0], expected_positions[:, 0])
        self.assertTrue(np.allclose(unwrapped_positions, expected_positions))

