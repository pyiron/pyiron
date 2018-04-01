import unittest
import os
from pyiron_base.core.settings.generic import Settings
s = Settings(config={'file': 'bandstructure.db',
                     'top_level_dirs': os.path.abspath(os.getcwd()),
                     'resource_paths': os.path.join(os.path.abspath(os.getcwd()), '../static')})

import numpy as np
from pyiron_atomistics.structure.atoms import CrystalStructure
from pyiron_dft.bandstructure import Bandstructure


class TestBandStructure(unittest.TestCase):
    def setUp(self):
        self.bcc = CrystalStructure("Fe", bravais_basis='bcc', lattice_constants=[2.7])
        self.fcc = CrystalStructure("Fe", bravais_basis='fcc', lattice_constants=[3.6])

    def test_bcc(self):
        b1, b2, b3 = np.transpose(np.linalg.inv(self.bcc.cell))
        bandstr = Bandstructure()
        bandstr.structure = self.bcc
        special_points = {"Gamma": [0, 0, 0],
                          "G'": b2,
                          "L": 0.5 * (b1 + b2 + b3),
                          "K": 1. / 8. * (3 * b1 + 6 * b2 + 3 * b3),
                          "U": 1. / 8. * (2 * b1 + 5 * b2 + 5 * b3),
                          "X": 0.5 * (b2 + b3),
                          "W": 0.25 * b1 + 0.75 * b2 + 0.5 * b3,
                          "X'": 0.5 * (b1 + 2 * b2 + b3),
                          "M": 0.5 * (b1 + b2),
                          "X1": 0.5 * b1,
                          "X2": 0.5 * b2
                          }
        path_dict = {'very_short': ['Gamma', 'X'],
                     'full': ['Gamma', 'X', 'L', 'W', 'Gamma']}
        self.assertEqual(path_dict['very_short'], bandstr.path_dict['very_short'])
        self.assertEqual(path_dict['full'], bandstr.path_dict['full'])
        for point in ['K', 'X1', 'X', 'L', "G'", "X'", 'X2', 'U', 'M', 'W', 'Gamma']:
            for coord in range(3):
                self.assertEqual(special_points[point][coord], bandstr.special_points[point][coord])

    def test_fcc(self):
        b1, b2, b3 = np.transpose(np.linalg.inv(self.fcc.cell))
        bandstr = Bandstructure()
        bandstr.structure = self.fcc
        special_points = {"Gamma": [0, 0, 0],
                          "G'": b2,
                          "L": 0.5 * (b1 + b2 + b3),
                          "K": 1. / 8. * (3 * b1 + 6 * b2 + 3 * b3),
                          "U": 1. / 8. * (2 * b1 + 5 * b2 + 5 * b3),
                          "X": 0.5 * (b2 + b3),
                          "W": 0.25 * b1 + 0.75 * b2 + 0.5 * b3,
                          "X'": 0.5 * (b1 + 2 * b2 + b3),
                          "M": 0.5 * (b1 + b2),
                          "X1": 0.5 * b1,
                          "X2": 0.5 * b2
                          }
        path_dict = {'very_short': ['L', 'Gamma', 'X'],
                     'full_CM': ["G'", "X'", 'K', 'Gamma', 'L'],
                     'full_20': ['G', 'X', 'U', 'Gamma', 'L'],
                     'full': ['Gamma', 'X', 'U', 'L', 'Gamma', 'K']}
        self.assertEqual(path_dict['very_short'], bandstr.path_dict['very_short'])
        self.assertEqual(path_dict['full_CM'], bandstr.path_dict['full_CM'])
        self.assertEqual(path_dict['full_20'], bandstr.path_dict['full_20'])
        self.assertEqual(path_dict['full'], bandstr.path_dict['full'])
        for point in ['K', 'X1', 'X', 'L', "G'", "X'", 'X2', 'U', 'M', 'W', 'Gamma']:
            for coord in range(3):
                self.assertEqual(special_points[point][coord], bandstr.special_points[point][coord])


if __name__ == '__main__':
    unittest.main()