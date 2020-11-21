# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
from pyiron.atomistics.structure.atoms import CrystalStructure


class TestAtoms(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        pass

    @classmethod
    def setUpClass(cls):
        pass

    def test_allow_ragged(self):
        a_0 = 4
        struct = CrystalStructure(elements='Al', lattice_constants=a_0, bravais_basis='fcc').repeat(10)
        self.assertEqual(
            struct.analyse.get_layers().tolist(), np.rint(2*struct.positions/a_0).astype(int).tolist()
        )


if __name__ == "__main__":
    unittest.main()
