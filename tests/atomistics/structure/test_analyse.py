# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
from pyiron.atomistics.structure.atoms import Atoms, CrystalStructure


class TestAtoms(unittest.TestCase):
    def test_allow_ragged(self):
        a_0 = 4
        struct = CrystalStructure('Al', lattice_constants=a_0, bravais_basis='fcc').repeat(10)
        layers = struct.analyse.get_layers().tolist()
        self.assertEqual(
            layers, np.rint(2*struct.positions/a_0).astype(int).tolist()
        )
        struct.append(Atoms(elements=['C'], positions=np.random.random((1,3))))
        self.assertEqual(
            layers, struct.analyse.get_layers(id_list=struct.select_index('Al')).tolist()
        )
        with self.assertRaises(ValueError):
            _ = struct.analyse.get_layers(id_list=[])


if __name__ == "__main__":
    unittest.main()
