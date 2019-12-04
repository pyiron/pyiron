# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron.vasp.potential import get_enmax_among_species


class TestPotential(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))

    def test_get_enmax_among_species(self):
        float_out = get_enmax_among_species(['Fe'], return_list=False)
        self.assertTrue(isinstance(float_out, float))

        tuple_out = get_enmax_among_species(['Fe'], return_list=True)
        self.assertTrue(isinstance(tuple_out, tuple))

        self.assertRaises(KeyError, get_enmax_among_species, species_lst=['X'])
        self.assertRaises(ValueError, get_enmax_among_species, species_lst=['Fe'], xc='FOO')


if __name__ == "__main__":
    unittest.main()
