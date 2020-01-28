# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron.atomistics.structure.atoms import ase_to_pyiron
from pyiron.atomistics.structure.pyscal import get_steinhardt_parameter_structure
from ase.build import bulk


class TestAtoms(unittest.TestCase):
    def test_bcc_cubic(self):
        if os.name != "nt":
            sysq, ind = get_steinhardt_parameter_structure(
                ase_to_pyiron(
                    bulk(
                        'Fe',
                        cubic=True
                    ).repeat([3, 3, 3])
                )
            )
            self.assertEqual(sum(sysq[0]), 1.9639610121239321)
            self.assertEqual(sum(ind), 53)

    def test_bcc_triclinic(self):
        if os.name != "nt":
            sysq, ind = get_steinhardt_parameter_structure(
                ase_to_pyiron(
                    bulk(
                        'Fe',
                        cubic=False
                    ).repeat([3, 3, 3])
                )
            )
            self.assertEqual(sum(sysq[0]), 2.1642207792631947)
            self.assertTrue(sum(ind) in [9, 18])

    def test_fcc_cubic(self):
        if os.name != "nt":
            sysq, ind = get_steinhardt_parameter_structure(
                ase_to_pyiron(
                bulk(
                    'Al',
                    cubic=True
                ).repeat([3, 3, 3])
                )
            )
            self.assertEqual(sum(sysq[0]), 20.621590627301256)
            self.assertTrue(sum(ind) in [36, 72])

    def test_fcc_triclinic(self):
        if os.name != "nt":
            sysq, ind = get_steinhardt_parameter_structure(
                ase_to_pyiron(
                    bulk(
                        'Al',
                        cubic=False
                    ).repeat([3, 3, 3])
                )
            )
            self.assertEqual(sum(sysq[0]), 4.393254735603798)
            self.assertTrue(sum(ind) in [12, 15])

    def test_hcp_triclinic(self):
        if os.name != "nt":
            sysq, ind = get_steinhardt_parameter_structure(
                ase_to_pyiron(
                    bulk(
                        'Mg',
                        cubic=False
                    ).repeat([3, 3, 3])
                )
            )
            self.assertEqual(sum(sysq[0]), 10.116252659874087)
            self.assertTrue(sum(ind) in [18, 36])
