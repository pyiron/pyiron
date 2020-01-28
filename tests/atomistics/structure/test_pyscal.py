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
            self.assertEqual(sysq, 0)
            self.assertEqual(ind, 0)

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
            self.assertEqual(sysq, 0)
            self.assertEqual(ind, 0)

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
            self.assertEqual(sysq, 0)
            self.assertEqual(ind, 0)

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
            self.assertEqual(sysq, 0)
            self.assertEqual(ind, 0)

    def test_hcp_cubic(self):
        if os.name != "nt":
            sysq, ind = get_steinhardt_parameter_structure(
             ase_to_pyiron(
                 bulk(
                     'Mg',
                     cubic=True
                 ).repeat([3, 3, 3])
             )
            )
            self.assertEqual(sysq, 0)
            self.assertEqual(ind, 0)

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
            self.assertEqual(sysq, 0)
            self.assertEqual(ind, 0)