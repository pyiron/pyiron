# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from ase.atom import Atom as ASEAtom
import unittest
from pyiron.atomistics.structure.atom import Atom, ase_to_pyiron
from pyiron.atomistics.structure.periodic_table import PeriodicTable


class TestAtom(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.Fe_atom = Atom("Fe")
        pse = PeriodicTable()
        al = pse.element("Al", spin=-1)
        cls.Al_atom = Atom(al)
        cls.Ne_atom = Atom("Ne", charge=-0.1, momentum=0.5, spin=-1)

    def test__init__(self):
        self.Fe_atom.rel = True
        self.assertEqual(self.Fe_atom.element.Abbreviation, "Fe")
        self.assertEqual(self.Fe_atom.position.tolist(), [0, 0, 0])
        self.assertEqual(self.Fe_atom.rel, True)

        self.assertEqual(Atom(Z=13).element.Abbreviation, "Al")
        self.assertRaises(ValueError, Atom)
        self.assertRaises(ValueError, Atom, 13)
        self.assertEqual(Atom("Si", (0, 0, 0)).position.tolist(), [0, 0, 0])

        self.assertEqual(self.Al_atom.symbol, "Al")
        self.assertEqual(self.Al_atom.element.tags["spin"], -1)

    def test_mass(self):
        self.assertEqual(round(Atom("Si", (0, 0, 0)).mass, 4), 28.085)

    def test_number(self):
        self.assertEqual(self.Fe_atom.number, 26)

    def test_z(self):
        self.assertEqual(self.Fe_atom.x, 0)

    def test_ase_compatibility(self):
        self.assertIsInstance(self.Fe_atom, ASEAtom)
        self.assertIsInstance(self.Al_atom, ASEAtom)
        self.assertIsInstance(self.Ne_atom, ASEAtom)
        self.assertEqual(self.Ne_atom.data["spin"], -1)
        self.assertEqual(self.Ne_atom.data["charge"], -0.1)
        self.assertEqual(self.Ne_atom.data["momentum"], 0.5)
        self.assertEqual(self.Ne_atom.spin, -1)
        self.assertEqual(self.Ne_atom.charge, -0.1)
        self.assertEqual(self.Ne_atom.momentum, 0.5)
        self.assertEqual(self.Al_atom.mass, 26.9815385)
        self.assertEqual(self.Fe_atom.mass, 55.845)
        self.assertEqual(self.Ne_atom.mass, 20.1797)
        self.assertEqual(self.Al_atom.symbol, "Al")
        self.assertEqual(self.Fe_atom.symbol, "Fe")
        self.assertEqual(self.Ne_atom.symbol, "Ne")
        self.assertEqual(self.Al_atom.number, 13)
        self.assertEqual(self.Fe_atom.number, 26)
        self.assertEqual(self.Ne_atom.number, 10)
        self.assertNotEqual(self.Ne_atom, self.Al_atom)
        self.assertEqual(self.Fe_atom, Atom("Fe"))
        self.assertNotEqual(self.Fe_atom, ASEAtom("Fe"))
        self.assertEqual(self.Fe_atom, ase_to_pyiron(ASEAtom("Fe")))
        self.assertEqual(Atom("Ne", charge=-0.1, momentum=0.5),
                         ase_to_pyiron(ASEAtom("Ne", charge=-0.1, momentum=0.5)))


if __name__ == "__main__":
    unittest.main()
