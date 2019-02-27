import unittest
from pyiron.atomistics.structure.atom import Atom
from pyiron.atomistics.structure.periodic_table import PeriodicTable


class TestAtom(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.Fe_atom = Atom("Fe")
        pse = PeriodicTable()
        al = pse.element("Al", spin=-1)
        cls.Al_atom = Atom(al)

    def test__init__(self):
        self.Fe_atom.rel = True
        self.assertEqual(self.Fe_atom.element.Abbreviation, "Fe")
        self.assertEqual(self.Fe_atom.position.tolist(), [0, 0, 0])
        self.assertEqual(self.Fe_atom.rel, True)

        self.assertEqual(Atom(Z=13).element.Abbreviation, "Al")
        self.assertRaises(ValueError, Atom)
        self.assertRaises(ValueError, Atom, 13)
        self.assertEqual(Atom('Si', (0, 0, 0)).position.tolist(), [0, 0, 0])

        self.assertEqual(self.Al_atom.symbol, 'Al')
        self.assertEqual(self.Al_atom.element.tags["spin"], -1)

    def test_mass(self):
        self.assertEqual(round(Atom('Si', (0, 0, 0)).mass, 4), 28.0855)

    def test_number(self):
        self.assertEqual(self.Fe_atom.number, 26)

    def test_z(self):
        self.assertEqual(self.Fe_atom.x, 0)


if __name__ == '__main__':
    unittest.main()