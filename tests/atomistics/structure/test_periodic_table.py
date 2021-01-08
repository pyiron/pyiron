# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron_atomistic.atomistics.structure.atoms import CrystalStructure
from pyiron_atomistic.atomistics.structure.periodic_table import PeriodicTable
from pyiron_base import Project


class TestPeriodicTable(unittest.TestCase):
    """
    define unittest snippets for all python files in the structure directory
    these tests should run fast (a few 10 ms)
    TODO: add write and load to h5 (e.g. addElement needs respective changes in read/load routines)
    """

    @classmethod
    def setUpClass(cls):
        cls.pse = PeriodicTable()

    def test_numbertechnic(self):
        el1 = self.pse.element(1)
        self.assertEqual(el1.Abbreviation, "H")

    def test_Element_by_Abbreviation(self):
        el1 = self.pse.element("Na")
        self.assertEqual(el1.Abbreviation, "Na")

    def test_Element_by_Index(self):
        el1 = self.pse.element(20)
        self.assertEqual(el1.Abbreviation, "Ca")

    def test_Abbreviation_range(self):
        self.assertEqual(len(self.pse.dataframe.Abbreviation[self.pse.Period < 4]), 18)

    def test_add_element_without_tags(self):
        fe_up = self.pse.add_element("Fe", "B_up")
        self.assertEqual(int(fe_up.MeltingPoint), 1808)

    def test_add_Abbreviation_bug(self):
        fe_up = self.pse.add_element("Fe", "B_up")
        self.assertEqual(fe_up.Abbreviation, "B_up")

    def test_add_element_tags(self):
        fe_up = self.pse.add_element(
            "Fe", "Fe_up", spin="up", pseudo_name="GGA", testtag="testtest"
        )
        self.assertEqual(fe_up.Abbreviation, "Fe_up")
        self.assertEqual(fe_up.tags["spin"], "up")
        self.assertEqual(fe_up.tags["pseudo_name"], "GGA")
        self.assertEqual(fe_up.tags["testtag"], "testtest")

    def test_atomic_mass(self):
        el1 = self.pse.element("Fe")
        self.assertAlmostEqual(el1.AtomicMass, 55.845, places=3)

    def test_group(self):
        el1 = self.pse.element("Fe")
        self.assertEqual(el1.Group, 8)

    def test_Period(self):
        el1 = self.pse.element("Fe")
        self.assertEqual(el1.Period, 4)

    def test_add_MeltingPoint(self):
        el1 = self.pse.element("Fe")
        self.assertEqual(int(el1.MeltingPoint), 1808)

    def test_set_item(self):
        el1 = self.pse.element("Fe")
        el1.MeltingPoint = 1900
        self.assertEqual(int(el1.MeltingPoint), 1900)

    def test_is_element(self):
        self.assertEqual(self.pse.is_element("Fe"), True)

    def test_Chemical_Element_to_and_from_hdf(self):
        ni_up = self.pse.add_element("Ni", "Ni_up", spin="up")
        pr = Project(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "test_periodic_table"
            )
        )
        basis = CrystalStructure(
            element=ni_up, bravais_basis="fcc", lattice_constant=3.7
        )
        ham = pr.create_job(pr.job_type.Lammps, "lammps_test_1")
        test_ham = pr.create_job(pr.job_type.Lammps, "lammps_test_1")
        ham.structure = basis
        ham.to_hdf()
        test_ham.from_hdf()
        self.assertEqual(
            test_ham["input/structure/species"][0], ham["input/structure/species"][0]
        )
        ham.remove()

    def test_Chemical_Element_to_and_from_hdf_with_None_Parent(self):
        pr = Project(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "test_periodic_table"
            )
        )
        basis = CrystalStructure(
            element="Ni", bravais_basis="fcc", lattice_constant=3.7
        )
        ham = pr.create_job(pr.job_type.Lammps, "lammps_test_2")
        test_ham = pr.create_job(pr.job_type.Lammps, "lammps_test_2")
        ham.structure = basis
        ham.to_hdf()
        test_ham.from_hdf()
        self.assertEqual(
            test_ham["input/structure/species"][0], ham["input/structure/species"][0]
        )
        ham.remove()

    def test_add_tags(self):
        tag_dic = {"a": "b", "c": "d", "e": "f"}
        fe_up = self.pse.add_element(
            "Fe", "Fe_up", spin="up", pseudo_name="GGA", testtag="testtest"
        )
        fe_up.add_tags(tag_dic)
        self.assertEqual(fe_up.tags["spin"], "up")
        self.assertEqual(fe_up.tags["a"], "b")
        self.assertEqual(fe_up.tags["c"], "d")
        self.assertEqual(fe_up.tags["e"], "f")

    def test__eq__(self):
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        o_1 = pse.element("O")
        o_2 = pse.element("O")
        h_1 = pse.element("H")
        self.assertNotEqual(o_up, o_1)
        self.assertNotEqual(o_up, o_2)
        self.assertEqual(o_1, o_2)
        self.assertNotEqual(o_1, h_1)

    def test_gt__(self):
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_up = pse.element("O_up")
        o_1 = pse.element("O")
        o_2 = pse.element("O")
        h_1 = pse.element("H")
        self.assertTrue(o_up > o_1)
        self.assertFalse(o_up < o_2)
        self.assertFalse(o_1 > o_2)
        self.assertFalse(o_1 < o_2)
        self.assertTrue(o_1 > h_1)
        self.assertTrue(o_up > h_1)

    def test_ge__(self):
        pse = PeriodicTable()
        pse.add_element("O", "O_up", spin="up")
        o_1 = pse.element("O")
        o_2 = pse.element("O")
        self.assertTrue(o_1 <= o_2)
        self.assertTrue(o_1 >= o_2)


if __name__ == "__main__":
    unittest.main()
