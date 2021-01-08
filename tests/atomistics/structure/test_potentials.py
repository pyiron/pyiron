# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron_atomistic.lammps.potential import LammpsPotentialFile
from pyiron_atomistic.vasp.potential import VaspPotential


class TestLammpsPotentialFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.kim = LammpsPotentialFile()
        cls.potential_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "../../static/lammps/potentials"
        )

    def test_find(self):
        Fe_lst = ["Fe_C_Becquart_eam", "Fe_C_Hepburn_Ackland_eam"]

        self.assertEqual(sorted(list(self.kim.find("Fe")["Name"])), sorted(Fe_lst))
        AlMg_lst = ["Al_Mg_Mendelev_eam"]
        self.assertEqual(sorted(list(self.kim.find({"Al", "Mg"})["Name"])), AlMg_lst)

    def test_pythonic_functions(self):
        self.assertEqual(
            list(self.kim.find("Fe")["Name"]), list(self.kim["Fe"].list()["Name"])
        )
        self.assertEqual(
            list(self.kim.find("Fe")["Name"]), list(self.kim.Fe.list()["Name"])
        )
        self.assertEqual(
            list(self.kim.find({"Al", "Mg"})["Name"]),
            list(self.kim["Al"]["Mg"].list()["Name"]),
        )
        self.assertEqual(
            list(self.kim.find({"Al", "Mg"})["Name"]),
            list(self.kim.Mg.Al.list()["Name"]),
        )


class TestVaspPotential(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vasp = VaspPotential()
        cls.potential_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "../../static/vasp/potentials"
        )

    def test_find(self):
        self.assertEqual(
            list(self.vasp.pbe.find("Fe")["Name"]),
            [
                "Fe-gga-pbe",
                "Fe_GW-gga-pbe",
                "Fe_pv-gga-pbe",
                "Fe_sv-gga-pbe",
                "Fe_sv_GW-gga-pbe",
            ],
        )
        self.assertEqual(
            sorted(list(self.vasp.pbe.find({"Fe", "C"})["Name"])),
            [
                "C-gga-pbe",
                "C_GW-gga-pbe",
                "C_GW_new-gga-pbe",
                "C_h-gga-pbe",
                "C_s-gga-pbe",
                "Fe-gga-pbe",
                "Fe_GW-gga-pbe",
                "Fe_pv-gga-pbe",
                "Fe_sv-gga-pbe",
                "Fe_sv_GW-gga-pbe",
            ],
        )

    def test_pythonic_functions(self):
        self.assertEqual(
            list(self.vasp.pbe.Fe.list()["Name"]),
            [
                "Fe-gga-pbe",
                "Fe_GW-gga-pbe",
                "Fe_pv-gga-pbe",
                "Fe_sv-gga-pbe",
                "Fe_sv_GW-gga-pbe",
            ],
        )
        self.assertEqual(
            sorted(list(self.vasp.pbe.Fe.C.list()["Name"])),
            [
                "C-gga-pbe",
                "C_GW-gga-pbe",
                "C_GW_new-gga-pbe",
                "C_h-gga-pbe",
                "C_s-gga-pbe",
                "Fe-gga-pbe",
                "Fe_GW-gga-pbe",
                "Fe_pv-gga-pbe",
                "Fe_sv-gga-pbe",
                "Fe_sv_GW-gga-pbe",
            ],
        )


if __name__ == "__main__":
    unittest.main()
