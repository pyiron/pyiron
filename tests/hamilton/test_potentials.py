import os
import unittest

from pyiron_base.core.settings.config.testing import ConfigTesting
from pyiron_base.core.settings.generic import Settings

config = ConfigTesting(sql_lite_database='./testing_random.db', path_project=str(os.getcwd()),
                       path_potentials='../../../static/potentials/')
s = Settings(config=config)

from pyiron_lammps.potential import LammpsPotentialFile
from pyiron_vasp.potential import VaspPotential


class TestOpenKimPotential(unittest.TestCase):
    def setUp(self):
        self.kim = LammpsPotentialFile()
        self.potential_path = os.path.abspath(s.path_potentials)

    def test_find(self):
        Fe_lst = ['Al_Fe_eam_fs',
                  'FeCuNi_eam_alloy',
                  'FeNiCr_Bonny_2013_ptDef_eam_alloy',
                  'FeNiCr_eam_alloy',
                  'Fe_2_eam_fs',
                  'Fe_2_eam_fs',
                  'Fe_5_eam_fs',
                  'Fe_C_Hepburn_Ackland_eam_fs',
                  'Fe_Mishin2006_eam_alloy',
                  'Fe_Ni_eam_alloy',
                  'Fe_P_eam_fs',
                  'Fe_eam_fs',
                  'Fe_eam_fs',
                  'V_Fe_eam_fs']

        self.assertEqual(sorted(list(self.kim.find("Fe")['Name'])), sorted(Fe_lst))
        AlMg_lst = ['Al_Mg_eam_fs', 'mg_al_set_eam_alloy']
        self.assertEqual(sorted(list(self.kim.find({"Al", "Mg"})['Name'])), AlMg_lst)

    def test_pythonic_functions(self):
        self.assertEqual(list(self.kim.find("Fe")['Name']), list(self.kim["Fe"].list()['Name']))
        self.assertEqual(list(self.kim.find("Fe")['Name']), list(self.kim.Fe.list()['Name']))
        self.assertEqual(list(self.kim.find({"Al", "Mg"})['Name']), list(self.kim["Al"]["Mg"].list()['Name']))
        self.assertEqual(list(self.kim.find({"Al", "Mg"})['Name']), list(self.kim.Mg.Al.list()['Name']))


class TestVaspPotential(unittest.TestCase):
    def setUp(self):
        self.vasp = VaspPotential()
        self.potential_path = os.path.abspath(s.path_potentials)

    def test_find(self):
        self.assertEqual(list(self.vasp.pbe.find('Fe')['Name']), ['Fe-gga-pbe', 'Fe_GW-gga-pbe', 'Fe_pv-gga-pbe',
                                                                  'Fe_sv-gga-pbe', 'Fe_sv_GW-gga-pbe'])
        self.assertEqual(sorted(list(self.vasp.pbe.find({'Fe', 'C'})['Name'])),
                         ['C-gga-pbe', 'C_GW-gga-pbe', 'C_GW_new-gga-pbe', 'C_h-gga-pbe', 'C_s-gga-pbe', 'Fe-gga-pbe',
                          'Fe_GW-gga-pbe', 'Fe_pv-gga-pbe', 'Fe_sv-gga-pbe', 'Fe_sv_GW-gga-pbe'])

    def test_pythonic_functions(self):
        self.assertEqual(list(self.vasp.pbe.Fe.list()['Name']), ['Fe-gga-pbe', 'Fe_GW-gga-pbe', 'Fe_pv-gga-pbe',
                                                                  'Fe_sv-gga-pbe', 'Fe_sv_GW-gga-pbe'])
        self.assertEqual(sorted(list(self.vasp.pbe.Fe.C.list()['Name'])),
                         ['C-gga-pbe', 'C_GW-gga-pbe', 'C_GW_new-gga-pbe', 'C_h-gga-pbe', 'C_s-gga-pbe', 'Fe-gga-pbe',
                          'Fe_GW-gga-pbe', 'Fe_pv-gga-pbe', 'Fe_sv-gga-pbe', 'Fe_sv_GW-gga-pbe'])


if __name__ == '__main__':
    unittest.main()
