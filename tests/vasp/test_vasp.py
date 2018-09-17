import unittest
from pyiron_base.core.settings.generic import Settings
import os

s = Settings(config={'sql_file': 'vasp.db',
                     'project_paths': os.path.abspath(os.getcwd()),
                     'resource_paths': os.path.join(os.path.abspath(os.getcwd()), '../static')})


from pyiron_atomistics.structure.atoms import CrystalStructure
from pyiron_vasp.vasp import Input, Output
from pyiron.project import Project
from pyiron_vasp.potential import VaspPotentialFile

__author__ = "surendralal"


class TestVasp(unittest.TestCase):
    """
    Tests the pyiron.objects.hamilton.dft.vasp.Vasp class
    """

    def setUp(self):
        self.project = Project('test_vasp')
        self.job = self.project.create_job("Vasp", "trial")

    @classmethod
    def tearDownClass(cls):
        s.close_connection()
        os.remove('vasp.db')

    def test_init(self):
        self.assertEqual(self.job.__name__, "Vasp")
        self.assertEqual(self.job._sorted_indices, None)
        self.assertIsInstance(self.job.input, Input)
        self.assertIsInstance(self.job._output_parser, Output)

    def test_calc_static(self):
        self.job.calc_static(electronic_steps=90, retain_charge_density=True, retain_electrostatic_potential=True)
        self.assertEqual(self.job.input.incar["IBRION"], -1)
        self.assertEqual(self.job.input.incar["NELM"], 90)
        self.assertEqual(self.job.input.incar["LVTOT"], True)
        self.assertEqual(self.job.input.incar["LCHARG"], True)

    def test_set_structure(self):
        self.assertEqual(self.job.structure, None)
        atoms = CrystalStructure("Pt", BravaisBasis="fcc", a=3.98)
        self.job.structure = atoms
        self.assertEqual(self.job.structure, atoms)
        self.job.structure = None
        self.assertEqual(self.job.structure, None)
        self.job.structure = atoms
        self.assertEqual(self.job.structure, atoms)

    def test_list_potenitals(self):
        self.assertRaises(ValueError, self.job.list_potentials)

    def tearDown(self):
        pass


class TestInput(unittest.TestCase):

    def setUp(self):
        pass


if __name__ == '__main__':
    unittest.main()