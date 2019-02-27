import unittest
import os
from pyiron.atomistics.structure.atoms import CrystalStructure
from pyiron.vasp.base import Input, Output
from pyiron.base.project.generic import Project
from pyiron.vasp.potential import VaspPotentialSetter

__author__ = "surendralal"


class TestVasp(unittest.TestCase):
    """
    Tests the pyiron.objects.hamilton.dft.vasp.Vasp class
    """

    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'test_vasp'))
        cls.job = cls.project.create_job("Vasp", "trial")

    @classmethod
    def tearDownClass(cls):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        self.assertEqual(self.job.__name__, "Vasp")
        self.assertEqual(self.job._sorted_indices, None)
        self.assertIsInstance(self.job.input, Input)
        self.assertIsInstance(self.job._output_parser, Output)
        self.assertIsInstance(self.job._potential, VaspPotentialSetter)
        self.assertTrue(self.job._compress_by_default)

    def test_potential(self):
        self.assertEqual(self.job.potential, self.job._potential)

    def test_plane_wave_cutoff(self):
        self.assertIsInstance(self.job.plane_wave_cutoff, (float, int))
        self.job.plane_wave_cutoff = 350
        self.assertEqual(self.job.input.incar["ENCUT"], 350)
        self.assertEqual(self.job.plane_wave_cutoff, 350)
        self.assertEqual(self.job.plane_wave_cutoff, self.job.encut)
        self.job.encut = 450
        self.assertEqual(self.job.encut, 450)
        self.assertEqual(self.job.input.incar["ENCUT"], 450)
        self.assertEqual(self.job.plane_wave_cutoff, 450)

    def test_exchange_correlation_functional(self):
        self.assertEqual(self.job.exchange_correlation_functional, "GGA")
        self.assertEqual(self.job.input.potcar["xc"], "GGA")
        self.job.exchange_correlation_functional = "LDA"
        self.assertEqual(self.job.exchange_correlation_functional, "LDA")
        self.assertEqual(self.job.input.potcar["xc"], "LDA")

    def test_get_nelect(self):
        atoms = CrystalStructure("Pt", BravaisBasis="fcc", a=3.98)
        self.job.structure = atoms
        self.assertEqual(self.job.get_nelect(), 10)

    def test_set_empty_states(self):
        atoms = CrystalStructure("Pt", BravaisBasis="fcc", a=3.98)
        self.job.structure = atoms
        self.job.set_empty_states(n_empty_states=10)
        self.assertEqual(self.job.input.incar["NBANDS"], 15)
        self.job.structure = atoms.repeat([3, 1, 1])
        self.job.set_empty_states(n_empty_states=10)
        self.assertEqual(self.job.input.incar["NBANDS"], 25)

    def test_calc_static(self):
        self.job.calc_static(electronic_steps=90, retain_charge_density=True, retain_electrostatic_potential=True)
        self.assertEqual(self.job.input.incar["IBRION"], -1)
        self.assertEqual(self.job.input.incar["NELM"], 90)
        self.assertEqual(self.job.input.incar["LVTOT"], True)
        self.assertEqual(self.job.input.incar["LCHARG"], True)

    def test_set_structure(self):
        atoms = CrystalStructure("Pt", BravaisBasis="fcc", a=3.98)
        self.job.structure = atoms
        self.assertEqual(self.job.structure, atoms)
        self.job.structure = None
        self.assertEqual(self.job.structure, None)
        self.job.structure = atoms
        self.assertEqual(self.job.structure, atoms)

    def test_list_potenitals(self):
        self.assertRaises(ValueError, self.job.list_potentials)


class TestInput(unittest.TestCase):

    def setUp(self):
        pass


if __name__ == '__main__':
    unittest.main()