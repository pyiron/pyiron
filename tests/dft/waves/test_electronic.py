import unittest
import os
import posixpath
import numpy as np

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.vasp.vasprun import Vasprun
from pyiron.dft.waves.dos import Dos
from pyiron.dft.waves.electronic import ElectronicStructure

"""
@author: surendralal

Unittests for the pyiron.objects.electronic module
"""


class TestElectronicStructure(unittest.TestCase):

    def setUp(self):
        self.es_list = list()
        self.es_obj = ElectronicStructure()
        file_list = ["vasprun_1.xml", "vasprun_2.xml"]
        for f in file_list:
            vp = Vasprun()
            self.assertIsInstance(vp.vasprun_dict, dict)
            direc = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../static/vasprun_samples"))
            filename = posixpath.join(direc, f)
            vp.from_file(filename)
            es = vp.get_electronic_structure()
            self.es_list.append(es)

    def test_init(self):
        for es in self.es_list:
            self.assertIsInstance(es, ElectronicStructure)
            self.assertIsInstance(es.eigenvalues, np.ndarray)
            self.assertIsInstance(es.occupancies, np.ndarray)
            self.assertEqual(np.shape(es.occupancy_matrix), np.shape(es.eigenvalue_matrix))
            if es.structure is not None:
                self.assertIsInstance(es.structure, Atoms)

    def test_add_kpoint(self):
        self.assertEqual(len(self.es_obj.kpoints), 0)
        self.es_obj.add_kpoint(value=[0., 0., 0.], weight=1.0)
        self.assertEqual(len(self.es_obj.kpoints), 1)

    def test_get_dos(self):
        for es in self.es_list:
            dos = es.get_dos()
            self.assertIsInstance(dos, Dos)

    def test_eigenvalues(self):
        for es in self.es_list:
            self.assertEqual(len(es.eigenvalues), np.product(np.shape(es.eigenvalue_matrix)))

    def test_occupancies(self):
        for es in self.es_list:
            self.assertEqual(len(es.occupancies), np.product(np.shape(es.occupancy_matrix)))

