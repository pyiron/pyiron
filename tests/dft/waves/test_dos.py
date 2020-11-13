# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import posixpath
import numpy as np

from pyiron.vasp.vasprun import Vasprun
from pyiron.dft.waves.dos import Dos, NoResolvedDosError

"""
@author: surendralal

Unittests for the pyiron.objects.electronic module
"""


class TestDos(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.es_list = list()
        file_list = ["vasprun_1.xml", "vasprun_2.xml"]
        for f in file_list:
            vp = Vasprun()
            direc = os.path.abspath(
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "../../static/vasp_test_files/vasprun_samples",
                )
            )
            filename = posixpath.join(direc, f)
            vp.from_file(filename)
            es = vp.get_electronic_structure()
            cls.es_list.append(es)

    def test_init_from_es(self):
        for es in self.es_list:
            if es.grand_dos_matrix is not None:
                self.assertIsInstance(es.grand_dos_matrix, np.ndarray)
            dos = Dos(es_obj=es, n_bins=100)
            self.assertIsInstance(dos.energies[0], np.ndarray)
            self.assertIsInstance(dos.t_dos[0], np.ndarray)
            self.assertIsInstance(dos.orbital_dict, dict)
            self.assertIsInstance(dos.n_bins, int)
            self.assertEqual(len(dos.energies), len(dos.t_dos))
            if es.grand_dos_matrix is None:
                self.assertRaises(
                    NoResolvedDosError, dos.get_spatially_resolved_dos, atom_indices=[0]
                )
                self.assertRaises(
                    NoResolvedDosError,
                    dos.get_orbital_resolved_dos,
                    orbital_indices=[0],
                )
                self.assertRaises(
                    NoResolvedDosError,
                    dos.get_spatial_orbital_resolved_dos,
                    atom_indices=[0],
                    orbital_indices=[0],
                )
            else:
                self.assertIsInstance(dos.es_obj.grand_dos_matrix, np.ndarray)

    def test_get_spatially_resolved_dos(self):
        for es in self.es_list:
            dos = Dos(es_obj=es, n_bins=100)
            if es.grand_dos_matrix is not None:
                self.assertIsInstance(es.grand_dos_matrix, np.ndarray)
                _, _, _, n_atoms, _ = np.shape(dos.es_obj.grand_dos_matrix)
                total_rdos = np.zeros_like(dos.t_dos)
                for i in range(n_atoms):
                    r_dos = dos.get_spatially_resolved_dos([i])
                    total_rdos += r_dos
                self.assertTrue(np.allclose(dos.t_dos, total_rdos))

    def test_get_orbital_resolved_dos(self):
        for es in self.es_list:
            dos = Dos(es_obj=es, n_bins=100)
            if es.grand_dos_matrix is not None:
                self.assertIsInstance(es.grand_dos_matrix, np.ndarray)
                _, _, _, _, n_orbitals = np.shape(dos.es_obj.grand_dos_matrix)
                total_rdos = np.zeros_like(dos.t_dos)
                for i in range(n_orbitals):
                    r_dos = dos.get_orbital_resolved_dos([i])
                    total_rdos += r_dos
                self.assertTrue(np.allclose(dos.t_dos, total_rdos))

    def test_get_spatial_orbital_resolved_dos(self):
        for es in self.es_list:
            dos = Dos(es_obj=es, n_bins=100)
            if es.grand_dos_matrix is not None:
                self.assertIsInstance(es.grand_dos_matrix, np.ndarray)
                _, _, _, n_atoms, n_orbitals = np.shape(dos.es_obj.grand_dos_matrix)
                atom_indices = np.arange(n_atoms)
                orbital_indices = np.arange(n_orbitals)
                r_dos = dos.get_spatial_orbital_resolved_dos(
                    atom_indices=atom_indices, orbital_indices=orbital_indices
                )
                self.assertTrue(np.allclose(dos.t_dos, r_dos))
