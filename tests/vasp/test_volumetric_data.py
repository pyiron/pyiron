# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import posixpath
import numpy as np
from pyiron_atomistic.vasp.volumetric_data import VaspVolumetricData


class TestVaspVolumetricData(unittest.TestCase):

    """
    Testing routines in the vasp/structure module.
    """

    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        vol_directory = os.path.join(
            cls.file_location, "../static/vasp_test_files/chgcar_samples"
        )
        file_list = os.listdir(vol_directory)
        cls.file_list = [posixpath.join(vol_directory, f) for f in file_list]
        cls.vd_obj = VaspVolumetricData()

    def test_read_vol_data(self):
        for chgcar_file in self.file_list:
            if chgcar_file.split("/")[-1] == "CHGCAR_spin":
                atoms, [total_data, diff_data] = self.vd_obj._read_vol_data(
                    chgcar_file, normalize=True
                )
                self.assertEqual(total_data.shape, diff_data.shape)
                self.assertEqual(total_data.shape, (28, 28, 28))
                self.assertEqual(
                    round(np.average(total_data) * atoms.get_volume(), 1), 8
                )
                self.vd_obj.from_file(chgcar_file, normalize=True)
                self.assertTrue(np.array_equal(total_data, self.vd_obj.total_data))
                self.assertTrue(np.array_equal(diff_data, self.vd_obj.diff_data))
                atoms, [total_data_old, diff_data_old] = self.vd_obj._read_vol_data_old(
                    chgcar_file, normalize=True
                )
                self.assertTrue(np.array_equal(total_data_old, self.vd_obj.total_data))
                self.assertTrue(np.array_equal(diff_data_old, self.vd_obj.diff_data))

            if chgcar_file.split("/")[-1] == "CHGCAR_no_spin":
                atoms, [total_data] = self.vd_obj._read_vol_data(
                    chgcar_file, normalize=True
                )
                self.assertEqual(total_data.shape, (28, 28, 28))
                self.assertEqual(
                    round(np.average(total_data) * atoms.get_volume(), 1), 8
                )
                _, [total_data] = self.vd_obj._read_vol_data(
                    chgcar_file, normalize=False
                )
                self.assertEqual(round(np.average(total_data), 4), 8)
                self.vd_obj.from_file(chgcar_file, normalize=False)
                self.assertTrue(np.array_equal(total_data, self.vd_obj.total_data))
                atoms, [total_data_old] = self.vd_obj._read_vol_data_old(
                    chgcar_file, normalize=False
                )
                self.assertTrue(np.array_equal(total_data_old, self.vd_obj.total_data))

            if chgcar_file.split("/")[-1] == "CHGCAR_empty":
                atoms, total_data = self.vd_obj._read_vol_data(
                    chgcar_file, normalize=True
                )
                self.assertIsNone(atoms)
                self.assertIsNone(total_data)
