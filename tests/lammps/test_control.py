# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from scipy.spatial.transform import Rotation
import unittest
from pyiron_atomistic.lammps.control import LammpsControl, LAMMPS_UNIT_CONVERSIONS


class TestLammps(unittest.TestCase):
    def setUp(self):
        self.lc = LammpsControl()

    def test_generate_seed_from_job(self):
        job_hash_dict = {
            "job_0_0": self.lc.generate_seed_from_job(job_name="job_0", seed=0),
            "job_0_1": self.lc.generate_seed_from_job(job_name="job_0", seed=1),
            "job_0_2": self.lc.generate_seed_from_job(job_name="job_0", seed=2),
            "job_1_0": self.lc.generate_seed_from_job(job_name="job_1", seed=0),
            "job_1_1": self.lc.generate_seed_from_job(job_name="job_1", seed=1),
            "job_1_2": self.lc.generate_seed_from_job(job_name="job_1", seed=2),
        }
        self.assertEqual(job_hash_dict["job_0_0"], 94639)
        self.assertEqual(job_hash_dict["job_0_1"], 84051)
        self.assertEqual(job_hash_dict["job_0_2"], 50062)
        self.assertEqual(job_hash_dict["job_1_0"], 84649)
        self.assertEqual(job_hash_dict["job_1_1"], 99268)
        self.assertEqual(job_hash_dict["job_1_2"], 45752)

    def test_mean(self):
        self.lc.measure_mean_value('energy_pot')
        self.assertEqual(
            self.lc['fix___mean_energy_pot'],
            'all ave/time 1 ${mean_repeat_times} ${thermotime} v_energy_pot'
        )
        self.lc.measure_mean_value('pressures')
        self.assertEqual(self.lc['variable___pressure_1'], 'equal pyy')
        self.lc.measure_mean_value('energy_tot', 2)
        self.assertEqual(
            self.lc['fix___mean_energy_tot'],
            'all ave/time 2 ${mean_repeat_times} ${thermotime} v_energy_tot'
        )
        self.lc.measure_mean_value('volume')
        self.lc.measure_mean_value('temperature')
        self.lc.measure_mean_value('positions')
        self.assertEqual(self.lc['compute___unwrap'], 'all property/atom xu yu zu')
        self.lc.measure_mean_value('forces')
        self.assertEqual(self.lc['variable___forces_0'], 'atom fx')
        self.lc.measure_mean_value('velocities')
        self.assertEqual(self.lc['variable___velocities_0'], 'atom vx')
        with self.assertWarns(Warning):
            self.lc.measure_mean_value('pe**2', name='pepe')
        with self.assertRaises(NotImplementedError):
            self.lc.measure_mean_value('something')

    def test_pressure_to_lammps(self):
        # Correct normalization without rotation. Note that we convert from GPa to bar for LAMMPS.
        no_rot = np.identity(3)
        cnv = LAMMPS_UNIT_CONVERSIONS[self.lc["units"]]["pressure"]
        self.assertTrue(
            np.isclose(self.lc.pressure_to_lammps(0.0, no_rot), 0.0)
        )
        self.assertTrue(
            np.isclose(self.lc.pressure_to_lammps(1.0, no_rot), 1.0*cnv)
        )
        for input_pressure in ([1.0, 2.0, 3.0],
                               [1.0, 2.0, 3.0, None, None, None],
                               [None, None, None, None, None, 2.0],
                               [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                               np.random.uniform(-1, 1, 6),
                               np.random.uniform(-1, 1, 6)):
            output_pressure = [p * cnv if p is not None else None
                               for p in input_pressure]
            out = self.lc.pressure_to_lammps(input_pressure, no_rot)
            for out_i, ref_i in zip(out, output_pressure):
                self.assertTrue(
                    (out_i is None and ref_i is None)
                    or
                    np.isclose(out_i, ref_i)
                )
        # Check if invalid input raises exceptions.
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps("foo", no_rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([], no_rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([1,2,3,4,5,6,7], no_rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([None, None, None, None, None, None], no_rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps(["foo", "bar"], no_rot)
        # With rotation.
        rot = Rotation.random().as_matrix()
        self.assertTrue(
            np.isclose(self.lc.pressure_to_lammps(0.0, rot), 0.0)
        )
        self.assertTrue(
            np.isclose(self.lc.pressure_to_lammps(1.0, rot), 1.0*cnv)
        )
        tmp = self.lc.pressure_to_lammps([1.0, 1.0, 1.0], rot)
        self.assertTrue(
            np.all(np.isclose(tmp[:3], [1.0*cnv, 1.0*cnv, 1.0*cnv]))
            and tmp[3] == tmp[4] == tmp[5] == None
        )
        tmp = self.lc.pressure_to_lammps([1.0, 1.0, 1.0, None, None, None], rot)
        self.assertTrue(
            np.all(np.isclose(tmp[:3], [1.0*cnv, 1.0*cnv, 1.0*cnv]))
            and tmp[3] == tmp[4] == tmp[5] == None
        )
        del tmp
        for input_pressure in ([1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                               [1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
                               [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                               [1.0, -2.0, 3.0, -4.0, -5.0, 6.0],
                               np.random.uniform(-1, 1, 6),
                               np.random.uniform(-1, 1, 6)):
            output_pressure = np.array([[input_pressure[0], input_pressure[3], input_pressure[4]],
                                        [input_pressure[3], input_pressure[1], input_pressure[5]],
                                        [input_pressure[4], input_pressure[5], input_pressure[2]]])
            output_pressure = rot.T @ output_pressure @ rot
            output_pressure = output_pressure[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]] * cnv
            out = self.lc.pressure_to_lammps(input_pressure, rot)
            self.assertTrue(np.all(np.isclose(out, output_pressure)))
        # Check if invalid input raises exceptions.
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([1.0], rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([1.0, None, None], rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([1.0, None, None, None, None, None], rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([1.0, 2.0, 3.0, None, None, None], rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([None, 1.0, 1.0, 0.0, 0.0, 0.0], rot)
        with self.assertRaises(ValueError):
            self.lc.pressure_to_lammps([None, 1.0, 1.0, 1.0, 2.0, 3.0], rot)

    def test_is_isotropic_hydrostatic(self):
        for isotropic_pressure in (
                [0, 0, 0, None, None, None],
                [1., 1., 1., None, None, None],
                [0, 0, 0, 0, 0, 0],
                [1., 1., 1., 0., 0., 0.]

        ):
            self.assertTrue(self.lc._is_isotropic_hydrostatic(isotropic_pressure))
        for nonisotropic_pressure in (
                [0, 0, 0, 1, 1, 1],
                [None, 0, 0, None, None, None],
                [0, 0, 0, 0, 0, None]
        ):
            self.assertFalse(self.lc._is_isotropic_hydrostatic(nonisotropic_pressure))

    def test_fix_move_linear_by_id(self):
        fix = self.lc.fix_move_linear_by_id
        good_ids = np.arange(4, dtype=int)
        good_velocities = [None, -0.004, 0.001]
        fix(good_ids, good_velocities)

        has_str = ['one', 2, 3]
        has_list = [[1, 2, 3], 4, 5]
        empty = []
        self.assertRaises(TypeError, fix, has_str, good_velocities)
        self.assertRaises(TypeError, fix, has_list, good_velocities)
        self.assertRaises(ValueError, fix, empty, good_velocities)
        self.assertRaises(TypeError, fix, [True, True, False], good_velocities)
        self.assertRaises(ValueError, fix, [-2, -1, 0], good_velocities)

        self.assertRaises(TypeError, fix, good_ids, has_str)
        self.assertRaises(TypeError, fix, good_ids, has_list)
        self.assertRaises(ValueError, fix, good_ids, empty)
        self.assertRaises(ValueError, fix, good_ids, np.arange(4))


if __name__ == "__main__":
    unittest.main()
