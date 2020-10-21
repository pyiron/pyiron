# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import unittest
from pyiron.lammps.control import LammpsControl


class TestLammps(unittest.TestCase):
    def test_generate_seed_from_job(self):
        lc = LammpsControl()
        job_hash_dict = {
            "job_0_0": lc.generate_seed_from_job(job_name="job_0", seed=0),
            "job_0_1": lc.generate_seed_from_job(job_name="job_0", seed=1),
            "job_0_2": lc.generate_seed_from_job(job_name="job_0", seed=2),
            "job_1_0": lc.generate_seed_from_job(job_name="job_1", seed=0),
            "job_1_1": lc.generate_seed_from_job(job_name="job_1", seed=1),
            "job_1_2": lc.generate_seed_from_job(job_name="job_1", seed=2),
        }
        self.assertEqual(job_hash_dict["job_0_0"], 94639)
        self.assertEqual(job_hash_dict["job_0_1"], 84051)
        self.assertEqual(job_hash_dict["job_0_2"], 50062)
        self.assertEqual(job_hash_dict["job_1_0"], 84649)
        self.assertEqual(job_hash_dict["job_1_1"], 99268)
        self.assertEqual(job_hash_dict["job_1_2"], 45752)

    def test_mean(self):
        lc = LammpsControl()
        lc.measure_mean_value('energy_pot')
        self.assertEqual(lc['fix___mean_energy_pot'], 'all ave/time 1 ${mean_repeat_times} ${thermotime} v_energy_pot')
        lc.measure_mean_value('pressures')
        self.assertEqual(lc['variable___pressure_1'], 'equal pyy')
        lc.measure_mean_value('energy_tot', 2)
        self.assertEqual(lc['fix___mean_energy_tot'], 'all ave/time 2 ${mean_repeat_times} ${thermotime} v_energy_tot')
        lc.measure_mean_value('volume')
        lc.measure_mean_value('temperature')
        lc.measure_mean_value('positions')
        self.assertEqual(lc['compute___unwrap'], 'all property/atom xu yu zu')
        lc.measure_mean_value('forces')
        self.assertEqual(lc['variable___forces_0'], 'atom fx')
        lc.measure_mean_value('velocities')
        self.assertEqual(lc['variable___velocities_0'], 'atom vx')
        with self.assertWarns(Warning):
            lc.measure_mean_value('pe**2', name='pepe')
        with self.assertRaises(NotImplementedError):
            lc.measure_mean_value('something')

    def test_pressure_to_lammps(self):
        lc = LammpsControl()
        # Correct normalization without rotation.
        no_rot = np.identity(3)
        self.assertEqual(lc.pressure_to_lammps(1.0, no_rot), [1.0, 1.0, 1.0, None, None, None])
        self.assertEqual(lc.pressure_to_lammps([1.0, 2.0, 3.0], no_rot), [1.0, 2.0, 3.0, None, None, None])
        self.assertEqual(lc.pressure_to_lammps([1.0, 2.0, 3.0, None, None, None], no_rot), [1.0, 2.0, 3.0, None, None, None])
        self.assertEqual(lc.pressure_to_lammps([None, None, None, None, None, 2.0], no_rot), [None, None, None, None, None, 2.0])        

if __name__ == "__main__":
    unittest.main()
