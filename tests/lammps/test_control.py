import unittest
from pyiron.lammps.control import LammpsControl


class TestLammps(unittest.TestCase):
    def test_generate_seed_from_job(self):
        lc = LammpsControl()
        job_hash_dict = {'job_0_0': lc.generate_seed_from_job(job_name='job_0', seed=0),
                         'job_0_1': lc.generate_seed_from_job(job_name='job_0', seed=1),
                         'job_0_2': lc.generate_seed_from_job(job_name='job_0', seed=2),
                         'job_1_0': lc.generate_seed_from_job(job_name='job_1', seed=0),
                         'job_1_1': lc.generate_seed_from_job(job_name='job_1', seed=1),
                         'job_1_2': lc.generate_seed_from_job(job_name='job_1', seed=2)}
        self.assertEqual(job_hash_dict['job_0_0'], 45625)
        self.assertEqual(job_hash_dict['job_0_1'], 822)
        self.assertEqual(job_hash_dict['job_0_2'], 71750)
        self.assertEqual(job_hash_dict['job_1_0'], 68885)
        self.assertEqual(job_hash_dict['job_1_1'], 56543)
        self.assertEqual(job_hash_dict['job_1_2'], 2668)
