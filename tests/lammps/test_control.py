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
        self.assertEqual(job_hash_dict['job_0_0'], 94639)
        self.assertEqual(job_hash_dict['job_0_1'], 84051)
        self.assertEqual(job_hash_dict['job_0_2'], 50062)
        self.assertEqual(job_hash_dict['job_1_0'], 84649)
        self.assertEqual(job_hash_dict['job_1_1'], 99268)
        self.assertEqual(job_hash_dict['job_1_2'], 45752)


if __name__ == '__main__':
    unittest.main()
