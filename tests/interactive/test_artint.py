import os
import numpy as np
import unittest
from pyiron.project import Project
from pyiron.atomistics.structure.atoms import Atoms
from pyiron_mpie.interactive.artint import ART


class TestARTInteractive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.art = ART(art_id = 0, direction = [1, 0, 0])
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'art'))
        cls.basis = Atoms(elements=8*['Fe'],
                           scaled_positions=np.random.random(24).reshape(-1, 3),
                           cell=2.6 * np.eye(3))
        job = cls.project.create_job(cls.project.job_type.AtomisticExampleJob, "job_single")
        job.server.run_mode.interactive = True
        job.structure = cls.basis
        cls.artint = cls.project.create_job('ARTInteractive', 'job_art')
        cls.artint.ref_job = job
        cls.artint.input.art_id = 0
        cls.artint.input.direction = np.ones(3)
        cls.artint.run()

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'art'))
        cls.project.remove_jobs_silently(recursive=True)

    def test_R(self):
        self.assertAlmostEqual(np.linalg.norm(self.art._R-np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]])), 0)

    def test_get_forces(self):
        f_in = np.random.random(24).reshape(-1, 3)
        f_in -= np.mean(f_in, axis=0)
        f_out = self.art.get_forces(f_in)
        self.assertAlmostEqual(np.sum(f_out), 0)
        self.assertLessEqual(f_in[0, 0]*f_out[0, 0], 0)
        self.assertGreaterEqual(f_in[0, 1]*f_out[0, 1], 0)
        self.assertGreaterEqual(f_in[0, 2]*f_out[0, 2], 0)

    def test_forces(self):
        self.assertEqual(self.artint.output.forces[-1].shape, (8,3))


if __name__ == '__main__':
    unittest.main()
