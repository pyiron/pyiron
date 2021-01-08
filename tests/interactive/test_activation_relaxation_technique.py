import os
import numpy as np
import unittest
from pyiron_atomistic.project import Project
from pyiron_atomistic.atomistics.structure.atoms import Atoms
from pyiron_atomistic.interactive.activation_relaxation_technique import ARTInteractive


class TestARTInteractive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.art = ARTInteractive(art_id = 0, direction = [1, 0, 0])
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'art'))

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'art'))
        cls.project.remove_jobs_silently(recursive=True)
        cls.project.remove(enable=True)

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
        self.art.fix_layer = True
        f_out = self.art.get_forces(f_in)
        self.assertAlmostEqual(np.sum(f_out[1:,0]), 0)
        self.art.fix_layer = False

    def test_run(self):
        basis = Atoms(elements=8*['Fe'],
                      scaled_positions=np.random.random(24).reshape(-1, 3),
                      cell=2.6 * np.eye(3))
        job = self.project.create_job(self.project.job_type.AtomisticExampleJob, "job_single")
        job.server.run_mode.interactive = True
        job.structure = basis
        artint = self.project.create_job('ART', 'job_art')
        artint.ref_job = job
        with self.assertRaises(AssertionError):
            artint.validate_ready_to_run()
        artint.input.art_id = 0
        artint.input.direction = np.ones(3)
        artint.run()
        self.assertEqual(artint.output.forces.shape, (1,8,3))
        artint.interactive_close()
        self.assertTrue(artint.status.finished)

    def test_errors(self):
        with self.assertRaises(ValueError):
            ARTInteractive(art_id = -1, direction = [1, 0, 0])
        with self.assertRaises(ValueError):
            ARTInteractive(art_id = 0, direction = [0, 0, 0])
        with self.assertRaises(ValueError):
            ARTInteractive(art_id = 0, direction = [1, 0, 0], gamma=-0.1)
        with self.assertRaises(ValueError):
            ARTInteractive(art_id = 0, direction = [1, 0, 0], fix_layer=True, non_art_id=[0])


if __name__ == '__main__':
    unittest.main()
