import unittest
import numpy as np
import pyscal.core as pc
from pyiron.project import Project
from pyiron_base import ProjectHDFio
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.atomistics.structure.atoms import Atoms, CrystalStructure
import os
import pyiron.atomistics.structure.pyscal as pas


class Testpyscal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, "test_job"))
        cls.job = AtomisticGenericJob(
            project=ProjectHDFio(project=cls.project, file_name="test_job"),
            job_name="test_job",
        )
        cls.job.structure = CrystalStructure(
            element="Al", bravais_basis="fcc", lattice_constants=4
        ).repeat(4)

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, "test_job"))
        project.remove_jobs_silently(recursive=True)
        project.remove(enable=True)

    def test_attributes(self):
        self.assertIsInstance(self.job.structure, Atoms)

    def test_simple_system(self):
        """
        Test a simple ase to pyscal conversion
        """
        sysp = pc.System()
        sysp.read_inputfile(self.job.structure, format="ase")
        self.assertEqual(len(sysp.atoms), 256)

    def test_steinhardt_parameters(self):
        """
        Test the calculation of Steinhardts parameters
        """
        perfect_vals = [0.00, 0.00, 0.190, 0.00, 0.575, 0.00, 0.404, 0.00,
                        0.013, 0.00, 0.600]

        qtest = np.random.randint(2, 13, size=2)

        qs, _ = pas.get_steinhardt_parameter_job(self.job, cutoff=0, n_clusters=2, q=qtest)
        for c, q in enumerate(qs):
            self.assertLess(np.abs(np.mean(q) - perfect_vals[qtest[c]-2]), 1E-3)

    def test_centrosymmetry(self):
        csm = pas.analyse_centro_symmetry(self.job.structure, num_neighbors=12)
        self.assertLess(np.mean(csm), 1E-5)

    def test_cna(self):
        cna = pas.analyse_cna_adaptive(self.job.structure)
        self.assertEqual(cna[1], len(self.job.structure))

        rand = np.random.randint(0, len(self.job.structure))

        cna = pas.analyse_cna_adaptive(self.job.structure, mode="numeric")
        self.assertEqual(cna[rand], 1)

        cna = pas.analyse_cna_adaptive(self.job.structure, mode="str")
        self.assertEqual(cna[rand], "fcc")

    def test_volume(self):
        vols = pas.analyse_voronoi_volume(self.job.structure)
        self.assertLess(np.abs(np.mean(vols) - 16.0), 1E-3)

if __name__ == "__main__":
    unittest.main()
