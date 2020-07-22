import unittest
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
from pyiron.project import Project
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.atomistics.structure.atoms import Atoms, CrystalStructure
import warnings
import numpy as np
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
        project.remove_jobs(recursive=True)
        project.remove(enable=True)

    def test_attributes(self):
        self.assertIsInstance(self.job.structure, Atoms)

    def test_simple_system(self):
        """
        Test a simple ase to pyscal conversion
        """
        sysp = pc.System()
        sysp.read_inputfile(self.job.structure, format="ase")
        assert len(sysp.atoms) == 256

    def test_steinhardt_parameters(self):
        """
        Test the calculation of Steinhardts parameters
        """
        perfect_vals = [0.00, 0.00, 0.190, 0.00, 0.575, 0.00, 0.404, 0.00,
                        0.013, 0.00, 0.600]

        qtest = np.random.randint(2, 13, size=2)

        qs, ind = pas.get_steinhardt_parameter_job(self.job, cutoff=0, n_clusters=2, q=qtest)
        for c, q in enumerate(qs):
            assert np.abs(np.mean(q) - perfect_vals[qtest[c]-2]) < 1E-3    


    def test_centrosymmetry(self):
        csm = pas.analyse_centro_symmetry(self.job.structure, num_neighbors=12)
        assert np.mean(csm) < 1E-5


if __name__ == "__main__":
    unittest.main()

