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


if __name__ == "__main__":
    unittest.main()

