import unittest
import numpy as np
import os
from pyiron.base.project import Project
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.objects.generic.hdfio import ProjectHDFio
from pyiron.lammps.lammps import Lammps


class TestLammps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, 'lammps'))
        cls.job = Lammps(project=ProjectHDFio(project=cls.project, file_name='lammps'), job_name='lammps')

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, 'lammps'))
        project.remove_jobs(recursive=True)
        project.remove(enable=True)

    def test_selective_dynamics(self):
        atoms = Atoms('Fe8', positions=np.zeros((8, 3)), cell=np.eye(3))
        atoms.add_tag(selective_dynamics=[True, True, True])
        self.job.structure = atoms
        self.job._set_selective_dynamics()
        self.assertFalse('group' in self.job.input.control._dataset["Parameter"])
        atoms.add_tag(selective_dynamics=None)
        atoms.selective_dynamics[1] = [True, True, False]
        atoms.selective_dynamics[2] = [True, False, True]
        atoms.selective_dynamics[3] = [False, True, True]
        atoms.selective_dynamics[4] = [False, True, False]
        atoms.selective_dynamics[5] = [False, False, True]
        atoms.selective_dynamics[6] = [True, False, False]
        atoms.selective_dynamics[7] = [False, False, False]
        self.job.structure = atoms
        self.job._set_selective_dynamics()
        self.assertTrue('group' in self.job.input.control._dataset["Parameter"])
        para_lst = np.array(self.job.input.control._dataset["Parameter"])
        self.assertEqual(len(para_lst[para_lst == 'group']), 7)


if __name__ == '__main__':
    unittest.main()
