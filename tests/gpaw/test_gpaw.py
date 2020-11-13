# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import numpy as np
from pyiron_base import Project, ProjectHDFio
from pyiron.gpaw.gpaw import Gpaw
from pyiron.atomistics.structure.atoms import Atoms


class TestGpaw(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, "gpaw"))
        atoms = Atoms("Fe1", positions=np.zeros((1, 3)), cell=np.eye(3))
        job = Gpaw(
            project=ProjectHDFio(project=cls.project, file_name="gpaw"),
            job_name="gpaw",
        )
        job.structure = atoms
        job.encut = 300
        job.set_kpoints([5, 5, 5])
        job.to_hdf()
        cls.job = Gpaw(
            project=ProjectHDFio(project=cls.project, file_name="gpaw"),
            job_name="gpaw",
        )
        cls.job.from_hdf()

    def test_encut(self):
        self.assertEqual(self.job.encut, 300)

    def test_kpoint_mesh(self):
        self.assertEqual(self.job.input["kpoints"], [5, 5, 5])
