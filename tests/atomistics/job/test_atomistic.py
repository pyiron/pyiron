# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
import os
from pyiron_atomistic.project import Project
from pyiron_base import ProjectHDFio
from pyiron_atomistic.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistic.atomistics.structure.atoms import Atoms, CrystalStructure
import warnings


class TestAtomisticGenericJob(unittest.TestCase):

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

    def test_get_displacements(self):
        n_steps = 10
        # increasing each position by 0.5 at each step
        positions = np.array([self.job.structure.positions + 0.5 * i for i in range(n_steps)])
        # constant cell
        cells = np.array([self.job.structure.cell] * n_steps)
        disp = self.job.output.get_displacements(self.job.structure, positions=positions, cells=cells)
        disp_ref = np.ones_like(positions) * 0.5
        disp_ref[0] *= 0.0
        self.assertTrue(np.allclose(disp, disp_ref))
        # varying cell
        cells = np.array([self.job.structure.cell * ((i+1) / 10) for i in range(n_steps)])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            disp = self.job.output.get_displacements(self.job.structure,
                                                     positions=positions, cells=cells)
            self.assertEqual(len(w), 1)
            self.assertIsInstance(w[-1].message, UserWarning)
        self.assertFalse(np.allclose(disp, disp_ref))
        dummy_struct = self.job.structure.copy()
        disp_ref = list()
        for pos, cell in zip(positions, cells):
            pos_init = dummy_struct.get_scaled_positions().copy()
            dummy_struct.set_cell(cell, scale_atoms=False)
            dummy_struct.positions = pos
            dummy_struct.center_coordinates_in_unit_cell()
            diff = dummy_struct.get_scaled_positions()-pos_init
            diff[diff >= 0.5] -= 1.0
            diff[diff <= -0.5] += 1.0
            disp_ref.append(np.dot(diff, cell))
        self.assertTrue(np.allclose(disp, disp_ref))


if __name__ == "__main__":
    unittest.main()
