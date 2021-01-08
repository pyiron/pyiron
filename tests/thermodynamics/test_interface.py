# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
import numpy as np
import pandas
from ase.build import bulk
from pyiron_atomistic.atomistics.structure.atoms import ase_to_pyiron
import pyiron_atomistic
from pyiron_base import Project
from pyiron_atomistic.atomistics.structure.atoms import Atoms
from pyiron_atomistic.thermodynamics.interfacemethod import half_velocity, \
    freeze_one_half, fix_iso, fix_z_dir, create_job_template, check_diamond


class TestHessianJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(
            os.path.join(cls.file_location, "interface")
        )
        cls.project.remove_jobs_silently(recursive=True)
        cls.structure = Atoms(
            positions=[[0, 0, 0], [1, 1, 1]], elements=["Fe", "Fe"], cell=2 * np.eye(3)
        )
        cls.temperature = 500.0
        cls.cpu_cores = 2
        potential = pandas.DataFrame({
            'Name': ['Fe Morse'],
            'Filename': [[]],
            'Model': ['Morse'],
            'Species': [['Fe']],
            'Config': [['atom_style full\n',
                        'pair_coeff 1 2 morse 0.019623 1.8860 3.32833\n']]
        })
        project_dict = {
            "job_type": "Lammps",
            "project": cls.project,
            "potential": potential,
            "cpu_cores": cls.cpu_cores
        }
        cls.job = create_job_template(
            job_name="test_fix_iso",
            structure=cls.structure,
            project_parameter=project_dict
        )
        cls.job.calc_md(
            temperature=cls.temperature,
            temperature_damping_timescale=100.0,
            pressure=0.0,
            pressure_damping_timescale=1000.0
        )

    def test_freeze_one_half(self):
        structure_freeze = freeze_one_half(self.structure.copy())
        self.assertTrue(all(structure_freeze.selective_dynamics[0]))
        self.assertFalse(all(structure_freeze.selective_dynamics[1]))

    def test_job_creation(self):
        job = self.job.copy()
        self.assertTrue(isinstance(job, pyiron_atomistic.lammps.lammps.Lammps))
        self.assertEqual(job.server.cores, self.cpu_cores)

    def test_fix_iso(self):
        job = self.job.copy()
        md_str = job.input.control["fix___ensemble"]
        job_fixed = fix_iso(job=job)
        self.assertEqual(job_fixed.input.control["fix___ensemble"], md_str + " couple xyz")

    def test_half_velocity(self):
        job = self.job.copy()
        job_half = half_velocity(job=job, temperature=self.temperature)
        self.assertEqual(float(job_half.input.control['velocity'].split()[2]), self.temperature)

    def test_fix_z_dir(self):
        job = self.job.copy()
        job.calc_md(
            temperature=self.temperature,
            temperature_damping_timescale=100.0,
            pressure=[0.0, 0.0, 0.0],
            pressure_damping_timescale=1000.0
        )
        md_str = "all npt temp " + str(self.temperature) + " " + str(self.temperature) + " 0.1 z 0.0 0.0 1.0"
        job_z = fix_z_dir(job=job)
        self.assertEqual(md_str, job_z.input.control["fix___ensemble"])

    def test_check_diamond(self):
        al_fcc = ase_to_pyiron(bulk("Al", cubic=True))
        fe_bcc = ase_to_pyiron(bulk("Fe", cubic=True))
        ti_hcp = ase_to_pyiron(bulk("Ti", orthorhombic=True))
        si_dia = ase_to_pyiron(bulk("Si", cubic=True))
        self.assertFalse(check_diamond(structure=al_fcc))
        self.assertFalse(check_diamond(structure=fe_bcc))
        self.assertFalse(check_diamond(structure=ti_hcp))
        self.assertTrue(check_diamond(structure=si_dia))
