# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import posixpath

from ase.io import read

from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.base.project.generic import Project

from pyiron.gaussian.gaussian import GaussianInput
from pyiron.gaussian.gaussian import Gaussian



class TestGaussian(unittest.TestCase):
        """
        Tests the pyiron.pyiron.gaussian Gaussian class
        """

        @classmethod
        def setUpClass(cls):
            cls.execution_path = os.path.dirname(os.path.abspath(__file__))
            cls.project = Project(os.path.join(cls.execution_path, "test_gaussian"))
            cls.job = cls.project.create_job("Gaussian", "trial")
            cls.job_complete = Gaussian(
                project=ProjectHDFio(project=cls.project, file_name="gaussian_complete"),
                job_name="gaussian_complete",
            )
            log_file = posixpath.join(
                cls.execution_path, "../static/gaussian_test_files/sp/input.log"
            )
            cls.job_complete.structure = read(log_file)

        @classmethod
        def tearDownClass(cls):
            cls.execution_path = os.path.dirname(os.path.abspath(__file__))
            project = Project(os.path.join(cls.execution_path, "test_gaussian"))
            project.remove_jobs(recursive=True)
            project.remove(enable=True)

        def setUp(self):
            self.job.structure = None

        def test_init(self):
            self.assertEqual(self.job.__name__, "Gaussian")
            self.assertIsInstance(self.job.input, GaussianInput)


        def test_input(self):
            self.assertIsInstance(self.job.input['lot'], (string))
            self.assertIsInstance(self.job.input['basis_set'], (string))
            self.assertIsInstance(self.job.input['spin_mult'], (int))
            self.assertIsInstance(self.job.input['charge'], (int))
            # self.assertIsInstance(self.job.plane_wave_cutoff, (float, int))
            self.job.input['lot'] = 'B3LYP'
            self.assertEqual(self.job.input['lot'], 'B3LYP')
            self.job.input['basis_set'] = 'sto-3g'
            self.assertEqual(self.job.input['basis_set'], 'sto-3g')
            self.job.input['spin_mult'] = 1
            self.assertEqual(self.job.input['spin_mult'], 2)
            self.job.input['charge'] = 1
            self.assertEqual(self.job.input['charge'], 1)


        def test_calc_static(self):
            self.job.calc_static(
                electronic_steps=90,
            )
            self.assertEqual(self.job.input['jobtype'], 'sp')
            self.assertEqual(self.job.input['settings'], {'SCF':['MaxCycle=90']})


        def test_calc_minimize(self):
            self.job.calc_minimize(
                electronic_steps=80,
                ionic_steps=90,
                algorithm=None,
                ionic_forces='tight',
            )
            self.assertEqual(self.job.input['jobtype'], 'opt(MaxCycles=90,tight)')
            self.assertEqual(self.job.input['settings'], {'SCF':['MaxCycle=80']})

        def test_set_structure(self):
            self.assertEqual(self.job.structure, None)
            atoms = create_atoms(elements=['H', 'H', 'O'], positions=[[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]])
            self.job.structure = atoms
            self.assertEqual(self.job.structure, atoms)
            self.job.structure = None
            self.assertEqual(self.job.structure, None)
            self.job.structure = atoms
            self.assertEqual(self.job.structure, atoms)


        def test_run_complete(self):
            self.job_complete.exchange_correlation_functional = "PBE"
            self.job_complete.set_occupancy_smearing(smearing="fermi", width=0.2)
            self.job_complete.calc_static()
            self.job_complete.set_convergence_precision(electronic_energy=1e-7)
            self.job_complete.write_electrostatic_potential = False
            self.assertEqual(self.job_complete.input.incar["SIGMA"], 0.2)
            self.assertEqual(self.job_complete.input.incar["LVTOT"], False)
            self.assertEqual(self.job_complete.input.incar["EDIFF"], 1e-7)
            file_directory = posixpath.join(
                self.execution_path, "../static/vasp_test_files/full_job_sample"
            )
            self.job_complete.restart_file_list.append(
                posixpath.join(file_directory, "vasprun.xml")
            )
            self.job_complete.restart_file_list.append(
                posixpath.join(file_directory, "OUTCAR")
            )
            self.job_complete.run(run_mode="manual")
            self.job_complete.status.collect = True
            self.job_complete.run()
            nodes = [
                "positions",
                "temperature",
                "energy_tot",
                "steps",
                "positions",
                "forces",
                "cells",
                "pressures",
            ]
            with self.job_complete.project_hdf5.open("output/generic") as h_gen:
                hdf_nodes = h_gen.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))
            nodes = [
                "energy_free",
                "energy_int",
                "energy_zero",
                "final_magmoms",
                "magnetization",
                "n_elect",
                "scf_dipole_mom",
                "scf_energy_free",
                "scf_energy_int",
                "scf_energy_zero",
            ]
            with self.job_complete.project_hdf5.open("output/generic/dft") as h_dft:
                hdf_nodes = h_dft.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))
            nodes = ["efermi", "eig_matrix", "k_points", "k_weights", "occ_matrix"]
            with self.job_complete.project_hdf5.open(
                "output/electronic_structure"
            ) as h_dft:
                hdf_nodes = h_dft.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))

            job_chg_den = self.job_complete.restart_from_charge_density(job_name="chg")
            self.assertEqual(job_chg_den.structure, self.job_complete.get_structure(-1))
            self.assertTrue(
                posixpath.join(self.job_complete.working_directory, "CHGCAR")
                in job_chg_den.restart_file_list
            )
            with job_chg_den.project_hdf5.open("output") as h_out:
                self.assertTrue(h_out.list_nodes() == [])
                self.assertTrue(h_out.list_groups() == [])

            with job_chg_den.project_hdf5.open("input") as h_in:
                self.assertFalse(h_in.list_nodes() == [])
                self.assertFalse(h_in.list_groups() == [])

            job_wave = self.job_complete.restart_from_wave_functions(job_name="wave")
            self.assertEqual(job_wave.structure, self.job_complete.get_structure(-1))
            self.assertTrue(
                posixpath.join(self.job_complete.working_directory, "WAVECAR")
                in job_wave.restart_file_list
            )
            with job_wave.project_hdf5.open("output") as h_out:
                self.assertTrue(h_out.list_nodes() == [])
                self.assertTrue(h_out.list_groups() == [])

            with job_wave.project_hdf5.open("input") as h_in:
                self.assertFalse(h_in.list_nodes() == [])
                self.assertFalse(h_in.list_groups() == [])

            job_chg_wave = self.job_complete.restart_from_wave_and_charge(
                job_name="chg_wave"
            )
            self.assertEqual(job_chg_wave.structure, self.job_complete.get_structure(-1))
            self.assertTrue(
                posixpath.join(self.job_complete.working_directory, "WAVECAR")
                in job_chg_wave.restart_file_list
            )
            self.assertTrue(
                posixpath.join(self.job_complete.working_directory, "CHGCAR")
                in job_chg_wave.restart_file_list
            )
            for key, val in job_chg_wave.restart_file_dict.items():
                self.assertTrue(key, val)
            with job_chg_wave.project_hdf5.open("output") as h_out:
                self.assertTrue(h_out.list_nodes() == [])
                self.assertTrue(h_out.list_groups() == [])

            with job_chg_wave.project_hdf5.open("input") as h_in:
                self.assertFalse(h_in.list_nodes() == [])
                self.assertFalse(h_in.list_groups() == [])
