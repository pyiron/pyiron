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
            self.job_complete.input['lot'] = 'B3LYP'
            self.job_complete.input['basis_set'] = '6-31+G(d)'
            self.job_complete.calc_static()
            file_directory = posixpath.join(
                self.execution_path, '../static/gaussian_test_files/sp'
            )
            self.job_complete.restart_file_list.append(
                posixpath.join(file_directory, "input.log")
            )
            self.job_complete.restart_file_list.append(
                posixpath.join(file_directory, "input.fchk")
            )
            self.job_complete.run(run_mode="manual")
            self.job_complete.status.collect = True
            self.job_complete.run()
            nodes = [
                "positions",
                "energy_tot",
            ]
            with self.job_complete.project_hdf5.open("output/generic") as h_gen:
                hdf_nodes = h_gen.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))

            nodes = [
                "charges",
                "dipole",
                "masses",
                "numbers",
                "positions",
            ]
            with self.job_complete.project_hdf5.open("output/structure") as h_gen:
                hdf_nodes = h_gen.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))

            nodes = [
                "alpha_orbital_e",
                "beta_orbital_e",
                "n_alpha_electrons",
                "n_basis_functions",
                "n_beta_electrons",
                "n_electrons",
                "scf_density",
                "spin_scf_density",
            ]
            with self.job_complete.project_hdf5.open("output/structure/dft") as h_gen:
                hdf_nodes = h_gen.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))
                
