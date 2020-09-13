# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import posixpath

import numpy as np
from molmod.units import *

from pyiron.atomistics.structure.atoms import Atoms
from pyiron_base import ProjectHDFio, Project

from pyiron.quickff.quickff import QuickFFInput
from pyiron.quickff.quickff import QuickFF

class TestQuickFF(unittest.TestCase):
        """
        Tests the pyiron.quickff.quickff QuickFF class
        """

        @classmethod
        def setUpClass(cls):
            cls.execution_path = os.path.dirname(os.path.abspath(__file__))
            cls.project = Project(os.path.join(cls.execution_path, "test_quickff"))
            cls.job = cls.project.create_job("QuickFF", "trial")
            cls.job_complete = QuickFF(
                project=ProjectHDFio(project=cls.project, file_name="quickff_complete"),
                job_name="quickff_complete",
            )

        @classmethod
        def tearDownClass(cls):
            cls.execution_path = os.path.dirname(os.path.abspath(__file__))
            project = Project(os.path.join(cls.execution_path, "test_quickff"))
            project.remove_jobs_silently(recursive=True)
            project.remove(enable=True)

        def setUp(self):
            if self.project.load('trial') is not None:
                self.project.remove_job('trial')
            self.job = self.project.create_job("QuickFF", 'trial')

            if self.project.load('quickff_complete') is not None:
                self.project.remove_job("quickff_complete")
            self.job_complete = QuickFF(
                project=ProjectHDFio(project=self.project, file_name="quickff_complete"),
                job_name="quickff_complete",
            )

        def test_init(self):
            self.assertEqual(self.job.__name__, "QuickFF")
            self.assertIsInstance(self.job.input, QuickFFInput)
            self.assertEqual(self.job.ffatypes, None)
            self.assertEqual(self.job.ffatype_ids, None)
            self.assertEqual(self.job.aiener, None)
            self.assertEqual(self.job.aigrad, None)
            self.assertEqual(self.job.aihess, None)
            self.assertEqual(self.job.fn_ei, None)
            self.assertEqual(self.job.fn_vdw, None)


        def test_input(self):
            self.assertIsInstance(self.job.input['fn_yaff'], (str))
            self.assertIsInstance(self.job.input['fn_charmm22_prm'], (type(None)))
            self.assertIsInstance(self.job.input['fn_charmm22_psf'], (type(None)))
            self.assertIsInstance(self.job.input['fn_sys'], (str))
            self.assertIsInstance(self.job.input['plot_traj'], (type(None)))
            self.assertIsInstance(self.job.input['xyz_traj'], (bool))
            self.assertIsInstance(self.job.input['fn_traj'], (type(None)))
            self.assertIsInstance(self.job.input['log_level'], (str))
            self.assertIsInstance(self.job.input['log_file'], (str))
            self.assertIsInstance(self.job.input['program_mode'], (str))
            self.assertIsInstance(self.job.input['only_traj'], (str))
            self.assertIsInstance(self.job.input['ffatypes'], (type(None)))
            self.assertIsInstance(self.job.input['ei'], (type(None)))
            self.assertIsInstance(self.job.input['ei_rcut'], (type(None)))
            self.assertIsInstance(self.job.input['vdw'], (type(None)))
            self.assertIsInstance(self.job.input['vdw_rcut'], (float))
            self.assertIsInstance(self.job.input['covres'], (type(None)))
            self.assertIsInstance(self.job.input['excl_bonds'], (type(None)))
            self.assertIsInstance(self.job.input['excl_bends'], (type(None)))
            self.assertIsInstance(self.job.input['excl_dihs'], (type(None)))
            self.assertIsInstance(self.job.input['excl_oopds'], (type(None)))
            self.assertIsInstance(self.job.input['do_hess_mass_weighting'], (bool))
            self.assertIsInstance(self.job.input['do_hess_negfreq_proj'], (bool))
            self.assertIsInstance(self.job.input['do_cross_svd'], (bool))
            self.assertIsInstance(self.job.input['pert_traj_tol'], (float))
            self.assertIsInstance(self.job.input['pert_traj_energy_noise'], (type(None)))
            self.assertIsInstance(self.job.input['cross_svd_rcond'], (float))
            self.assertIsInstance(self.job.input['do_bonds'], (bool))
            self.assertIsInstance(self.job.input['do_bends'], (bool))
            self.assertIsInstance(self.job.input['do_dihedrals'], (bool))
            self.assertIsInstance(self.job.input['do_oops'], (bool))
            self.assertIsInstance(self.job.input['do_cross_ASS'], (bool))
            self.assertIsInstance(self.job.input['do_cross_ASA'], (bool))
            self.assertIsInstance(self.job.input['do_cross_DSS'], (bool))
            self.assertIsInstance(self.job.input['do_cross_DSD'], (bool))
            self.assertIsInstance(self.job.input['do_cross_DAA'], (bool))
            self.assertIsInstance(self.job.input['do_cross_DAD'], (bool))
            self.assertIsInstance(self.job.input['consistent_cross_rvs'], (bool))
            self.assertIsInstance(self.job.input['remove_dysfunctional_cross'], (bool))
            self.assertIsInstance(self.job.input['bond_term'], (str))
            self.assertIsInstance(self.job.input['bend_term'], (str))
            self.assertIsInstance(self.job.input['do_squarebend'], (bool))
            self.assertIsInstance(self.job.input['do_bendclin'], (bool))
            self.assertIsInstance(self.job.input['do_sqoopdist_to_oopdist'], (bool))

            self.job.set_ei(posixpath.join(self.execution_path, "../static/quickff_test_files/pars_mbisgauss.txt"))
            self.assertEqual(self.job.input['ei'], 'pars_mbisgauss.txt')
            self.assertEqual(self.job.fn_ei, posixpath.join(self.execution_path, "../static/quickff_test_files/pars_mbisgauss.txt"))
            self.job.set_vdw(posixpath.join(self.execution_path, "../static/quickff_test_files/pars_vdw.txt"))
            self.assertEqual(self.job.input['vdw'], 'pars_vdw.txt')
            self.assertEqual(self.job.fn_vdw, posixpath.join(self.execution_path, "../static/quickff_test_files/pars_vdw.txt"))


        def test_set_structure(self):
            self.assertEqual(self.job.structure, None)
            self.job.read_abinitio(posixpath.join(self.execution_path, "../static/quickff_test_files/input.fchk"))
            self.assertIsInstance(self.job.structure, (Atoms))

            self.assertIsInstance(self.job.aiener, (float))
            self.assertTrue(
                np.array_equal(self.job.aigrad.shape, (len(self.job.structure), 3))
            )
            self.assertTrue(
                np.array_equal(self.job.aihess.shape, (len(self.job.structure), 3, len(self.job.structure), 3))
            )

            self.job.detect_ffatypes(ffatype_level='low')
            ffatypes = self.job.ffatypes
            ffatype_ids = self.job.ffatype_ids

            full_list = [ffatypes[i] for i in ffatype_ids]
            self.job.detect_ffatypes(ffatypes=full_list)
            self.assertCountEqual(self.job.ffatypes, ffatypes)
            self.assertCountEqual(self.job.ffatype_ids, ffatype_ids)

            rules =[
                ('H', '1'),
                ('C', '6'),
            ]
            self.job.detect_ffatypes(ffatype_rules=rules)
            self.assertCountEqual(self.job.ffatypes, ffatypes)
            self.assertCountEqual(self.job.ffatype_ids, ffatype_ids)

            self.assertRaises(IOError, self.job.detect_ffatypes, ffatype_rules=rules, ffatype_level='high')

        def test_run_complete(self):
            self.job_complete.read_abinitio(
                posixpath.join(self.execution_path, "../static/quickff_test_files/input.fchk")
            )

            self.job_complete.set_ei(
                posixpath.join(self.execution_path, "../static/quickff_test_files/pars_mbisgauss.txt")
            )

            self.job_complete.detect_ffatypes(ffatype_level='low')

            file_directory = posixpath.join(
                self.execution_path, '../static/quickff_test_files/'
            )
            self.job_complete.restart_file_list.append(
                posixpath.join(file_directory, "system.chk")
            )
            self.job_complete.restart_file_list.append(
                posixpath.join(file_directory, "pars_cov.txt")
            )

            self.job_complete.run(run_mode="manual")
            self.job_complete.status.collect = True
            self.job_complete.run()
            nodes = [
                'bend',
                'bond',
                'cross',
                'oopdist',
                'torsion',
            ]
            with self.job_complete.project_hdf5.open("output/generic") as h_gen:
                hdf_nodes = h_gen.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))

            nodes = [
                'bonds',
                'ffatype_ids',
                'ffatypes',
                'numbers',
                'pos',
                'rvecs'
            ]
            with self.job_complete.project_hdf5.open("output/system") as h_gen:
                hdf_nodes = h_gen.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))


if __name__ == "__main__":
    unittest.main()
