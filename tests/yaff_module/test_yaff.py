# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os,sys
import posixpath

from molmod.units import *

from pyiron import ase_to_pyiron
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.base.project.generic import Project
from pyiron.project import Project as Pr

from pyiron.yaff.yaff import YaffInput
from pyiron.yaff.yaff import Yaff

class TestYaff(unittest.TestCase):
        """
        Tests the pyiron.yaff.yaff Yaff class
        """

        @classmethod
        def setUpClass(cls):
            cls.execution_path = os.path.dirname(os.path.abspath(__file__))
            cls.project = Project(os.path.join(cls.execution_path, "test_yaff"))
            cls.job = cls.project.create_job("Yaff", "trial")
            cls.job_complete = Yaff(
                project=ProjectHDFio(project=cls.project, file_name="yaff_complete"),
                job_name="yaff_complete",
            )

        @classmethod
        def tearDownClass(cls):
            cls.execution_path = os.path.dirname(os.path.abspath(__file__))
            project = Project(os.path.join(cls.execution_path, "test_yaff"))
            project.remove_jobs(recursive=True)
            project.remove(enable=True)

        def setUp(self):
            self.project.remove_job("trial")
            self.job = self.project.create_job("Yaff", "trial")

            self.project.remove_job("yaff_complete")
            self.job_complete = Yaff(
                project=ProjectHDFio(project=self.project, file_name="yaff_complete"),
                job_name="yaff_complete",
            )

        def test_init(self):
            self.assertEqual(self.job.__name__, "Yaff")
            self.assertIsInstance(self.job.input, YaffInput)
            self.assertEqual(self.job.ffatypes, None)
            self.assertEqual(self.job.ffatype_ids, None)
            self.assertEqual(self.job.enhanced, None)


        def test_input(self):
            self.assertIsInstance(self.job.input['rcut'], (float))
            self.assertIsInstance(self.job.input['alpha_scale'], (float))
            self.assertIsInstance(self.job.input['gcut_scale'], (float))
            self.assertIsInstance(self.job.input['smooth_ei'], (bool))
            self.assertIsInstance(self.job.input['gpos_rms'], (float))
            self.assertIsInstance(self.job.input['dpos_rms'], (float))
            self.assertIsInstance(self.job.input['grvecs_rms'], (float))
            self.assertIsInstance(self.job.input['drvecs_rms'], (float))
            self.assertIsInstance(self.job.input['hessian_eps'], (float))
            self.assertIsInstance(self.job.input['timestep'], (float))
            self.assertIsInstance(self.job.input['temp'], (type(None)))
            self.assertIsInstance(self.job.input['press'], (type(None)))
            self.assertIsInstance(self.job.input['timecon_thermo'], (float))
            self.assertIsInstance(self.job.input['timecon_baro'], (float))
            self.assertIsInstance(self.job.input['nsteps'], (int))
            self.assertIsInstance(self.job.input['h5step'], (int))
            self.job.input['temp'] = 300*kelvin
            self.assertEqual(self.job.input['temp'], 300*kelvin)
            self.job.input['press'] = 1e5*pascal
            self.assertEqual(self.job.input['press'], 1e5*pascal)
            self.job.ffatypes = ['H', 'N']
            self.assertEqual(self.job.ffatypes, ['H', 'N'])
            self.job.ffatype_ids = [0,1]
            self.assertEqual(self.job.ffatype_ids, [0,1])
            self.job.enhanced = {}
            self.assertEqual(self.job.enhanced, {})


        def test_set_mtd(self):
            assert True

        def test_set_us(self):
            assert True


        def test_calc_static(self):
            self.job.calc_static()
            self.assertEqual(self.job.input['jobtype'], 'sp')

        def test_calc_minimize(self):
            self.job.calc_minimize(gpos_tol=1e-8, dpos_tol=1e-6, grvecs_tol=1e-8, drvecs_tol=1e-6, max_iter=2000, n_print=10)
            self.assertEqual(self.job.input['jobtype'], 'opt')
            self.assertEqual(self.job.input['gpos_rms'], 1e-8)
            self.assertEqual(self.job.input['dpos_rms'], 1e-6)
            self.assertEqual(self.job.input['grvecs_rms'], 1e-8)
            self.assertEqual(self.job.input['drvecs_rms'], 1e-6)
            self.assertEqual(self.job.input['nsteps'], 2000)
            self.assertEqual(self.job.input['h5step'], 10)

        def test_calc_minimize_cell(self):
            self.job.calc_minimize(cell=True)
            self.assertEqual(self.job.input['jobtype'], 'opt_cell')

        def test_calc_md_nve(self):
            self.job.calc_md(time_step=1.0*femtosecond)
            self.assertEqual(self.job.input['jobtype'], 'nve')
            self.assertEqual(self.job.input['timestep'], 1.0*femtosecond)

        def test_calc_md_nvt(self):
            self.job.calc_md(temperature=300*kelvin, time_step=1.0*femtosecond,timecon_thermo=100.0*femtosecond
            )
            self.assertEqual(self.job.input['jobtype'], 'nvt')
            self.assertEqual(self.job.input['temp'], 300*kelvin)
            self.assertEqual(self.job.input['timestep'], 1.0*femtosecond)
            self.assertEqual(self.job.input['timecon_thermo'], 100.0*femtosecond)

        def test_calc_md_npt(self):
            self.job.calc_md(temperature=300*kelvin, pressure=1*bar, time_step=1.0*femtosecond,
                        timecon_thermo=100.0*femtosecond,timecon_baro=1000.0*femtosecond
            )
            self.assertEqual(self.job.input['jobtype'], 'npt')
            self.assertEqual(self.job.input['temp'], 300*kelvin)
            self.assertEqual(self.job.input['press'], 1*bar)
            self.assertEqual(self.job.input['timestep'], 1.0*femtosecond)
            self.assertEqual(self.job.input['timecon_thermo'], 100.0*femtosecond)
            self.assertEqual(self.job.input['timecon_baro'], 1000.0*femtosecond)

        def test_set_structure(self):
            self.assertEqual(self.job.structure, None)
            self.job.load_chk(posixpath.join(self.execution_path, "../static/yaff_test_files/sp/input.chk"))
            self.assertIsInstance(self.job.structure, (Atoms))

            ffatypes = self.job.ffatypes
            ffatype_ids = self.job.ffatype_ids

            self.job.detect_ffatypes(ffatype_level='low')
            self.assertCountEqual(self.job.ffatypes, ffatypes)
            self.assertCountEqual(self.job.ffatype_ids, ffatype_ids)

            full_list = [ffatypes[i] for i in ffatype_ids]
            self.job.detect_ffatypes(ffatypes=full_list)
            self.assertCountEqual(self.job.ffatypes, ffatypes)
            self.assertCountEqual(self.job.ffatype_ids, ffatype_ids)

            rules =[
                ('H', '1 & =1%8'),
                ('O', '8 & =2%1'),
            ]
            self.job.detect_ffatypes(ffatype_rules=rules)
            self.assertCountEqual(self.job.ffatypes, ffatypes)
            self.assertCountEqual(self.job.ffatype_ids, ffatype_ids)

            self.assertRaises(IOError, self.job.detect_ffatypes, ffatype_rules=rules, ffatype_level='high')

        '''
        def test_run_sp_complete(self):
            self.job_complete.structure = ase_to_pyiron(read_gaussian_out(
                posixpath.join(self.execution_path, "../static/gaussian_test_files/sp/input.log")
            ))
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

        def test_run_bsse_complete(self):
            self.job_complete.structure = ase_to_pyiron(read_gaussian_out(
                posixpath.join(self.execution_path, "../static/gaussian_test_files/bsse/input.log")
            ))
            self.job_complete.input['lot'] = 'B3LYP'
            self.job_complete.input['basis_set'] = '6-31+G(d)'
            self.job_complete.calc_minimize()
            file_directory = posixpath.join(
                self.execution_path, '../static/gaussian_test_files/bsse'
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

            nodes = [
                "energy_tot_corrected",
                "bsse_correction",
                "sum_of_fragments",
                "complexation_energy_raw",
                "complexation_energy_corrected",
            ]
            with self.job_complete.project_hdf5.open("output/structure/bsse") as h_gen:
                hdf_nodes = h_gen.list_nodes()
                self.assertTrue(all([node in hdf_nodes for node in nodes]))


        def test_run_empdisp_complete(self):
            self.job_complete.structure = ase_to_pyiron(read_gaussian_out(
                posixpath.join(self.execution_path, "../static/gaussian_test_files/empdisp/input.log")
            ))
            self.job_complete.input['lot'] = 'B3LYP'
            self.job_complete.input['basis_set'] = '6-311+G*'
            self.job_complete.input['settings'] = {'EmpiricalDispersion':['GD3']}
            self.job_complete.calc_static()

            file_directory = posixpath.join(
                self.execution_path, '../static/gaussian_test_files/empdisp'
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

            import io
            from contextlib import redirect_stdout

            with io.StringIO() as buf, redirect_stdout(buf):
                self.job_complete.log()
                output = buf.getvalue()

            self.assertTrue("R6Disp" in output)
        '''


if __name__ == "__main__":
    unittest.main()
