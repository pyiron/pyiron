# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import posixpath

import numpy as np
from molmod.units import *

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.base.project.generic import Project

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
        if self.project.load('trial') is not None:
            self.project.remove_job("trial")
        self.job = self.project.create_job("Yaff", "trial")

        if self.project.load('yaff_complete') is not None:
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
        self.assertIsInstance(self.job.input['rcut'], float)
        self.assertIsInstance(self.job.input['alpha_scale'], float)
        self.assertIsInstance(self.job.input['gcut_scale'], float)
        self.assertIsInstance(self.job.input['smooth_ei'], bool)
        self.assertIsInstance(self.job.input['gpos_rms'], float)
        self.assertIsInstance(self.job.input['dpos_rms'], float)
        self.assertIsInstance(self.job.input['grvecs_rms'], float)
        self.assertIsInstance(self.job.input['drvecs_rms'], float)
        self.assertIsInstance(self.job.input['hessian_eps'], float)
        self.assertIsInstance(self.job.input['timestep'], float)
        self.assertIsInstance(self.job.input['temp'], (type(None)))
        self.assertIsInstance(self.job.input['press'], (type(None)))
        self.assertIsInstance(self.job.input['timecon_thermo'], float)
        self.assertIsInstance(self.job.input['timecon_baro'], float)
        self.assertIsInstance(self.job.input['nsteps'], int)
        self.assertIsInstance(self.job.input['h5step'], int)
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
        self.job.set_mtd([('torsion', [0, 1, 2, 3])], 1*kjmol, 5.*deg, 40, stride=20, temp=300.)
        np.testing.assert_array_equal(self.job.enhanced['ickinds'], np.array(['torsion'], dtype='S22'))
        np.testing.assert_array_equal(self.job.enhanced['icindices'], np.array([[1, 2, 3, 4]]))  #  plumed starts counting from 1
        np.testing.assert_array_equal(self.job.enhanced['height'], np.array([1*kjmol]))
        np.testing.assert_array_equal(self.job.enhanced['sigma'], np.array([5*deg]))
        self.assertEqual(self.job.enhanced['pace'], 40)
        self.assertEqual(self.job.enhanced['file'], 'HILLS')
        self.assertEqual(self.job.enhanced['file_colvar'], 'COLVAR')
        self.assertEqual(self.job.enhanced['stride'], 20)
        self.assertEqual(self.job.enhanced['temp'], 300.)

    def test_set_us(self):
        self.job.set_us([('torsion', [0, 1, 2, 3])], 10*kjmol, 5.*deg, fn_colvar='COLVAR', temp=300.)
        np.testing.assert_array_equal(self.job.enhanced['ickinds'], np.array(['torsion'], dtype='S22'))
        np.testing.assert_array_equal(self.job.enhanced['icindices'], np.array([[1, 2, 3, 4]]))  #  plumed starts counting from 1
        np.testing.assert_array_equal(self.job.enhanced['kappa'], np.array([10*kjmol]))
        np.testing.assert_array_equal(self.job.enhanced['loc'], np.array([5*deg]))
        self.assertEqual(self.job.enhanced['file_colvar'], 'COLVAR')
        self.assertEqual(self.job.enhanced['stride'], 10)
        self.assertEqual(self.job.enhanced['temp'], 300.)

    def test_calc_static(self):
        self.job.calc_static()
        self.assertEqual(self.job.input['jobtype'], 'sp')

    def test_calc_minimize(self):
        self.job.calc_minimize(
            gpos_tol=1e-8,
            dpos_tol=1e-6,
            grvecs_tol=1e-8,
            drvecs_tol=1e-6,
            max_iter=2000,
            n_print=10
        )
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
        self.job.calc_md(
            temperature=300*kelvin,
            time_step=1.0*femtosecond,
            timecon_thermo=100.0*femtosecond
        )
        self.assertEqual(self.job.input['jobtype'], 'nvt')
        self.assertEqual(self.job.input['temp'], 300*kelvin)
        self.assertEqual(self.job.input['timestep'], 1.0*femtosecond)
        self.assertEqual(self.job.input['timecon_thermo'], 100.0*femtosecond)

    def test_calc_md_npt(self):
        self.job.calc_md(
            temperature=300*kelvin,
            pressure=1*bar,
            time_step=1.0*femtosecond,
            timecon_thermo=100.0*femtosecond,
            timecon_baro=1000.0*femtosecond
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
        self.assertIsInstance(self.job.structure, Atoms)

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

    def test_run_sp_complete(self):
        self.job_complete.load_chk(
            posixpath.join(self.execution_path, "../static/yaff_test_files/sp/system.chk")
        )

        ffpars = """
        BONDHARM:UNIT  K kjmol/A**2
        BONDHARM:UNIT  R0 A
        BONDHARM:PARS  H O  4.9657739952e+03  9.6881765966e-01

        BENDAHARM:UNIT  K kjmol/rad**2
        BENDAHARM:UNIT  THETA0 deg
        BENDAHARM:PARS  H O H  3.0369893169e+02  9.6006623652e+01

        FIXQ:UNIT Q0 e
        FIXQ:UNIT P e
        FIXQ:UNIT R angstrom
        FIXQ:SCALE 1 1.0
        FIXQ:SCALE 2 1.0
        FIXQ:SCALE 3 1.0
        FIXQ:DIELECTRIC 1.0
        FIXQ:ATOM   H  0.4505087957  0.7309000000
        FIXQ:ATOM   O -0.9012059960  1.1325000000
        """

        self.job_complete.input['ffpars'] = ffpars
        self.job_complete.calc_static()
        file_directory = posixpath.join(
            self.execution_path, '../static/yaff_test_files/sp'
        )
        self.job_complete.restart_file_list.append(
            posixpath.join(file_directory, "output.h5")
        )

        self.job_complete.run(run_mode="manual")
        self.job_complete.status.collect = True
        self.job_complete.run()
        nodes = [
            "energy_pot",
            "positions",
            "cells",
        ]
        with self.job_complete.project_hdf5.open("output/generic") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

        nodes = [
            "ffatype_ids",
            "ffatypes",
            "masses",
            "numbers",
            "positions",
        ]
        with self.job_complete.project_hdf5.open("output/structure") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

    def test_run_minimize_complete(self):
        self.job_complete.load_chk(
            posixpath.join(self.execution_path, "../static/yaff_test_files/opt/opt.chk")
        )

        ffpars = """
        BONDHARM:UNIT  K kjmol/A**2
        BONDHARM:UNIT  R0 A
        BONDHARM:PARS  H O  4.9657739952e+03  9.6881765966e-01

        BENDAHARM:UNIT  K kjmol/rad**2
        BENDAHARM:UNIT  THETA0 deg
        BENDAHARM:PARS  H O H  3.0369893169e+02  9.6006623652e+01

        FIXQ:UNIT Q0 e
        FIXQ:UNIT P e
        FIXQ:UNIT R angstrom
        FIXQ:SCALE 1 1.0
        FIXQ:SCALE 2 1.0
        FIXQ:SCALE 3 1.0
        FIXQ:DIELECTRIC 1.0
        FIXQ:ATOM   H  0.4505087957  0.7309000000
        FIXQ:ATOM   O -0.9012059960  1.1325000000
        """

        self.job_complete.input['ffpars'] = ffpars
        self.job_complete.calc_minimize()
        file_directory = posixpath.join(
            self.execution_path, '../static/yaff_test_files/opt'
        )
        self.job_complete.restart_file_list.append(
            posixpath.join(file_directory, "output.h5")
        )

        self.job_complete.run(run_mode="manual")
        self.job_complete.status.collect = True
        self.job_complete.run()
        nodes = [
            "energy_pot",
            "positions",
            "cells",
            "steps",
            "volume",
        ]
        with self.job_complete.project_hdf5.open("output/generic") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

        nodes = [
            "ffatype_ids",
            "ffatypes",
            "masses",
            "numbers",
            "positions",
        ]
        with self.job_complete.project_hdf5.open("output/structure") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

    def test_run_md_complete(self):
        self.job_complete.load_chk(
            posixpath.join(self.execution_path, "../static/yaff_test_files/md/system.chk")
        )

        ffpars = """
        BONDHARM:UNIT  K kjmol/A**2
        BONDHARM:UNIT  R0 A
        BONDHARM:PARS  H O  4.9657739952e+03  9.6881765966e-01

        BENDAHARM:UNIT  K kjmol/rad**2
        BENDAHARM:UNIT  THETA0 deg
        BENDAHARM:PARS  H O H  3.0369893169e+02  9.6006623652e+01

        FIXQ:UNIT Q0 e
        FIXQ:UNIT P e
        FIXQ:UNIT R angstrom
        FIXQ:SCALE 1 1.0
        FIXQ:SCALE 2 1.0
        FIXQ:SCALE 3 1.0
        FIXQ:DIELECTRIC 1.0
        FIXQ:ATOM   H  0.4505087957  0.7309000000
        FIXQ:ATOM   O -0.9012059960  1.1325000000
        """

        self.job_complete.input['ffpars'] = ffpars
        self.job_complete.calc_md(temperature=300*kelvin, nsteps=1000, time_step=1.0*femtosecond)
        file_directory = posixpath.join(
            self.execution_path, '../static/yaff_test_files/md'
        )
        self.job_complete.restart_file_list.append(
            posixpath.join(file_directory, "output.h5")
        )

        self.job_complete.run(run_mode="manual")
        self.job_complete.status.collect = True
        self.job_complete.run()
        nodes = [
            "energy_pot",
            "energy_kin",
            "energy_tot",
            "energy_cons",
            "positions",
            "cells",
            "steps",
            "volume",
            "temperature",
            "time",
        ]
        with self.job_complete.project_hdf5.open("output/generic") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

        nodes = [
            "ffatype_ids",
            "ffatypes",
            "masses",
            "numbers",
            "positions",
        ]
        with self.job_complete.project_hdf5.open("output/structure") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

    def test_run_metadynamics_complete(self):
        self.job_complete.load_chk(
            posixpath.join(self.execution_path, "../static/yaff_test_files/metadynamics/system.chk")
        )

        ffpars = """
        # BONDHARM
        #---------
        BONDHARM:UNIT  K kjmol/A**2
        BONDHARM:UNIT  R0 A
        BONDHARM:PARS      H1_n        N3  4.1645524453e+03  1.0246931776e+00

        # BENDAHARM
        #----------
        BENDAHARM:UNIT  K kjmol/rad**2
        BENDAHARM:UNIT  THETA0 deg
        BENDAHARM:PARS      H1_n        N3      H1_n  3.3770180977e+02  1.0417265544e+02

        # SQOOPDIST
        #----------
        SQOOPDIST:UNIT  K kjmol/A**4
        SQOOPDIST:UNIT  D0 A**2
        SQOOPDIST:PARS      H1_n      H1_n      H1_n        N3  2.9152308037e+01  5.6225686721e+00

        # Cross
        #------
        Cross:UNIT  KSS kjmol/angstrom**2
        Cross:UNIT  KBS0 kjmol/(angstrom*rad)
        Cross:UNIT  KBS1 kjmol/(angstrom*rad)
        Cross:UNIT  R0 angstrom
        Cross:UNIT  R1 angstrom
        Cross:UNIT  THETA0 deg
        Cross:PARS      H1_n        N3      H1_n  -4.0297219155e+01   1.2044509670e+02   1.2044509670e+02  1.0217511358e+00  1.0217511358e+00  1.0594536883e+02

        #Fixed charges
        #---------------
        FIXQ:UNIT Q0 e
        FIXQ:UNIT P e
        FIXQ:UNIT R angstrom
        FIXQ:SCALE 1 1.0
        FIXQ:SCALE 2 1.0
        FIXQ:SCALE 3 1.0
        FIXQ:DIELECTRIC 1.0

        # Atomic parameters
        # ----------------------------------------------------
        # KEY        label  Q_0A              R_A
        # ----------------------------------------------------
        FIXQ:ATOM       N3  0.0000000000  1.1039000000
        FIXQ:ATOM     H1_n  0.0000000000  0.7309000000
        # Bond parameters
        # ----------------------------------------------------
        # KEY         label0   label1           P_AB
        # ----------------------------------------------------
        FIXQ:BOND      H1_n        N3   0.3816559688
        """

        self.job_complete.input['ffpars'] = ffpars
        self.job_complete.calc_md(temperature=300*kelvin, nsteps=1000, time_step=0.5*femtosecond, n_print=1)
        self.job_complete.set_mtd([('torsion', [0, 1, 2, 3])], 1*kjmol, 5.*deg, 40)
        file_directory = posixpath.join(
            self.execution_path, '../static/yaff_test_files/metadynamics'
        )
        self.job_complete.restart_file_list.append(
            posixpath.join(file_directory, "output.h5")
        )
        self.job_complete.restart_file_list.append(
            posixpath.join(file_directory, "COLVAR")
        )

        self.job_complete.run(run_mode="manual")
        self.job_complete.status.collect = True
        self.job_complete.run()

        nodes = [
            "bias",
            "cv",
            "time",
        ]
        with self.job_complete.project_hdf5.open("output/enhanced") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

        nodes = [
            "energy_pot",
            "energy_kin",
            "energy_tot",
            "energy_cons",
            "positions",
            "cells",
            "steps",
            "volume",
            "temperature",
            "time",
        ]
        with self.job_complete.project_hdf5.open("output/generic") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))

        nodes = [
            "ffatype_ids",
            "ffatypes",
            "masses",
            "numbers",
            "positions",
        ]
        with self.job_complete.project_hdf5.open("output/structure") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))


if __name__ == "__main__":
    unittest.main()
