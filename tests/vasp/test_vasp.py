# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import posixpath
from pyiron.atomistics.structure.atoms import CrystalStructure
from pyiron.vasp.base import Input, Output
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.base.project.generic import Project
from pyiron.vasp.potential import VaspPotentialSetter
from pyiron.vasp.vasp import Vasp
from pyiron.vasp.structure import read_atoms

__author__ = "Sudarsan Surendralal"


class TestVasp(unittest.TestCase):
    """
    Tests the pyiron.objects.hamilton.dft.vasp.Vasp class
    """

    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, "test_vasp"))
        cls.job = cls.project.create_job("Vasp", "trial")
        cls.job_complete = Vasp(
            project=ProjectHDFio(project=cls.project, file_name="vasp_complete"),
            job_name="vasp_complete",
        )
        poscar_file = posixpath.join(
            cls.execution_path, "../static/vasp_test_files/full_job_sample/POSCAR"
        )
        cls.job_complete.structure = read_atoms(poscar_file, species_from_potcar=True)

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, "test_vasp"))
        project.remove_jobs(recursive=True)
        project.remove(enable=True)

    def setUp(self):
        self.job.structure = None

    def test_init(self):
        self.assertEqual(self.job.__name__, "Vasp")
        self.assertEqual(self.job._sorted_indices, None)
        self.assertIsInstance(self.job.input, Input)
        self.assertIsInstance(self.job._output_parser, Output)
        self.assertIsInstance(self.job._potential, VaspPotentialSetter)
        self.assertTrue(self.job._compress_by_default)
        self.assertEqual(self.job.get_eddrmm_handling(), "not_converged")

    def test_eddrmm(self):
        self.job.set_eddrmm_handling("ignore")
        self.assertEqual(self.job.get_eddrmm_handling(), "ignore")
        self.job.set_eddrmm_handling("restart")
        self.assertEqual(self.job.get_eddrmm_handling(), "restart")
        self.job.set_eddrmm_handling()
        self.assertEqual(self.job.get_eddrmm_handling(), "not_converged")
        self.assertRaises(ValueError, self.job.set_eddrmm_handling, status="blah")

    def test_potential(self):
        self.assertEqual(self.job.potential, self.job._potential)

    def test_plane_wave_cutoff(self):
        self.assertIsInstance(self.job.plane_wave_cutoff, (float, int, type(None))
        # self.assertIsInstance(self.job.plane_wave_cutoff, (float, int))
        self.job.plane_wave_cutoff = 350
        self.assertEqual(self.job.input.incar["ENCUT"], 350)
        self.assertEqual(self.job.plane_wave_cutoff, 350)
        self.assertEqual(self.job.plane_wave_cutoff, self.job.encut)
        self.job.encut = 450
        self.assertEqual(self.job.encut, 450)
        self.assertEqual(self.job.input.incar["ENCUT"], 450)
        self.assertEqual(self.job.plane_wave_cutoff, 450)

    def test_exchange_correlation_functional(self):
        self.assertEqual(self.job.exchange_correlation_functional, "GGA")
        self.assertEqual(self.job.input.potcar["xc"], "GGA")
        self.job.exchange_correlation_functional = "LDA"
        self.assertEqual(self.job.exchange_correlation_functional, "LDA")
        self.assertEqual(self.job.input.potcar["xc"], "LDA")

    def test_get_nelect(self):
        atoms = CrystalStructure("Pt", BravaisBasis="fcc", a=3.98)
        self.job.structure = atoms
        self.assertEqual(self.job.get_nelect(), 10)

    def test_set_empty_states(self):
        atoms = CrystalStructure("Pt", BravaisBasis="fcc", a=3.98)
        self.job.structure = atoms
        self.job.set_empty_states(n_empty_states=10)
        self.assertEqual(self.job.input.incar["NBANDS"], 15)
        self.job.structure = atoms.repeat([3, 1, 1])
        self.job.set_empty_states(n_empty_states=10)
        self.assertEqual(self.job.input.incar["NBANDS"], 25)

    def test_calc_static(self):
        self.job.calc_static(
            electronic_steps=90,
            retain_charge_density=True,
            retain_electrostatic_potential=True,
        )
        self.assertEqual(self.job.input.incar["IBRION"], -1)
        self.assertEqual(self.job.input.incar["NELM"], 90)
        self.assertEqual(self.job.input.incar["LVTOT"], True)
        self.assertEqual(self.job.input.incar["LCHARG"], True)

    def test_set_structure(self):
        self.assertEqual(self.job.structure, None)
        atoms = CrystalStructure("Pt", BravaisBasis="fcc", a=3.98)
        self.job.structure = atoms
        self.assertEqual(self.job.structure, atoms)
        self.job.structure = None
        self.assertEqual(self.job.structure, None)
        self.job.structure = atoms
        self.assertEqual(self.job.structure, atoms)

    def test_list_potenitals(self):
        self.assertRaises(ValueError, self.job.list_potentials)

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


if __name__ == "__main__":
    unittest.main()
