import unittest
import numpy as np
import os
import posixpath
from pyiron.base.project.generic import Project
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.lammps.lammps import Lammps
from pyiron.lammps.base import LammpsStructure
import ase.units as units
import shutil


class TestLammps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, 'lammps'))
        cls.job = Lammps(project=ProjectHDFio(project=cls.project, file_name='lammps'), job_name='lammps')
        cls.job_water = Lammps(project=ProjectHDFio(project=cls.project, file_name='lammps_water'),
                               job_name='lammps_water')

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

    def test_structure_atomic(self):
        atoms = Atoms('Fe1', positions=np.zeros((1, 3)), cell=np.eye(3))
        lmp_structure = LammpsStructure()
        lmp_structure._el_eam_lst = ['Fe']
        lmp_structure.structure = atoms
        self.assertEqual(lmp_structure._dataset['Value'], ['File for LAMMPS',
                                                           'atoms',
                                                           'atom types',
                                                           '',
                                                           '1.000000000000000 xlo xhi',
                                                           '1.000000000000000 ylo yhi',
                                                           '1.000000000000000 zlo zhi',
                                                           '',
                                                           '',
                                                           '',
                                                           '55.845001',
                                                           '',
                                                           '',
                                                           '',
                                                           '1 0.000000000000000 0.000000000000000 0.000000000000000',
                                                           ''])

    def test_structure_charge(self):
        atoms = Atoms('Fe1', positions=np.zeros((1, 3)), cell=np.eye(3))
        atoms.add_tag(charge=2.0)
        lmp_structure = LammpsStructure()
        lmp_structure.atom_type = 'charge'
        lmp_structure._el_eam_lst = ['Fe']
        lmp_structure.structure = atoms
        self.assertEqual(lmp_structure._dataset['Value'], ['File for LAMMPS',
                                                           'atoms',
                                                           'atom types',
                                                           '',
                                                           '1.000000000000000 xlo xhi',
                                                           '1.000000000000000 ylo yhi',
                                                           '1.000000000000000 zlo zhi',
                                                           '',
                                                           '',
                                                           '',
                                                           '55.845001',
                                                           '',
                                                           '',
                                                           '',
                                                           '1 2.000000 0.000000000000000 0.000000000000000 0.000000000000000',
                                                           ''])

    def test_avilable_versions(self):
        self.job.executable = os.path.join(self.execution_path, '../static/lammps/bin/run_lammps_2018.03.16.sh')
        self.assertTrue([2018, 3, 16] == self.job._get_executable_version_number())
        self.job.executable = os.path.join(self.execution_path, '../static/lammps/bin/run_lammps_2018.03.16_mpi.sh')
        self.assertTrue([2018, 3, 16] == self.job._get_executable_version_number())

    def test_lammps_water(self):
        density = 1.0e-24  # g/A^3
        n_mols = 27
        mol_mass_water = 18.015  # g/mol
        # Determining the supercell size size
        mass = mol_mass_water * n_mols / units.mol  # g
        vol_h2o = mass / density  # in A^3
        a = vol_h2o ** (1. / 3.)  # A
        # Constructing the unitcell
        n = int(round(n_mols ** (1. / 3.)))
        dx = 0.7
        r_O = [0, 0, 0]
        r_H1 = [dx, dx, 0]
        r_H2 = [-dx, dx, 0]
        unit_cell = (a / n) * np.eye(3)
        water = Atoms(elements=['H', 'H', 'O'], positions=[r_H1, r_H2, r_O], cell=unit_cell)
        water.set_repeat([n, n, n])
        self.job_water.structure = water
        self.job_water.potential = 'H2O_tip3p'
        self.job_water.calc_md(temperature=350, initial_temperature=350, time_step=1, n_ionic_steps=1000, n_print=200)
        file_directory = posixpath.join(self.execution_path, "../static/lammps_test_files")
        self.job_water.restart_file_list.append(posixpath.join(file_directory, "dump.out"))
        self.job_water.restart_file_list.append(posixpath.join(file_directory, "log.lammps"))
        self.job_water.run(run_mode="manual")
        self.job_water.status.collect = True
        self.job_water.run()
        nodes = ["positions", "temperature", "temperatures", "energy_tot", "steps", "positions", "forces", "cells",
                 "pressures"]
        with self.job_water.project_hdf5.open("output/generic") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))
        self.assertTrue(np.array_equal(self.job_water["output/generic/positions"].shape, (6, 81, 3)))
        self.assertTrue(np.array_equal(self.job_water["output/generic/positions"].shape,
                                       self.job_water["output/generic/forces"].shape))
        self.assertEqual(len(self.job_water["output/generic/steps"]), 6)
        self.assertTrue(np.array_equal(self.job_water["output/generic/temperatures"],
                                       self.job_water["output/generic/temperature"]))


if __name__ == '__main__':
    unittest.main()
