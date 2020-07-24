# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
import os
from pyiron.base.project.generic import Project
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.lammps.lammps import Lammps
from pyiron.lammps.base import LammpsStructure, UnfoldingPrism
import ase.units as units


class TestLammps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, "lammps"))
        cls.job = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="lammps",
        )
        cls.job_water = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps_water"),
            job_name="lammps_water",
        )
        cls.job_water_dump = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps_water_dump"),
            job_name="lammps_water_dump",
        )
        cls.job_dump = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps_dump_static"),
            job_name="lammps_dump_static",
        )
        cls.job_vcsgc_input = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps_vcsgc_input"),
            job_name="lammps_vcsgc_input",
        )
        cls.minimize_job = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="minimize_lammps",
        )
        cls.minimize_control_job = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="minimize_control_lammps",
        )
        cls.job_read_restart = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="read_restart",
        )
        cls.job_average = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="average",
        )
        cls.job_fail = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="fail",
        )

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, "lammps"))
        project.remove_jobs(recursive=True)
        project.remove(enable=True)

    def test_selective_dynamics(self):
        atoms = Atoms("Fe8", positions=np.zeros((8, 3)), cell=np.eye(3))
        atoms.add_tag(selective_dynamics=[True, True, True])
        self.job.structure = atoms
        self.job._set_selective_dynamics()
        self.assertFalse("group" in self.job.input.control._dataset["Parameter"])
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
        self.assertTrue(
            "group___constraintx" in self.job.input.control._dataset["Parameter"]
        )
        self.assertTrue(
            "group___constrainty" in self.job.input.control._dataset["Parameter"]
        )
        self.assertTrue(
            "group___constraintz" in self.job.input.control._dataset["Parameter"]
        )
        self.assertTrue(
            "group___constraintxy" in self.job.input.control._dataset["Parameter"]
        )
        self.assertTrue(
            "group___constraintyz" in self.job.input.control._dataset["Parameter"]
        )
        self.assertTrue(
            "group___constraintxz" in self.job.input.control._dataset["Parameter"]
        )
        self.assertTrue(
            "group___constraintxyz" in self.job.input.control._dataset["Parameter"]
        )

    def test_structure_atomic(self):
        atoms = Atoms("Fe1", positions=np.zeros((1, 3)), cell=np.eye(3))
        lmp_structure = LammpsStructure()
        lmp_structure._el_eam_lst = ["Fe"]
        lmp_structure.structure = atoms
        self.assertEqual(
            lmp_structure._dataset["Value"],
            [
                "Start File for LAMMPS",
                "1 atoms",
                "1 atom types",
                "",
                "0. 1.000000000000000 xlo xhi",
                "0. 1.000000000000000 ylo yhi",
                "0. 1.000000000000000 zlo zhi",
                "",
                "Masses",
                "",
                "1 55.845000",
                "",
                "Atoms",
                "",
                "1 1 0.000000000000000 0.000000000000000 0.000000000000000",
                "",
            ],
        )

    def test_structure_charge(self):
        atoms = Atoms("Fe1", positions=np.zeros((1, 3)), cell=np.eye(3))
        atoms.add_tag(charge=2.0)
        lmp_structure = LammpsStructure()
        lmp_structure.atom_type = "charge"
        lmp_structure._el_eam_lst = ["Fe"]
        lmp_structure.structure = atoms
        self.assertEqual(
            lmp_structure._dataset["Value"],
            [
                "Start File for LAMMPS",
                "1 atoms",
                "1 atom types",
                "",
                "0. 1.000000000000000 xlo xhi",
                "0. 1.000000000000000 ylo yhi",
                "0. 1.000000000000000 zlo zhi",
                "",
                "Masses",
                "",
                "1 55.845000",
                "",
                "Atoms",
                "",
                "1 1 2.000000 0.000000000000000 0.000000000000000 0.000000000000000",
                "",
            ],
        )

    def test_avilable_versions(self):
        self.job.executable = os.path.abspath(
            os.path.join(
                self.execution_path,
                "..",
                "static",
                "lammps",
                "bin",
                "run_lammps_2018.03.16.sh",
            )
        )
        self.assertTrue([2018, 3, 16] == self.job._get_executable_version_number())
        self.job.executable = os.path.abspath(
            os.path.join(
                self.execution_path,
                "..",
                "static",
                "lammps",
                "bin",
                "run_lammps_2018.03.16_mpi.sh",
            )
        )
        self.assertTrue([2018, 3, 16] == self.job._get_executable_version_number())

    def test_lammps_water(self):
        density = 1.0e-24  # g/A^3
        n_mols = 27
        mol_mass_water = 18.015  # g/mol
        # Determining the supercell size size
        mass = mol_mass_water * n_mols / units.mol  # g
        vol_h2o = mass / density  # in A^3
        a = vol_h2o ** (1.0 / 3.0)  # A
        # Constructing the unitcell
        n = int(round(n_mols ** (1.0 / 3.0)))
        dx = 0.7
        r_O = [0, 0, 0]
        r_H1 = [dx, dx, 0]
        r_H2 = [-dx, dx, 0]
        unit_cell = (a / n) * np.eye(3)
        water = Atoms(
            elements=["H", "H", "O"], positions=[r_H1, r_H2, r_O], cell=unit_cell, pbc=True)
        water.set_repeat([n, n, n])
        self.job_water.structure = water
        with self.assertWarns(UserWarning):
            self.job_water.potential = "H2O_tip3p"
        with self.assertRaises(ValueError):
            self.job_water.calc_md(temperature=[0, 100])
        with self.assertRaises(ValueError):
            self.job_water.calc_md(pressure=0)
        with self.assertRaises(ValueError):
            self.job_water.calc_md(temperature=[0, 100, 200])
        self.job_water.calc_md(
            temperature=350,
            initial_temperature=350,
            time_step=1,
            n_ionic_steps=1000,
            n_print=200,
        )
        file_directory = os.path.join(
            self.execution_path, "..", "static", "lammps_test_files"
        )
        self.job_water.restart_file_list.append(
            os.path.join(file_directory, "dump.out")
        )
        self.job_water.restart_file_list.append(
            os.path.join(file_directory, "log.lammps")
        )
        self.job_water.run(run_mode="manual")
        self.job_water.status.collect = True
        self.job_water.run()
        nodes = [
            "positions",
            "temperature",
            "energy_tot",
            "energy_pot",
            "steps",
            "time",
            "positions",
            "forces",
            "cells",
            "pressures",
            "unwrapped_positions",
        ]
        with self.job_water.project_hdf5.open("output/generic") as h_gen:
            hdf_nodes = h_gen.list_nodes()
            self.assertTrue(all([node in hdf_nodes for node in nodes]))
        self.assertTrue(
            np.array_equal(self.job_water["output/generic/positions"].shape, (6, 81, 3))
        )
        self.assertTrue(
            np.array_equal(
                self.job_water["output/generic/positions"].shape,
                self.job_water["output/generic/forces"].shape,
            )
        )
        self.assertEqual(len(self.job_water["output/generic/steps"]), 6)

    def test_dump_parser_water(self):
        density = 1.0e-24  # g/A^3
        n_mols = 27
        mol_mass_water = 18.015  # g/mol
        # Determining the supercell size size
        mass = mol_mass_water * n_mols / units.mol  # g
        vol_h2o = mass / density  # in A^3
        a = vol_h2o ** (1.0 / 3.0)  # A
        # Constructing the unitcell
        n = int(round(n_mols ** (1.0 / 3.0)))
        dx = 0.7
        r_O = [0, 0, 0]
        r_H1 = [dx, dx, 0]
        r_H2 = [-dx, dx, 0]
        unit_cell = (a / n) * np.eye(3)
        unit_cell[0][1] += 0.01
        water = Atoms(
            elements=["H", "H", "O"], positions=[r_H1, r_H2, r_O], cell=unit_cell, pbc=True)
        water.set_repeat([n, n, n])
        self.job_water_dump.structure = water
        with self.assertWarns(UserWarning):
            self.job_water_dump.potential = "H2O_tip3p"
        self.job_water_dump.calc_md(
            temperature=350,
            initial_temperature=350,
            time_step=1,
            n_ionic_steps=1000,
            n_print=200,
            pressure=0,
        )
        self.assertFalse('nan' in self.job_water_dump.input.control['fix___ensemble'])
        file_directory = os.path.join(
            self.execution_path, "..", "static", "lammps_test_files"
        )
        self.job_water_dump.restart_file_list.append(
            os.path.join(file_directory, "log.lammps")
        )
        self.job_water_dump.restart_file_list.append(
            os.path.join(file_directory, "dump.out")
        )
        self.job_water_dump.run(run_mode="manual")
        self.job_water_dump.status.collect = True
        self.job_water_dump.run()
        positions = np.loadtxt(os.path.join(file_directory, "positions_water.dat"))
        positions = positions.reshape(len(positions), -1, 3)
        forces = np.loadtxt(os.path.join(file_directory, "forces_water.dat"))
        forces = forces.reshape(len(forces), -1, 3)
        self.assertTrue(
            np.allclose(
                self.job_water_dump["output/generic/unwrapped_positions"], positions
            )
        )
        self.assertTrue(
            np.allclose(self.job_water_dump["output/generic/forces"], forces)
        )
        self.job_water_dump.write_traj(filename="test.xyz",
                                       file_format="xyz")
        atom_indices = self.job_water_dump.structure.select_index("H")
        snap_indices = [1, 3, 4]
        orig_pos = self.job_water_dump.output.positions
        self.job_water_dump.write_traj(filename="test.xyz",
                                       file_format="xyz",
                                       atom_indices=atom_indices,
                                       snapshot_indices=snap_indices)
        self.job_water_dump.write_traj(filename="test.xyz",
                                       file_format="xyz",
                                       atom_indices=atom_indices,
                                       snapshot_indices=snap_indices,
                                       overwrite_positions=np.zeros_like(orig_pos))
        self.assertRaises(ValueError, self.job_water_dump.write_traj, filename="test.xyz",
                          file_format="xyz",
                          atom_indices=atom_indices,
                          snapshot_indices=snap_indices,
                          overwrite_positions=np.zeros_like(orig_pos)[:-1])

        self.job_water_dump.write_traj(filename="test.xyz",
                                       file_format="xyz",
                                       atom_indices=atom_indices,
                                       snapshot_indices=snap_indices,
                                       overwrite_positions=np.zeros_like(orig_pos),
                                       overwrite_cells=self.job_water_dump.trajectory()._cells)
        self.job_water_dump.write_traj(filename="test.xyz",
                                       file_format="xyz",
                                       atom_indices=atom_indices,
                                       snapshot_indices=snap_indices,
                                       overwrite_positions=np.zeros_like(orig_pos)[:-1],
                                       overwrite_cells=self.job_water_dump.trajectory()._cells[:-1])
        self.assertRaises(ValueError, self.job_water_dump.write_traj, filename="test.xyz",
                          file_format="xyz",
                          atom_indices=atom_indices,
                          snapshot_indices=snap_indices,
                          overwrite_positions=np.zeros_like(orig_pos),
                          overwrite_cells=self.job_water_dump.trajectory()._cells[:-1])
        os.remove("test.xyz")
        self.assertTrue(np.array_equal(self.job_water_dump.trajectory()._positions,
                                       orig_pos))
        self.assertTrue(np.array_equal(self.job_water_dump.trajectory(stride=2)._positions,
                                       orig_pos[::2]))
        self.assertTrue(np.array_equal(
            self.job_water_dump.trajectory(atom_indices=atom_indices,
                                           snapshot_indices=snap_indices)._positions,
            orig_pos[snap_indices][:, atom_indices, :]))

        nx, ny, nz = orig_pos.shape
        random_array = np.random.rand(nx, ny, nz)
        random_cell = np.random.rand(nx, 3, 3)
        self.assertTrue(np.array_equal(
            self.job_water_dump.trajectory(atom_indices=atom_indices,
                                           snapshot_indices=snap_indices,
                                           overwrite_positions=random_array)._positions,
            random_array[snap_indices][:, atom_indices, :]))
        self.assertTrue(np.array_equal(
            self.job_water_dump.trajectory(atom_indices=atom_indices,
                                           snapshot_indices=snap_indices,
                                           overwrite_positions=random_array,
                                           overwrite_cells=random_cell)._cells,
            random_cell[snap_indices]))
        self.assertIsInstance(self.job_water_dump.get_structure(-1), Atoms)
        # Test for clusters
        with self.job_water_dump.project_hdf5.open("output/generic") as h_out:
            h_out["cells"] = None
        self.assertTrue(np.array_equal(
            self.job_water_dump.trajectory(atom_indices=atom_indices,
                                           snapshot_indices=snap_indices)._positions,
            orig_pos[snap_indices][:, atom_indices, :]))
        with self.job_water_dump.project_hdf5.open("output/generic") as h_out:
            h_out["cells"] = np.array([np.zeros((3, 3))] * len(h_out["positions"]))
        self.assertTrue(np.array_equal(
            self.job_water_dump.trajectory(atom_indices=atom_indices,
                                           snapshot_indices=snap_indices)._positions,
            orig_pos[snap_indices][:, atom_indices, :]))

    def test_dump_parser(self):
        structure = Atoms(
            elements=2 * ["Fe"],
            cell=2.78 * np.eye(3),
            positions=2.78 * np.outer(np.arange(2), np.ones(3)) * 0.5,
        )
        self.job_dump.structure = structure
        self.job_dump.potential = self.job_dump.list_potentials()[0]
        file_directory = os.path.join(
            self.execution_path, "..", "static", "lammps_test_files"
        )
        self.job_dump.collect_dump_file(cwd=file_directory, file_name="dump_static.out")
        self.assertTrue(
            np.array_equal(self.job_dump["output/generic/forces"].shape, (1, 2, 3))
        )
        self.assertTrue(
            np.array_equal(self.job_dump["output/generic/positions"].shape, (1, 2, 3))
        )
        self.assertTrue(
            np.array_equal(self.job_dump["output/generic/cells"].shape, (1, 3, 3))
        )
        self.assertTrue(
            np.array_equal(self.job_dump["output/generic/indices"].shape, (1, 2))
        )

    def test_vcsgc_input(self):
        unit_cell = Atoms(
            elements=['Al', 'Al', 'Al', 'Mg'],
            positions=[
                [0., 0., 0.],
                [0., 2., 2.],
                [2., 0., 2.],
                [2., 2., 0.]
            ],
            cell=4 * np.eye(3)
        )
        self.job_vcsgc_input.structure = unit_cell
        self.job_vcsgc_input.potential = self.job_vcsgc_input.list_potentials()[0]
        symbols = self.job_vcsgc_input.input.potential.get_element_lst()

        bad_element = {s: 0. for s in symbols}
        bad_element.update({'X': 1.})  # Non-existant chemical symbol
        self.assertRaises(
            ValueError, self.job_vcsgc_input.calc_vcsgc, mu=bad_element, temperature_mc=300.
        )

        self.assertRaises(
            ValueError, self.job_vcsgc_input.calc_vcsgc, target_concentration=bad_element, temperature_mc=300.
        )

        bad_conc = {s: 0. for s in symbols}
        bad_conc['Al'] = 0.99
        self.assertRaises(
            ValueError, self.job_vcsgc_input.calc_vcsgc, target_concentration=bad_conc, temperature_mc=300.
        )

        self.assertRaises(
            ValueError, self.job_vcsgc_input.calc_vcsgc, window_moves=-1, temperature_mc=300.
        )
        self.assertRaises(
            ValueError, self.job_vcsgc_input.calc_vcsgc, window_moves=1.1, temperature_mc=300.
        )

        self.assertRaises(
            ValueError, self.job_vcsgc_input.calc_vcsgc, window_size=0.3, temperature_mc=300.
        )

        mu = {s: 0. for s in symbols}
        mu[symbols[0]] = 1.
        self.assertRaises(
            ValueError, self.job_vcsgc_input.calc_vcsgc, mu=mu, temperature_mc=None, temperature=None
        )


        args = dict(
            mu=mu,
            target_concentration=None,
            kappa=1000.0,
            mc_step_interval=100,
            swap_fraction=0.1,
            temperature_mc=None,
            window_size=None,
            window_moves=None,
            seed=1,
            temperature=300.0,
        )
        input_string = 'all sgcmc {0} {1} {2} {3} randseed {4}'.format(
            args['mc_step_interval'],
            args['swap_fraction'],
            args['temperature'],
            ' '.join([str(args['mu'][symbol] - args['mu'][symbols[0]]) for symbol in symbols[1:]]),
            args['seed']
        )
        self.job_vcsgc_input.calc_vcsgc(**args)
        self.assertEqual(self.job_vcsgc_input.input.control['fix___vcsgc'], input_string)

        args['temperature_mc'] = 100.,
        input_string = 'all sgcmc {0} {1} {2} {3} randseed {4}'.format(
            args['mc_step_interval'],
            args['swap_fraction'],
            args['temperature_mc'],
            ' '.join([str(args['mu'][symbol] - args['mu'][symbols[0]]) for symbol in symbols[1:]]),
            args['seed']
        )
        self.job_vcsgc_input.calc_vcsgc(**args)
        self.assertEqual(self.job_vcsgc_input.input.control['fix___vcsgc'], input_string)

        conc = {s: 0. for s in symbols}
        conc[symbols[0]] = 0.5
        conc[symbols[-1]] = 0.5
        args['target_concentration'] = conc
        input_string += ' variance {0} {1}'.format(
            args['kappa'],
            ' '.join([str(conc[symbol]) for symbol in symbols[1:]])
        )
        self.job_vcsgc_input.calc_vcsgc(**args)
        self.assertEqual(self.job_vcsgc_input.input.control['fix___vcsgc'], input_string)

        args['window_moves'] = 10
        input_string += ' window_moves {0}'.format(args['window_moves'])
        self.job_vcsgc_input.calc_vcsgc(**args)
        self.assertEqual(self.job_vcsgc_input.input.control['fix___vcsgc'], input_string)

        args['window_size'] = 0.75
        input_string += ' window_size {0}'.format(args['window_size'])
        self.job_vcsgc_input.calc_vcsgc(**args)
        self.assertEqual(self.job_vcsgc_input.input.control['fix___vcsgc'], input_string)

    def test_calc_minimize_input(self):
        # Ensure that defaults match control defaults
        atoms = Atoms("Fe8", positions=np.zeros((8, 3)), cell=np.eye(3))
        self.minimize_control_job.structure = atoms
        self.minimize_control_job.input.control.calc_minimize()

        self.minimize_job.sturcture = atoms
        self.minimize_job._prism = UnfoldingPrism(atoms.cell)
        self.minimize_job.calc_minimize()
        for k in self.job.input.control.keys():
            self.assertEqual(self.minimize_job.input.control[k], self.minimize_control_job.input.control[k])

        # Ensure that pressure inputs are being parsed OK
        self.minimize_control_job.calc_minimize(pressure=0)
        self.assertEqual(
            self.minimize_control_job.input.control['fix___ensemble'],
            "all box/relax x 0.0 y 0.0 z 0.0 couple none"
        )

        self.minimize_control_job.calc_minimize(pressure=[1, 2, None, 0., 0., None])
        self.assertEqual(
            self.minimize_control_job.input.control['fix___ensemble'],
            "all box/relax x 10000.0 y 20000.0 xy 0.0 xz 0.0 couple none"
        )

    def test_read_restart_file(self):
        self.job_read_restart.read_restart_file()
        self.assertIsNone(self.job_read_restart['dimension'])

    def test_write_restart(self):
        self.job_read_restart.write_restart_file()
        self.assertEqual(self.job_read_restart.input.control['write_restart'], 'restart.out')

    def test_average(self):
        a_0 = 2.855312531
        atoms = Atoms("Fe2", positions=[3*[0], 3*[0.5*a_0]], cell=a_0*np.eye(3))
        self.job_average.structure = atoms
        self.job_average.potential = 'Fe_C_Becquart_eam'
        file_directory = os.path.join(
            self.execution_path, "..", "static", "lammps_test_files"
        )
        self.job_average.collect_dump_file(cwd=file_directory, file_name="dump_average.out")
        self.job_average.collect_output_log(cwd=file_directory, file_name="log_average.lammps")

    def test_validate(self):
        with self.assertRaises(ValueError):
            self.job_fail.validate_ready_to_run()
        a_0 = 2.855312531
        atoms = Atoms("Fe2", positions=[3 * [0], 3 * [0.5 * a_0]], cell=a_0 * np.eye(3), pbc=False)
        self.job_fail.structure = atoms
        with self.assertRaises(ValueError):
            self.job_fail.validate_ready_to_run()
        self.job_fail.potential = self.job_fail.list_potentials()[-1]
        self.job_fail.validate_ready_to_run()
        self.job_fail.structure.positions[0, 0] -= 2.855
        with self.assertRaises(ValueError):
            self.job_fail.validate_ready_to_run()
        self.job_fail.structure.pbc = True
        self.job_fail.validate_ready_to_run()
        self.job_fail.structure.pbc = [True, True, False]
        self.job_fail.validate_ready_to_run()
        self.job_fail.structure.pbc = [False, True, True]
        with self.assertRaises(ValueError):
            self.job_fail.validate_ready_to_run()


if __name__ == "__main__":
    unittest.main()
