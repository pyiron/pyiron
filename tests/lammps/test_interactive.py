# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import numpy as np
import os
from pyiron_base import Project, ProjectHDFio
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.lammps.lammps import Lammps


class InteractiveLibrary(object):
    def __init__(self):
        self._command = []

    def command(self, command_in):
        self._command.append(command_in)

    def scatter_atoms(self, *args):
        self._command.append(" ".join([str(arg) for arg in args]))


class TestLammpsInteractive(unittest.TestCase):
    def setUp(self):
        self.job._interactive_library = InteractiveLibrary()
        self.minimize_job._interactive_library = InteractiveLibrary()
        self.minimize_control_job._interactive_library = InteractiveLibrary()

    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, "lammps"))

        structure = Atoms(
            symbols="Fe2",
            positions=np.outer(np.arange(2), np.ones(3)),
            cell=2 * np.eye(3),
        )

        cls.job = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="lammps",
        )
        cls.job.server.run_mode.interactive = True
        cls.job.structure = structure

        cls.minimize_job = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="minimize_lammps",
        )
        cls.minimize_control_job = Lammps(
            project=ProjectHDFio(project=cls.project, file_name="lammps"),
            job_name="minimize_control_lammps",
        )
        # cls.control_job.server.run_mode.interactive = True  # Fails if we then call, e.g. `calc_minimize`

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, "lammps"))
        project.remove_jobs_silently(recursive=True)
        project.remove(enable=True)

    def test_interactive_cells_setter(self):
        self.job.interactive_cells_setter(np.eye(3))
        self.assertEqual(
            self.job._interactive_library._command[-1],
            "change_box all x final 0 1.000000 y final 0 1.000000 z final 0 1.000000 units box",
        )

    def test_interactive_positions_setter(self):
        self.job.interactive_positions_setter(np.arange(6).reshape(2, 3))
        self.assertTrue(self.job._interactive_library._command[0].startswith("x 1 3"))
        self.assertEqual(
            self.job._interactive_library._command[1], "change_box all remap"
        )

    def test_interactive_execute(self):
        self.job._interactive_lammps_input()
        self.assertEqual(
            self.job._interactive_library._command,
            [
                "fix ensemble all nve",
                "variable dumptime equal 100",
                "variable thermotime equal 100",
                "thermo_style custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol",
                "thermo_modify format float %20.15g",
                "thermo ${thermotime}",
            ],
        )

    def test_calc_minimize_input(self):
        # Ensure defaults match control
        atoms = Atoms("Fe8", positions=np.zeros((8, 3)), cell=np.eye(3))
        self.minimize_control_job.structure = atoms
        self.minimize_control_job.input.control.calc_minimize()
        self.minimize_control_job._interactive_lammps_input()
        self.minimize_job.structure = atoms
        self.minimize_job.calc_minimize()
        self.minimize_job._interactive_lammps_input()

        self.assertEqual(
            self.minimize_control_job._interactive_library._command,
            self.minimize_job._interactive_library._command
        )

        # Ensure that pressure inputs are being parsed OK
        self.minimize_job.calc_minimize(pressure=0)
        self.minimize_job._interactive_lammps_input()
        self.assertTrue(("fix ensemble all box/relax x 0.0 y 0.0 z 0.0 couple none" in
                         self.minimize_job._interactive_library._command))

        self.minimize_job.calc_minimize(pressure=[1, 2, None, 0., 0., None])
        self.minimize_job._interactive_lammps_input()
        self.assertTrue(("fix ensemble all box/relax x 10000.0 y 20000.0 xy 0.0 xz 0.0 couple none" in
                         self.minimize_job._interactive_library._command))


if __name__ == "__main__":
    unittest.main()
