# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import numpy as np
import unittest
import warnings
from pyiron.project import Project
from pyiron.atomistics.structure.periodic_table import PeriodicTable
from pyiron.atomistics.structure.atoms import Atoms


class InteractiveLibrary(object):
    def __init__(self):
        self.command = []

    def write(self, line):
        self.command.append(line)

    def flush(self):
        return None


class TestSphinx(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "../static/sphinx"))
        cls.sphinx = cls.project.create_job("Sphinx", "job_sphinx")
        cls.sphinx.structure = Atoms(elements=['Fe']*2, scaled_positions=[3*[0.0], 3*[0.5]], cell=2.6*np.eye(3))
        cls.sphinx.structure.set_initial_magnetic_moments(np.ones(2))
        cls.current_dir = os.path.abspath(os.getcwd())
        cls.sphinx._create_working_directory()
        cls.sphinx.input["VaspPot"] = False
        cls.sphinx.load_default_groups()
        cls.sphinx.write_input()
        cls.sphinx.version = "2.6.1"
        cls.sphinx.server.run_mode.interactive = True

    def setUp(self):
        self.sphinx._interactive_library = InteractiveLibrary()

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/input.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/Fe_GGA.atomicdata",
            )
        )
        os.rmdir(
            os.path.join(
                cls.file_location, "../static/sphinx/job_sphinx_hdf5/job_sphinx"
            )
        )
        os.rmdir(os.path.join(cls.file_location, "../static/sphinx/job_sphinx_hdf5"))

    def test_interactive_cells_setter(self):
        with self.assertRaises(NotImplementedError):
            self.sphinx.interactive_cells_setter(np.eye(3))

    def test_interactive_pipe_write(self):
        with self.assertRaises(TypeError) as c:
            self.sphinx._interactive_pipe_write(self.sphinx)

    def test_interactive_positions_setter(self):
        self.sphinx.interactive_positions_setter(np.zeros((2,3)))
        self.assertEqual(self.sphinx._interactive_library.command,
                         ['set structure\n', '0.0\n', '0.0\n', '0.0\n', '0.0\n', '0.0\n', '0.0\n'])

    def test_interactive_spin_constraints_setter(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.sphinx.interactive_spin_constraints_setter(np.zeros(2))
            self.assertEqual(len(w), 1)
        self.sphinx.fix_spin_constraint = True
        self.sphinx.interactive_spin_constraints_setter(np.zeros(2))
        self.assertEqual(self.sphinx._interactive_library.command, ['set spinconstraint\n', '0.0\n', '0.0\n'])


if __name__ == "__main__":
    unittest.main()
