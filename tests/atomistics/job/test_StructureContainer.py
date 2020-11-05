# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron.project import Project


class TestStructureContainer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lattice_constant = 3.5
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "structure_testing"))
        cls.basis = cls.project.create_structure(
            element="Fe", bravais_basis="fcc", lattice_constant=cls.lattice_constant
        )
        cls.structure_container = cls.project.create_job(
            cls.project.job_type.StructureContainer, "structure_container"
        )
        cls.structure_container.structure = cls.basis
        cls.structure_container.run()

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, "structure_testing"))
        ham = project.load(project.get_job_ids()[0])
        ham.remove()
        project.remove(enable=True)

    def test_container(self):
        structure_container = self.project.load(self.project.get_job_ids()[0])
        self.assertEqual(structure_container.job_id, self.project.get_job_ids()[0])
        self.assertEqual(structure_container.job_name, "structure_container")
        self.assertTrue(
            "atomistics/job/structure_testing/"
            in structure_container.project_hdf5.project_path
        )
        self.assertTrue(structure_container.status.finished)
        self.assertEqual(structure_container.structure, self.basis)


if __name__ == "__main__":
    unittest.main()
