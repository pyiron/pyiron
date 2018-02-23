import os
from pyironbase.core.settings.config.testing import ConfigTesting
from pyironbase.core.settings.generic import Settings
import unittest

config = ConfigTesting(sql_lite_database='./structure_testing.db', path_project=str(os.getcwd()),
                       path_potentials='../../../static/potentials/')
s = Settings(config=config)

from pyiron.project import Project


class TestStructureContainer(unittest.TestCase):
    def setUp(self):
        self.lattice_constant = 3.5
        self.project = Project('structure_testing')
        self.basis = self.project.create_structure(element="Fe", bravais_basis='fcc',
                                                   lattice_constant=self.lattice_constant)
        self.structure_container = self.project.create_job("StructureContainer", "structure_container")
        self.structure_container.structure = self.basis

    @classmethod
    def tearDownClass(cls):
        project = Project('structure_testing')
        ham = project.load(1)
        ham.remove()
        project.remove()
        s.close_connection()
        os.remove('structure_testing.db')

    def test_container(self):
        structure_container = self.project.load(1)
        self.assertEqual(structure_container.job_id, 1)
        self.assertEqual(structure_container.job_name, 'structure_container')
        self.assertEqual(structure_container.project_hdf5.project_path, 'structure_testing/')
        self.assertTrue(structure_container.status.finished)
        self.assertEqual(structure_container.structure, self.basis)

if __name__ == '__main__':
    unittest.main()
