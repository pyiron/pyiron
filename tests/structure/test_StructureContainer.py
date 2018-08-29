import unittest
import os
from pyiron_base.core.settings.generic import Settings
s = Settings(config={'sql_file': 'container.db',
                     'project_paths': os.path.abspath(os.getcwd()),
                     'resource_paths': os.path.join(os.path.abspath(os.getcwd()), '../static')})

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
        project.remove(enable=True)
        s.close_connection()
        os.remove('container.db')

    def test_container(self):
        structure_container = self.project.load(1)
        self.assertEqual(structure_container.job_id, 1)
        self.assertEqual(structure_container.job_name, 'structure_container')
        self.assertEqual(structure_container.project_hdf5.project_path, 'structure_testing/')
        self.assertTrue(structure_container.status.finished)
        self.assertEqual(structure_container.structure, self.basis)

if __name__ == '__main__':
    unittest.main()
