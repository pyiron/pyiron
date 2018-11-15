import unittest
import os
from pyiron.project import Project


class TestStructureContainer(unittest.TestCase):
    def setUp(self):
        self.lattice_constant = 3.5
        self.file_location = os.path.dirname(os.path.abspath(__file__))
        self.project = Project(os.path.join(self.file_location, 'structure_testing'))
        self.basis = self.project.create_structure(element="Fe", bravais_basis='fcc',
                                                   lattice_constant=self.lattice_constant)
        self.structure_container = self.project.create_job(self.project.job_type.StructureContainer,
                                                           "structure_container")
        self.structure_container.structure = self.basis

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, 'structure_testing'))
        ham = project.load(project.get_job_ids()[0])
        ham.remove()
        project.remove(enable=True)

    def test_container(self):
        structure_container = self.project.load(self.project.get_job_ids()[0])
        self.assertEqual(structure_container.job_id, self.project.get_job_ids()[0])
        self.assertEqual(structure_container.job_name, 'structure_container')
        self.assertTrue('tests/structure/structure_testing/' in structure_container.project_hdf5.project_path)
        self.assertTrue(structure_container.status.finished)
        self.assertEqual(structure_container.structure, self.basis)

if __name__ == '__main__':
    unittest.main()
