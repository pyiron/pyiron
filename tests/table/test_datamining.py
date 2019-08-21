# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron.project import Project
from pyiron.table.funct import _get_majority


class TestDatamining(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, 'table'))

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, 'table'))
        project.remove_jobs(recursive=True)
        project.remove(enable=True)

    def test_get_majority(self):
        lst = [1, 1, 2]
        majority = _get_majority(lst=lst, minority=False)
        self.assertTrue(majority == 1)
        majority, minority = _get_majority(lst=lst, minority=True)
        self.assertTrue(majority == 1)
        self.assertTrue(minority == [2])

    def test_vasp_import(self):
        def filter_job_type(job):
            return job.__name__ == 'Vasp'

        def get_alat(job):
            return 2 * job['input/structure/cell/cell'][0, 1]

        self.project.import_from_path(path=os.path.join(self.execution_path,
                                                        '../static/vasp_test_files/full_job_sample'),
                                      recursive=False)
        table = self.project.create_table()
        table.filter_function = filter_job_type
        _ = table.add.get_sigma
        _ = table.add.get_total_number_of_atoms
        _ = table.add.get_elements
        _ = table.add.get_convergence_check
        _ = table.add.get_number_of_species
        _ = table.add.get_number_of_ionic_steps
        _ = table.add.get_ismear
        _ = table.add.get_encut
        _ = table.add.get_n_kpts
        _ = table.add.get_number_of_final_electronic_steps
        _ = table.add.get_majority_species
        _ = table.add.get_job_name
        _ = table.add.get_job_id
        _ = table.add.get_energy_tot
        _ = table.add.get_energy_free
        _ = table.add.get_energy_int
        _ = table.add.get_equilibrium_parameters
        _ = table.add.get_magnetic_structure
        table.add['alat'] = get_alat
        table.run()
        df = table.get_dataframe()
        self.assertEqual(df['Number_of_atoms'].values[0], 2)
        self.assertEqual(df['Fe'].values[0], 2)
        self.assertTrue(df['Convergence'].values[0])
        self.assertEqual(df['Number_of_species'].values[0], 1)
        self.assertEqual(df['Number_of_ionic_steps'].values[0], 1)
        self.assertEqual(df['encut'].values[0], 250)
        self.assertEqual(df['n_kpts'].values[0], 4)
        self.assertEqual(df['majority_element'].values[0], 'Fe')
        self.assertEqual(df['minority_element_list'].values[0], '[]')
        self.assertEqual(df['job_name'].values[0], 'full_job_sample')
        self.assertEqual(df['energy_tot'].values[0], -17.7331698)
        self.assertEqual(df['energy_free'].values[0], -17.73798679)
        self.assertEqual(df['energy_int'].values[0], -17.72353582)
        self.assertEqual(df['alat'].values[0], 0.0)
        self.assertEqual(df['magnetic_structure'].values[0], 'ferro-magnetic')


if __name__ == '__main__':
    unittest.main()
