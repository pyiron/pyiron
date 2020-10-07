# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron.project import Project
from pyiron.table.funct import get_majority


class TestDatamining(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.execution_path, "table"))
        cls.project.create_table()

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(cls.execution_path, "table"))
        project.remove_jobs_silently(recursive=True)
        project.remove(enable=True)

    def test_get_majority(self):
        lst = [1, 1, 2]
        majority = get_majority(lst=lst, minority=False)
        self.assertTrue(majority == 1)
        majority, minority = get_majority(lst=lst, minority=True)
        self.assertTrue(majority == 1)
        self.assertTrue(minority == [2])

    def test_vasp_import(self):
        def filter_job_type(job):
            return job.__name__ == "Vasp"

        def filter_job_type_db(job_table):
            return job_table.hamilton == "Vasp"

        self.project.import_from_path(
            path=os.path.join(
                self.execution_path, "../static/vasp_test_files/full_job_sample"
            ),
            recursive=False,
        )
        table = self.project.create_table(delete_existing_job=True)
        table.filter_function = filter_job_type
        add_funtions(table)
        table.run()
        df = table.get_dataframe()

        table_new = self.project.create_table("table_new")
        # only use database filtering not job based filtering
        table_new.db_filter_function = filter_job_type_db
        add_funtions(table_new)
        table_new.run()
        df_new = table_new.get_dataframe()
        self.assertTrue(df.equals(df_new))
        self.assertEqual(df["Number_of_atoms"].values[0], 2)
        self.assertEqual(df["Fe"].values[0], 2)
        self.assertTrue(df["Convergence"].values[0])
        self.assertEqual(df["Number_of_species"].values[0], 1)
        self.assertEqual(df["Number_of_ionic_steps"].values[0], 1)
        self.assertEqual(df["encut"].values[0], 250)
        self.assertEqual(df["n_kpts"].values[0], 4)
        self.assertEqual(df["majority_element"].values[0], "Fe")
        self.assertEqual(df["minority_element_list"].values[0], "[]")
        self.assertEqual(df["job_name"].values[0], "full_job_sample")
        self.assertEqual(df["energy_tot"].values[0], -17.7331698)
        self.assertEqual(df["energy_free"].values[0], -17.7379867884)
        self.assertEqual(df["energy_int"].values[0], -17.72353582)
        self.assertEqual(df["alat"].values[0], 0.0)
        self.assertEqual(df["magnetic_structure"].values[0], "ferro-magnetic")
        self.assertEqual(df["avg. plane waves"].values[0], 196.375)
        self.assertEqual(df["energy_tot_wo_kin_corr"].values[0], -17.6003698)
        self.assertEqual(df["volume"].values[0], 21.95199999999999)
        # print(self.project.job_table())
        table_loaded = self.project["table_new"]
        df_loaded = table_loaded.get_dataframe()
        self.assertTrue(df.equals(df_loaded))


def get_alat(job):
    return 2 * job["input/structure/cell/cell"][0, 1]


def add_funtions(table):
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
    _ = table.add.get_average_waves
    _ = table.add.get_ekin_error
    _ = table.add.get_volume
    table.add["alat"] = get_alat


if __name__ == "__main__":
    unittest.main()
