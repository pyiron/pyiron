# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron.base.project.generic import Project


class TestChildids(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "testing_childids"))

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, "testing_childids"))
        project.remove_jobs_silently(recursive=True)
        project.remove(enable=True)

    def test_childids(self):
        ham_master_1 = self.project.create_job("ScriptJob", "master")
        ham_master_1.save()
        ham_child_1_1 = self.project.create_job("ScriptJob", "child_1_1")
        ham_child_1_1.master_id = ham_master_1.job_id
        ham_child_1_1.save()
        self.assertEqual(ham_master_1.child_ids, [ham_child_1_1.job_id])
        ham_child_1_2 = self.project.create_job("ScriptJob", "child_1_2")
        ham_child_1_2.master_id = ham_master_1.job_id
        ham_child_1_2.save()
        self.assertEqual(
            ham_master_1.child_ids, [ham_child_1_1.job_id, ham_child_1_2.job_id]
        )
        ham_child_1_3 = self.project.create_job("ScriptJob", "child_1_3")
        ham_child_1_3.master_id = ham_master_1.job_id
        ham_child_1_3.save()
        self.assertEqual(
            ham_master_1.child_ids,
            [ham_child_1_1.job_id, ham_child_1_2.job_id, ham_child_1_3.job_id],
        )
        ham_child_1_4 = self.project.create_job("ScriptJob", "child_1_4")
        ham_child_1_4.master_id = ham_master_1.job_id
        ham_child_1_4.save()
        self.assertEqual(
            ham_master_1.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
            ],
        )
        ham_child_1_5 = self.project.create_job("ScriptJob", "child_1_5")
        ham_child_1_5.master_id = ham_master_1.job_id
        ham_child_1_5.save()
        self.assertEqual(
            ham_master_1.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
                ham_child_1_5.job_id,
            ],
        )
        ham_child_1_6 = self.project.create_job("ScriptJob", "child_1_6")
        ham_child_1_6.master_id = ham_master_1.job_id
        ham_child_1_6.save()
        self.assertEqual(
            ham_master_1.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
                ham_child_1_5.job_id,
                ham_child_1_6.job_id,
            ],
        )
        ham_child_1_7 = self.project.create_job("ScriptJob", "child_1_7")
        ham_child_1_7.master_id = ham_master_1.job_id
        ham_child_1_7.save()
        self.assertEqual(
            ham_master_1.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
                ham_child_1_5.job_id,
                ham_child_1_6.job_id,
                ham_child_1_7.job_id,
            ],
        )
        ham_child_1_8 = self.project.create_job("ScriptJob", "child_1_8")
        ham_child_1_8.master_id = ham_master_1.job_id
        ham_child_1_8.save()
        self.assertEqual(
            ham_master_1.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
                ham_child_1_5.job_id,
                ham_child_1_6.job_id,
                ham_child_1_7.job_id,
                ham_child_1_8.job_id,
            ],
        )
        ham_child_1_9 = self.project.create_job("ScriptJob", "child_1_9")
        ham_child_1_9.master_id = ham_master_1.job_id
        ham_child_1_9.save()
        self.assertEqual(
            ham_master_1.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
                ham_child_1_5.job_id,
                ham_child_1_6.job_id,
                ham_child_1_7.job_id,
                ham_child_1_8.job_id,
                ham_child_1_9.job_id,
            ],
        )
        ham_master_2 = self.project.create_job("ScriptJob", "master_2")
        ham_master_2.save()
        ham_child_2_1 = self.project.create_job("ScriptJob", "child_2_1")
        ham_child_2_1.master_id = ham_master_2.job_id
        ham_child_2_1.save()
        self.assertEqual(ham_master_2.child_ids, [ham_child_2_1.job_id])
        self.assertEqual(
            ham_master_1.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
                ham_child_1_5.job_id,
                ham_child_1_6.job_id,
                ham_child_1_7.job_id,
                ham_child_1_8.job_id,
                ham_child_1_9.job_id,
            ],
        )
        ham_master_1_reload = self.project.load(ham_master_1.job_id)
        ham_master_2_reload = self.project.load(ham_master_2.job_id)
        self.assertEqual(ham_master_2_reload.child_ids, [ham_child_2_1.job_id])
        self.assertEqual(
            ham_master_1_reload.child_ids,
            [
                ham_child_1_1.job_id,
                ham_child_1_2.job_id,
                ham_child_1_3.job_id,
                ham_child_1_4.job_id,
                ham_child_1_5.job_id,
                ham_child_1_6.job_id,
                ham_child_1_7.job_id,
                ham_child_1_8.job_id,
                ham_child_1_9.job_id,
            ],
        )
