import os
from pyiron_base.core.settings.config.testing import ConfigTesting
from pyiron_base.core.settings.generic import Settings
import unittest

config = ConfigTesting(sql_lite_database='./testing_childids.db', path_project=str(os.getcwd()),
                       path_potentials='../../../static/potentials/')
s = Settings(config=config)

from pyiron_base.project import Project


class TestChildids(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Project('testing_childids')

    @classmethod
    def tearDownClass(cls):
        project = Project('testing_childids')
        project.remove_jobs(recursive=True)
        project.remove()
        s.close_connection()
        os.remove('testing_childids.db')

    def test_childids(self):
        ham_master_1 = self.project.create_job('ScriptJob', "master")
        ham_master_1.save()
        self.assertEqual(ham_master_1.job_id, 1)
        ham_child_1_1 = self.project.create_job('ScriptJob', "child_1_1")
        ham_child_1_1.master_id = ham_master_1.job_id
        ham_child_1_1.save()
        self.assertEqual(ham_child_1_1.job_id, 2)
        self.assertEqual(ham_master_1.child_ids, [2])
        ham_child_1_2 = self.project.create_job('ScriptJob', "child_1_2")
        ham_child_1_2.master_id = ham_master_1.job_id
        ham_child_1_2.save()
        self.assertEqual(ham_child_1_2.job_id, 3)
        self.assertEqual(ham_master_1.child_ids, [2, 3])
        ham_child_1_3 = self.project.create_job('ScriptJob', "child_1_3")
        ham_child_1_3.master_id = ham_master_1.job_id
        ham_child_1_3.save()
        self.assertEqual(ham_child_1_3.job_id, 4)
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4])
        ham_child_1_4 = self.project.create_job('ScriptJob', "child_1_4")
        ham_child_1_4.master_id = ham_master_1.job_id
        ham_child_1_4.save()
        self.assertEqual(ham_child_1_4.job_id, 5)
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4, 5])
        ham_child_1_5 = self.project.create_job('ScriptJob', "child_1_5")
        ham_child_1_5.master_id = ham_master_1.job_id
        ham_child_1_5.save()
        self.assertEqual(ham_child_1_5.job_id, 6)
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4, 5, 6])
        ham_child_1_6 = self.project.create_job('ScriptJob', "child_1_6")
        ham_child_1_6.master_id = ham_master_1.job_id
        ham_child_1_6.save()
        self.assertEqual(ham_child_1_6.job_id, 7)
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4, 5, 6, 7])
        ham_child_1_7 = self.project.create_job('ScriptJob', "child_1_7")
        ham_child_1_7.master_id = ham_master_1.job_id
        ham_child_1_7.save()
        self.assertEqual(ham_child_1_7.job_id, 8)
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4, 5, 6, 7, 8])
        ham_child_1_8 = self.project.create_job('ScriptJob', "child_1_8")
        ham_child_1_8.master_id = ham_master_1.job_id
        ham_child_1_8.save()
        self.assertEqual(ham_child_1_8.job_id, 9)
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4, 5, 6, 7, 8, 9])
        ham_child_1_9 = self.project.create_job('ScriptJob', "child_1_9")
        ham_child_1_9.master_id = ham_master_1.job_id
        ham_child_1_9.save()
        self.assertEqual(ham_child_1_9.job_id, 10)
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4, 5, 6, 7, 8, 9, 10])
        ham_master_2 = self.project.create_job('ScriptJob', "master_2")
        ham_master_2.save()
        self.assertEqual(ham_master_2.job_id, 11)
        ham_child_2_1 = self.project.create_job('ScriptJob', "child_2_1")
        ham_child_2_1.master_id = ham_master_2.job_id
        ham_child_2_1.save()
        self.assertEqual(ham_child_2_1.job_id, 12)
        self.assertEqual(ham_master_2.child_ids, [12])
        self.assertEqual(ham_master_1.child_ids, [2, 3, 4, 5, 6, 7, 8, 9, 10])
        ham_master_1_reload = self.project.load(1)
        ham_master_2_reload = self.project.load(11)
        self.assertEqual(ham_master_2_reload.child_ids, [12])
        self.assertEqual(ham_master_1_reload.child_ids, [2, 3, 4, 5, 6, 7, 8, 9, 10])
