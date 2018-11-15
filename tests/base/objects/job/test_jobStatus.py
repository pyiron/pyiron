import os
from datetime import datetime
from base.core.settings.database import DatabaseAccess
from base.objects.job.jobstatus import JobStatus
import unittest


class TestJobStatus(unittest.TestCase):
    def setUp(self):
        self.jobstatus = JobStatus()
        self.database = DatabaseAccess('sqlite:///test_job_status.db', 'simulation')
        par_dict = {'chemicalformula': 'H',
                    'computer': 'localhost',
                    'hamilton': 'Test',
                    'hamversion': '0.1',
                    'job': 'testing',
                    'parentid': 0,
                    'project': 'database.testing',
                    'projectpath': '/TESTING',
                    'status': 'initialized',
                    'timestart': datetime(2016, 5, 2, 11, 31, 4, 253377),
                    'timestop': datetime(2016, 5, 2, 11, 31, 4, 371165),
                    'totalcputime': 0.117788,
                    'username': 'Test'}
        self.job_id = self.database.add_item_dict(par_dict)
        self.jobstatus_database = JobStatus(db=self.database, job_id=self.job_id)

    def doCleanups(self):
        self.database.conn.close()
        os.remove('test_job_status.db')

    def test_initialized(self):
        self.assertTrue(self.jobstatus.initialized)
        self.jobstatus.string = 'finished'
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertTrue(self.jobstatus.finished)
        self.jobstatus.initialized = True
        self.assertTrue(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_appended(self):
        self.jobstatus.appended = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertTrue(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_created(self):
        self.jobstatus.created = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertTrue(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_submitted(self):
        self.jobstatus.submitted = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertTrue(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_running(self):
        self.jobstatus.running = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertTrue(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_aborted(self):
        self.jobstatus.aborted = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertTrue(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_collect(self):
        self.jobstatus.collect = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertTrue(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_suspended(self):
        self.jobstatus.suspended = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertTrue(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_refresh(self):
        self.jobstatus.refresh = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertTrue(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_busy(self):
        self.jobstatus.busy = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertTrue(self.jobstatus.busy)
        self.assertFalse(self.jobstatus.finished)

    def test_finished(self):
        self.jobstatus.finished = True
        self.assertFalse(self.jobstatus.initialized)
        self.assertFalse(self.jobstatus.appended)
        self.assertFalse(self.jobstatus.created)
        self.assertFalse(self.jobstatus.submitted)
        self.assertFalse(self.jobstatus.running)
        self.assertFalse(self.jobstatus.aborted)
        self.assertFalse(self.jobstatus.collect)
        self.assertFalse(self.jobstatus.suspended)
        self.assertFalse(self.jobstatus.refresh)
        self.assertFalse(self.jobstatus.busy)
        self.assertTrue(self.jobstatus.finished)

    def test_string(self):
        self.jobstatus.string = 'initialized'
        self.assertTrue(self.jobstatus.initialized)
        self.assertEqual(str(self.jobstatus), 'initialized')
        self.assertEqual(self.jobstatus.string, 'initialized')
        self.jobstatus.string = 'appended'
        self.assertTrue(self.jobstatus.appended)
        self.assertEqual(str(self.jobstatus), 'appended')
        self.assertEqual(self.jobstatus.string, 'appended')
        self.jobstatus.string = 'created'
        self.assertTrue(self.jobstatus.created)
        self.assertEqual(str(self.jobstatus), 'created')
        self.assertEqual(self.jobstatus.string, 'created')
        self.jobstatus.string = 'submitted'
        self.assertTrue(self.jobstatus.submitted)
        self.assertEqual(str(self.jobstatus), 'submitted')
        self.assertEqual(self.jobstatus.string, 'submitted')
        self.jobstatus.string = 'running'
        self.assertTrue(self.jobstatus.running)
        self.assertEqual(str(self.jobstatus), 'running')
        self.assertEqual(self.jobstatus.string, 'running')
        self.jobstatus.string = 'aborted'
        self.assertTrue(self.jobstatus.aborted)
        self.assertEqual(str(self.jobstatus), 'aborted')
        self.assertEqual(self.jobstatus.string, 'aborted')
        self.jobstatus.string = 'collect'
        self.assertTrue(self.jobstatus.collect)
        self.assertEqual(str(self.jobstatus), 'collect')
        self.assertEqual(self.jobstatus.string, 'collect')
        self.jobstatus.string = 'suspended'
        self.assertTrue(self.jobstatus.suspended)
        self.assertEqual(str(self.jobstatus), 'suspended')
        self.assertEqual(self.jobstatus.string, 'suspended')
        self.jobstatus.string = 'refresh'
        self.assertTrue(self.jobstatus.refresh)
        self.assertEqual(str(self.jobstatus), 'refresh')
        self.assertEqual(self.jobstatus.string, 'refresh')
        self.jobstatus.string = 'busy'
        self.assertTrue(self.jobstatus.busy)
        self.assertEqual(str(self.jobstatus), 'busy')
        self.assertEqual(self.jobstatus.string, 'busy')
        self.jobstatus.string = 'finished'
        self.assertTrue(self.jobstatus.finished)
        self.assertEqual(str(self.jobstatus), 'finished')
        self.assertEqual(self.jobstatus.string, 'finished')

    def test_database_connection(self):
        current_status = self.database.get_item_by_id(self.job_id)["status"]
        self.assertTrue(self.jobstatus_database.initialized)
        self.assertEqual(current_status, str(self.jobstatus_database))
        self.jobstatus_database.created = True
        new_status = self.database.get_item_by_id(self.job_id)["status"]
        self.assertTrue(self.jobstatus_database.created)
        self.assertNotEqual(current_status, str(self.jobstatus_database))
        self.assertEqual(new_status, str(self.jobstatus_database))
        self.database.item_update({'status': 'finished'}, self.job_id)
        finished_status = self.database.get_item_by_id(self.job_id)["status"]
        self.assertTrue(self.jobstatus_database.finished)
        self.assertNotEqual(current_status, str(self.jobstatus_database))
        self.assertNotEqual(new_status, str(self.jobstatus_database))
        self.assertEqual(finished_status, str(self.jobstatus_database))


if __name__ == '__main__':
    unittest.main()