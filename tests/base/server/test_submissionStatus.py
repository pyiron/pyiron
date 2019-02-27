from datetime import datetime
from pyiron.base.database.generic import DatabaseAccess
from pyiron.base.master.submissionstatus import SubmissionStatus
import unittest
import os


class TestSubmissionStatus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sub_status = SubmissionStatus()
        cls.database = DatabaseAccess('sqlite:///test_sub_status.db', 'simulation')
        par_dict = {'chemicalformula': 'H',
                    'computer': 'localhost#1#3',
                    'hamilton': 'Test',
                    'hamversion': '0.1',
                    'job': 'testing',
                    'parentid': 0,
                    'project': 'database.testing',
                    'projectpath': '/TESTING',
                    'status': 'suspended',
                    'timestart': datetime(2016, 5, 2, 11, 31, 4, 253377),
                    'timestop': datetime(2016, 5, 2, 11, 31, 4, 371165),
                    'totalcputime': 0.117788,
                    'username': 'Test'}
        cls.job_id = cls.database.add_item_dict(par_dict)
        cls.sub_status_database = SubmissionStatus(db=cls.database, job_id=cls.job_id)

    @classmethod
    def tearDownClass(cls):
        os.remove('test_sub_status.db')

    def test_submit_next(self):
        before = self.sub_status.submitted_jobs
        self.sub_status.submit_next()
        self.assertEqual(before + 1, self.sub_status.submitted_jobs)
        before = self.sub_status_database.submitted_jobs
        self.sub_status_database.submit_next()
        self.assertEqual(before + 1, self.sub_status_database.submitted_jobs)

    def test_refresh(self):
        self.sub_status.submitted_jobs = 5
        self.sub_status.refresh()
        self.assertEqual(5, self.sub_status.submitted_jobs)
        self.sub_status_database.submitted_jobs = 5
        self.sub_status_database.refresh()
        self.assertEqual(5, self.sub_status_database.submitted_jobs)
