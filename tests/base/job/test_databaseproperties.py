import unittest
import datetime
import os
from pyiron.base.project.generic import Project
from pyiron.base.job.core import DatabaseProperties


class TestDatabaseProperties(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database_entry = {'id': 150,
                              'parentid': None,
                              'masterid': None,
                              'projectpath': '/Users/jan/PyIron_data/projects/',
                              'project': '2019-02-14-database-properties/test/',
                              'job': 'vasp',
                              'subjob': '/vasp',
                              'chemicalformula': 'Fe2',
                              'status': 'finished',
                              'hamilton': 'Vasp',
                              'hamversion': '5.4',
                              'username': 'pyiron',
                              'computer': 'pyiron@MacBook-Pro-4.local#1',
                              'timestart': datetime.datetime(2019, 2, 14, 8, 4, 7, 248427),
                              'timestop': datetime.datetime(2019, 2, 14, 8, 4, 8, 366365),
                              'totalcputime': 1.0}
        cls.database_property = DatabaseProperties(job_dict=cls.database_entry)

    def test_properties(self):
        self.assertEqual(self.database_property.id, 150)
        self.assertEqual(self.database_property.parentid, None)
        self.assertEqual(self.database_property.masterid, None)
        self.assertEqual(self.database_property.projectpath, '/Users/jan/PyIron_data/projects/')
        self.assertEqual(self.database_property.project, '2019-02-14-database-properties/test/')
        self.assertEqual(self.database_property.job, 'vasp')
        self.assertEqual(self.database_property.subjob, '/vasp')
        self.assertEqual(self.database_property.chemicalformula, 'Fe2')
        self.assertEqual(self.database_property.status, 'finished')
        self.assertEqual(self.database_property.hamilton, 'Vasp')
        self.assertEqual(self.database_property.hamversion, '5.4')
        self.assertEqual(self.database_property.username, 'pyiron')
        self.assertEqual(self.database_property.computer, 'pyiron@MacBook-Pro-4.local#1')
        self.assertEqual(self.database_property.timestart, datetime.datetime(2019, 2, 14, 8, 4, 7, 248427))
        self.assertEqual(self.database_property.timestop, datetime.datetime(2019, 2, 14, 8, 4, 8, 366365))
        self.assertEqual(self.database_property.totalcputime, 1.0)

    def test_dir(self):
        self.assertEqual(sorted(list(self.database_entry.keys())), sorted(dir(self.database_property)))

    def test_bool(self):
        self.assertTrue(bool(self.database_property))
        self.assertFalse(bool(DatabaseProperties()))
        with self.assertRaises(AttributeError):
            _ = DatabaseProperties().job


class DatabasePropertyIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, 'database_prop'))
        cls.ham = cls.project.create_job("ExampleJob", "job_test_run")
        cls.ham.run()

    @classmethod
    def tearDownClass(cls):
        project = Project(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'database_prop'))
        ham = project.load(project.get_job_ids()[0])
        ham.remove()
        project.remove(enable=True)

    def test_properties(self):
        job_db_entry_dict = self.ham.project.db.get_item_by_id(self.ham.job_id)
        self.assertIsNotNone(job_db_entry_dict)
        self.assertEqual(self.ham.database_entry.id, job_db_entry_dict['id'])
        self.assertEqual(self.ham.database_entry.parentid, job_db_entry_dict['parentid'])
        self.assertEqual(self.ham.database_entry.masterid, job_db_entry_dict['masterid'])
        self.assertEqual(self.ham.database_entry.projectpath, self.project.root_path)
        self.assertEqual(self.ham.database_entry.project, self.project.project_path)
        self.assertEqual(self.ham.database_entry.job, 'job_test_run')
        self.assertEqual(self.ham.database_entry.subjob, '/job_test_run')
        self.assertEqual(self.ham.database_entry.status, 'finished')
        self.assertEqual(self.ham.database_entry.hamilton, "ExampleJob")
        self.assertEqual(self.ham.database_entry.hamversion, '0.3')
        self.assertEqual(self.ham.database_entry.username, 'pyiron')

    def test_inspect_job(self):
        job_inspect = self.project.inspect(self.ham.job_name)
        self.assertIsNotNone(job_inspect)
        self.assertEqual(job_inspect.database_entry.parentid, None)
        self.assertEqual(job_inspect.database_entry.masterid, None)
        self.assertEqual(job_inspect.database_entry.projectpath, self.project.root_path)
        self.assertEqual(job_inspect.database_entry.project, self.project.project_path)
        self.assertEqual(job_inspect.database_entry.job, 'job_test_run')
        self.assertEqual(job_inspect.database_entry.subjob, '/job_test_run')
        self.assertEqual(job_inspect.database_entry.status, 'finished')
        self.assertEqual(job_inspect.database_entry.hamilton, "ExampleJob")
        self.assertEqual(job_inspect.database_entry.hamversion, '0.3')
        self.assertEqual(job_inspect.database_entry.username, 'pyiron')

if __name__ == '__main__':
    unittest.main()
