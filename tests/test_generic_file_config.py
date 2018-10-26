from pyiron_base.core.settings.generic import Settings
from pathlib2 import Path
import os
import sys
import unittest


class TestConfigSettingsStatic(unittest.TestCase):
    def setUp(self):
        self.resource_path = os.path.dirname(os.path.abspath(__file__))
        self.project_path = os.path.dirname(os.path.abspath(__file__)) + '/'
        self.file_config = Settings()

    # def test_file_db_connection_name(self):
    #     self.assertEqual(self.file_config.db_connection_name, 'DEFAULT')
    #
    # def test_file_db_connection_string(self):
    #     self.assertEqual(self.file_config.db_connection_string, 'sqlite:///' + self.resource_path + '/sqlite.db')
    #
    # def test_file_db_connection_table(self):
    #     self.assertEqual(self.file_config.db_connection_table, 'jobs_pyiron')

    # def test_file_db_translate_dict(self):
    #     self.assertEqual(self.file_config.db_translate_dict,
    #                      {'DEFAULT': {self.project_path: self.project_path}})

    # def test_file_db_name(self):
    #     self.assertEqual(self.file_config.db_name, 'DEFAULT')

    def test_file_top_path(self):
        self.assertEqual(self.file_config.top_path(self.project_path + '/test'), self.project_path)

    def test_file_resource_paths(self):
        self.assertEqual(self.file_config.resource_paths, [self.resource_path])

    def test_file_login_user(self):
        self.assertEqual(self.file_config.login_user, 'pyiron')


if __name__ == '__main__':
    unittest.main()
