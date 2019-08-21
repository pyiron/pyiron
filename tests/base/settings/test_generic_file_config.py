# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.base.settings.generic import Settings
import os
import unittest


class TestConfigSettingsStatic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                         '../../static')).replace('\\', '/')
        cls.project_path = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
        cls.file_config = Settings()

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
        self.assertTrue(self.file_config.top_path(self.project_path + '/test') in self.project_path)

    def test_file_resource_paths(self):
        self.assertTrue(any([path for path in self.file_config.resource_paths if path in self.resource_path]))

    def test_file_login_user(self):
        self.assertEqual(self.file_config.login_user, 'pyiron')


if __name__ == '__main__':
    unittest.main()
