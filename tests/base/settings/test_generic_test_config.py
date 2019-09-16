# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pathlib2 import Path
import os
import unittest
from pyiron.base.settings.generic import Settings


class TestConfigSettingsStatic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.resource_path = (
            Path(__file__)
            .expanduser()
            .resolve()
            .absolute()
            .as_posix()
            .replace("\\", "/")
        )
        cls.test_config = Settings(
            config={
                "sql_file": "sqlite.db",
                "project_paths": os.path.join(cls.resource_path, "../../../../.."),
                "resource_paths": os.path.join(cls.resource_path, "../../../../static"),
            }
        )

    def test_get_config_from_environment(self):
        config = self.test_config.get_config_from_environment(environment={"PYIRONSQLFILE": '/a/b/c',
                                                                  "SYSTEM": 'linux'},
                                                              config={'user': 'pyiron'})
        self.assertEqual(config['sql_file'], '/a/b/c')
        self.assertEqual(config['user'], 'pyiron')
        self.assertEqual(len(config), 2)

    # def test_db_connection_name(self):
    #     self.assertEqual(self.test_config.db_connection_name, 'test')
    #
    # def test_db_connection_string(self):
    #     self.assertEqual(self.test_config.db_connection_string, 'sqlite:///sqlite.db')
    #
    # def test_db_connection_table(self):
    #     self.assertEqual(self.test_config.db_connection_table, 'jobs_pyiron')

    # def test_db_translate_dict(self):
    #     self.assertEqual(self.test_config.db_translate_dict,
    #                      {'test': {self.resource_path + '/': self.resource_path + '/'}})

    # def test_db_name(self):
    #     self.assertEqual(self.test_config.db_name, 'test')

    # def test_top_path(self):
    #     self.assertEqual(self.test_config.top_path(self.resource_path + '/test'),
    #                      self.resource_path + '/')

    def test_resource_paths(self):
        self.assertEqual(
            self.test_config.resource_paths,
            [
                os.path.abspath(
                    os.path.join(self.resource_path, "../../../../static")
                ).replace("\\", "/")
            ],
        )

    def test_login_user(self):
        self.assertEqual(self.test_config.login_user, "pyiron")


if __name__ == "__main__":
    unittest.main()
