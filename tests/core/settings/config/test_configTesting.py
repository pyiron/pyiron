import unittest
from pyiron_base.core.settings.config.testing import ConfigTesting
from pyiron_base.core.settings.generic import Settings


class TestConfigDefault(unittest.TestCase):
    def setUp(self):
        self.login_user = 'test'
        self.sql_lite_database = './config_testing.db'
        self.path_bin = '/path/to/binaries'
        self.path_potentials = '/path/to/potentials'
        self.path_project = '/path/to/testing/project/'
        self.code = '/path/to/source/code'
        self.config = ConfigTesting(login_user='test',
                                    sql_lite_database='./config_testing.db',
                                    path_bin='/path/to/binaries',
                                    path_potentials='/path/to/potentials',
                                    path_project='/path/to/testing/project',
                                    path_sourcecode='/path/to/source/code')

    def test_pyiron_envs(self):
        self.assertEqual(self.config.pyiron_envs,
                         {'top_level_dirs': {self.path_project: self.path_project},
                          'system': 'default',
                          'table_name': 'jobs_' + self.login_user,
                          'file': self.sql_lite_database,
                          'type': 'SQLite',
                          'user': self.login_user,
                          'resource_paths': ['/path/to/potentials', '/path/to/binaries']})

    def test_path_bin(self):
        self.assertEqual(self.config.path_bin, self.path_bin)

    def test_path_potentials(self):
        self.assertEqual(self.config.path_potentials, self.path_potentials)

    def test_path_pyiron(self):
        self.assertEqual(self.config.path_pyiron, self.code)


class TestConfigSettings(unittest.TestCase):
    def setUp(self):
        self.login_user = 'test'
        self.sql_lite_database = './config_testing.db'
        self.path_bin = '/path/to/binaries'
        self.path_potentials = '/path/to/potentials'
        self.path_project = '/path/to/testing/project/'
        self.code = '/path/to/source/code'
        config = ConfigTesting(login_user='test',
                                    sql_lite_database='./config_testing.db',
                                    path_bin='/path/to/binaries',
                                    path_potentials='/path/to/potentials',
                                    path_project='/path/to/testing/project',
                                    path_sourcecode='/path/to/source/code')
        self.config = Settings(config=config)

    def test_login_user(self):
        self.assertEqual(self.config.login_user, self.login_user)

    def test_path_bin(self):
        self.assertEqual(self.config.path_bin, self.path_bin)

    def test_path_potentials(self):
        self.assertEqual(self.config.path_potentials, self.path_potentials)

    def test_path_pyiron(self):
        self.assertEqual(self.config.path_pyiron, self.code)

    def test_open_connection(self):
        # print(self.config.open_connection())
        # self.fail()
        pass

    def test_top_path(self):
        self.assertEqual(self.config.top_path(self.path_project), self.path_project)


if __name__ == '__main__':
    unittest.main()