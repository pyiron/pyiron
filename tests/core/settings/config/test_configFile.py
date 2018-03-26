import os
from pyiron_base.core.settings.config.parser import ConfigFile
import unittest


class TestConfigFile(unittest.TestCase):
    def setUp(self):
        self.config = ConfigFile(config_file='pyiron_test_config')

    def test_login_user(self):
        self.assertEqual(self.config.login_user, 'pyiron')

    # def test_pyiron_envs(self):
    #     self.assertEqual(self.config.pyiron_envs,
    #                     [{'top_level_dirs':
    #                           {'/home/pyiron/myData/projects/': '/home/pyiron/myData/projects/',
    #                            '/home/pyiron/myPrograms/PyIron/': '/home/pyiron/myPrograms/PyIron/'},
    #                       'system': 'sqlite',
    #                       'table_name': 'jobs_pyiron',
    #                       'file': '/home/pyiron/myData/pyiron.db?timeout=10',
    #                       'type': 'SQLite'},
    #                      {'host': 'localhost',
    #                       'top_level_dirs': {'/home/pyiron/myData/projects/': '/home/pyiron/myData/projects/',
    #                                          '/home/pyiron/myPrograms/PyIron/': '/home/pyiron/myPrograms/PyIron/'},
    #                       'table_name': 'jobs_pyiron',
    #                       'user': 'postgres',
    #                       'database': 'mdb',
    #                       'password': 'postgres',
    #                       'type': 'Postgres',
    #                       'system': 'postgres'},
    #                      {'host': '127.0.0.1',
    #                       'top_level_dirs': {'/home/pyiron/on/remote/server/': '/home/pyiron/myData/projects/'},
    #                       'table_name': 'jobs_pyiron',
    #                       'user': 'postgres',
    #                       'database': 'mdb',
    #                       'password': 'postgres',
    #                       'type': 'Postgres',
    #                       'system': 'remote'}]
    #                     )

    def test_path_bin(self):
        self.assertEqual(self.config.path_bin, '/home/pyiron/myData/bin')

    def test_path_potentials(self):
        self.assertEqual(self.config.path_potentials, '/home/pyiron/myPrograms/PyIron/static/potentials')

    # def test_path_pyiron(self):
    #     if os.name == 'nt':
    #         posix_path = os.path.realpath(__file__).replace('\\', '/')
    #     else:
    #         posix_path = os.path.realpath(__file__)
    #     self.assertEqual(self.config.path_pyiron, '/'.join(posix_path.split('/')[:-5]))


if __name__ == '__main__':
    unittest.main()
