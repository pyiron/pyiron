import getpass
import os
from pyiron_base.core.settings.config.default import ConfigDefault
import unittest


class TestConfigDefault(unittest.TestCase):
    def setUp(self):
        self.config = ConfigDefault()
        self.user_home = os.path.expanduser('~').replace('\\', '/')
        self.pyiron_path = self.config.path_pyiron
        self.top_level_dir_1 = self.pyiron_path + '/examples/'
        self.top_level_dir_2 = self.pyiron_path + '/tests/'
        self.user_name = getpass.getuser()

    def test_pyiron_envs(self):
        self.assertEqual(self.config.pyiron_envs,
                         {'top_level_dirs': {self.top_level_dir_1: self.top_level_dir_1,
                                             self.top_level_dir_2: self.top_level_dir_2},
                          'system': 'default',
                          'user': self.user_name,
                          'resource_paths': [self.pyiron_path + '/static'],
                          'table_name': 'jobs_' + self.user_name,
                          'file': self.user_home + '/PyIron_data/pyiron.db',
                          'type': 'SQLite'})

    def test_path_bin(self):
        self.assertEqual(self.config.path_bin, self.user_home + '/PyIron_data/bin')

    def test_path_potentials(self):
        self.assertEqual(self.config.path_potentials, self.pyiron_path + '/static/potentials')

    # def test_path_pyiron(self):
    #     sep = os.path.sep
    #     new_path = sep.join(os.path.realpath(__file__).split(sep)[:-5])
    #     self.assertEqual(self.config.path_pyiron + 'base', new_path.replace('\\', '/'))
    #     # self.assertEqual(self.config.path_pyiron, '/'.join(posixpath.realpath(__file__).split('/')[:-5]))


if __name__ == '__main__':
    unittest.main()
