import getpass
import os
from pyiron_base.core.settings.config.default import ConfigDefault
from pyiron_base.core.settings.generic import Settings
import unittest


class TestConfigSettings(unittest.TestCase):
    def setUp(self):
        self.config = Settings(config=ConfigDefault())
        if os.name == 'nt':
            self.user_home = os.path.expanduser('~').replace('\\', '/')
        else:
            self.user_home = os.path.expanduser('~')
        self.pyiron_path = self.config.path_pyiron
        self.top_level_dir_1 = self.pyiron_path + '/examples/'
        self.top_level_dir_2 = self.pyiron_path + '/tests/'
        self.user_name = getpass.getuser()

    def test_login_user(self):
        self.assertEqual(self.config.login_user, self.user_name)

    # def test_path_pyiron(self):
    #     if os.name == 'nt':
    #         posix_path = os.path.realpath(__file__).replace('\\', '/')
    #     else:
    #         posix_path = os.path.realpath(__file__)
    #     self.assertEqual(self.config.path_pyiron + 'base', '/'.join(posix_path.split('/')[:-4]))

    def test_path_bin(self):
        self.assertEqual(self.config.path_bin, self.user_home + '/PyIron_data/bin')

    def test_path_potentials(self):
        self.assertEqual(self.config.path_potentials, self.pyiron_path + '/static/potentials')

    def test_open_connection(self):
        # print(self.config.open_connection())
        # self.fail()
        pass

    def test_top_path(self):
        self.assertEqual(self.config.top_path(self.top_level_dir_1), self.top_level_dir_1)
        self.assertEqual(self.config.top_path(self.top_level_dir_2), self.top_level_dir_2)


if __name__ == '__main__':
    unittest.main()
