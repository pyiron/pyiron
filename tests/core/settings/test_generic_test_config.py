import os
from pyiron_base.core.settings.generic import Settings
import unittest


class TestConfigSettingsStatic(unittest.TestCase):
    def setUp(self):
        self.test_config = Settings(config={'file': 'genericsettings.db',
                                            'top_level_dirs': os.path.abspath(os.getcwd()),
                                            'resource_paths': os.path.abspath(os.getcwd())})

    def test_db_connection_name(self):
        self.assertEqual(self.test_config.db_connection_name, 'test')

    def test_db_connection_string(self):
        self.assertEqual(self.test_config.db_connection_string, 'sqlite:///genericsettings.db')

    def test_db_connection_table(self):
        self.assertEqual(self.test_config.db_connection_table, 'jobs_pyiron')

    def test_db_translate_dict(self):
        self.assertEqual(self.test_config.db_translate_dict,
                         {'test': {os.path.abspath(os.getcwd()) + '/': os.path.abspath(os.getcwd()) + '/'}})

    def test_db_name(self):
        self.assertEqual(self.test_config.db_name, 'test')

    def test_top_path(self):
        self.assertEqual(self.test_config.top_path(os.path.abspath(os.getcwd()) + '/test'),
                         os.path.abspath(os.getcwd()) + '/')

    def test_resource_paths(self):
        self.assertEqual(self.test_config.resource_paths, [os.path.abspath(os.getcwd())])

    def test_login_user(self):
        self.assertEqual(self.test_config.login_user, 'pyiron')


if __name__ == '__main__':
    unittest.main()
