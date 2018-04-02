from pathlib2 import Path
import unittest
from pyiron_base.core.settings.generic import Settings


class TestConfigSettingsStatic(unittest.TestCase):
    def setUp(self):
        self.resource_path = Path('.').expanduser().resolve().absolute().as_posix()
        self.test_config = Settings(config={'file': 'genericsettings.db',
                                            'top_level_dirs': self.resource_path,
                                            'resource_paths': self.resource_path})

    def test_db_connection_name(self):
        self.assertEqual(self.test_config.db_connection_name, 'test')

    def test_db_connection_string(self):
        self.assertEqual(self.test_config.db_connection_string, 'sqlite:///genericsettings.db')

    def test_db_connection_table(self):
        self.assertEqual(self.test_config.db_connection_table, 'jobs_pyiron')

    def test_db_translate_dict(self):
        self.assertEqual(self.test_config.db_translate_dict,
                         {'test': {self.resource_path + '/': self.resource_path + '/'}})

    def test_db_name(self):
        self.assertEqual(self.test_config.db_name, 'test')

    def test_top_path(self):
        self.assertEqual(self.test_config.top_path(self.resource_path + '/test'),
                         self.resource_path + '/')

    def test_resource_paths(self):
        self.assertEqual(self.test_config.resource_paths, [self.resource_path])

    def test_login_user(self):
        self.assertEqual(self.test_config.login_user, 'pyiron')


if __name__ == '__main__':
    unittest.main()
