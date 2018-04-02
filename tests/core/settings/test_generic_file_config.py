from pyiron_base.core.settings.generic import Settings
from pathlib2 import Path
import os
import sys
import unittest


class TestConfigSettingsStatic(unittest.TestCase):
    def setUp(self):
        if sys.version_info.major < 3 and os.name == 'nt':
            # In Python 2.7 on Windows for pathlib2 it is required that the directories exist, so we create them
            if not os.path.exists(os.path.expanduser('~/pyiron/resources')):
                os.makedirs(os.path.expanduser('~/pyiron/resources'))
            if not os.path.exists(os.path.expanduser('~/pyiron/projects')):
                os.makedirs(os.path.expanduser('~/pyiron/projects'))
        self.user_path = Path('~').expanduser().resolve().absolute().as_posix()
        self.resource_path = Path('~/pyiron/resources').expanduser().resolve().absolute().as_posix()
        self.project_path = Path('~/pyiron/projects').expanduser().resolve().absolute().as_posix() + '/'
        self.file_config = Settings()

    def test_file_db_connection_name(self):
        self.assertEqual(self.file_config.db_connection_name, 'DEFAULT')

    def test_file_db_connection_string(self):
        self.assertEqual(self.file_config.db_connection_string, 'sqlite:///' + self.resource_path + '/sqlite.db')

    def test_file_db_connection_table(self):
        self.assertEqual(self.file_config.db_connection_table, 'jobs_pyiron')

    def test_file_db_translate_dict(self):
        self.assertEqual(self.file_config.db_translate_dict,
                         {'DEFAULT': {self.project_path: self.project_path}})

    def test_file_db_name(self):
        self.assertEqual(self.file_config.db_name, 'DEFAULT')

    def test_file_top_path(self):
        self.assertEqual(self.file_config.top_path(self.project_path + '/test'), self.project_path)

    def test_file_resource_paths(self):
        self.assertEqual(self.file_config.resource_paths, [self.resource_path])

    def test_file_login_user(self):
        self.assertEqual(self.file_config.login_user, 'pyiron')


if __name__ == '__main__':
    unittest.main()
