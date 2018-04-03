import os
from pyiron_base.core.settings.generic import Settings
import unittest

s = Settings(config={'sql_file': 'projectpath.db',
                     'project_paths': os.path.abspath(os.getcwd()),
                     'resource_paths': os.path.abspath(os.getcwd())})

from pyiron_base.core.project.path import ProjectPath


class TestProjectPath(unittest.TestCase):
    def setUp(self):
        if os.name == 'nt':
            self.current_dir = os.getcwd().replace('\\', '/')
        else:
            self.current_dir = os.getcwd()
        self.project_path = ProjectPath(path=self.current_dir)
        self.project_path = self.project_path.open('test_project_path')

    def test_open(self):
        with self.project_path.open('test_open') as test_open:
            self.assertEqual(test_open.path, self.current_dir + '/test_project_path/test_open/')
        self.project_path.removedirs('test_open')

    def test_close(self):
        with self.project_path.open('test_close') as test_close:
            self.assertEqual(test_close.path, self.current_dir + '/test_project_path/test_close/')
        self.assertEqual(self.project_path.path, self.current_dir + '/test_project_path/')
        self.project_path.removedirs('test_close')

    def test_copy(self):
        with self.project_path.open('test_copy') as test_copy:
            copied_path = test_copy.copy()
            self.assertEqual(copied_path.path, test_copy.path)
        self.project_path.removedirs('test_copy')

    def test_removedirs(self):
        self.project_path = self.project_path.open('test_removedirs')
        self.project_path = self.project_path.open('..')
        self.assertTrue('test_removedirs' in self.project_path.listdir())
        self.project_path.removedirs('test_removedirs')
        self.project_path.close()
        self.assertFalse('test_removedirs' in self.project_path.listdir())


if __name__ == '__main__':
    unittest.main()
