import os
from pyiron_base.core.settings.generic import Settings
import unittest

s = Settings(config={'sql_file': 'genericpath.db',
                     'project_paths': os.path.abspath(os.getcwd()),
                     'resource_paths': os.path.abspath(os.getcwd())})

from pyiron_base.project import Project
from pyiron_base.core.project.path import GenericPath


class TestGenericPath(unittest.TestCase):
    def setUp(self):
        self.current_dir = os.getcwd().replace('\\', '/')
        self.path_project = GenericPath(root_path=self.current_dir,
                                        project_path='project/path/')

    def test_root_path(self):
        self.assertEqual(self.path_project.root_path, self.current_dir + '/')

    def test_project_path(self):
        self.assertEqual(self.path_project.project_path, 'project/path/')

    def test_path(self):
        self.assertEqual(self.path_project.path, self.current_dir + '/project/path/')


class TestProject(unittest.TestCase):
    def setUp(self):
        self.current_dir = os.getcwd()
        self.project = Project('sub_folder')

    def tearDown(self):
        self.project.remove()

    def test_repr(self):
        self.assertEqual("{'nodes': [], 'groups': []}", self.project.__repr__())
        pr_down_one = self.project['..']
        pr_down_two = self.project['../..']
        pr_down_twice = self.project['..']['..']
        self.assertEqual(pr_down_two.__repr__(), pr_down_twice.__repr__())
        self.assertEqual(str(sorted([directory for directory in os.listdir('.')
                                     if not os.path.isfile(os.path.join('.', directory))])),
                         pr_down_one.__repr__()['groups'])
        self.assertEqual(str(sorted([directory for directory in os.listdir('..')
                                     if not os.path.isfile(os.path.join('..', directory))])),
                         pr_down_two.__repr__()['groups'])
        self.assertEqual(pr_down_two.__repr__()['groups'], pr_down_twice.__repr__()['groups'])

if __name__ == '__main__':
    unittest.main()