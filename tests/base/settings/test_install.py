import os
import shutil
import subprocess
from pyiron.base.settings.install import command_line
import unittest


class TestInstall(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.install_script = os.path.join(cls.execution_path, '../../../pyiron/base/settings/install.py')

    @classmethod
    def tearDownClass(cls):
        execution_path = os.path.dirname(os.path.abspath(__file__))
        shutil.rmtree(os.path.join(execution_path, 'resources'))
        shutil.rmtree(os.path.join(execution_path, 'project'))
        os.remove(os.path.join(execution_path, 'config'))

    def test_install(self):
        command_line(['-c', os.path.join(self.execution_path, 'config'),
                      '-r', os.path.join(self.execution_path, 'resources'),
                      '-p', os.path.join(self.execution_path, 'project')])

        with open(os.path.join(self.execution_path, 'config'), 'r') as f:
            content = f.readlines()
        self.assertEqual(content[0], '[DEFAULT]\n')
        self.assertIn('PROJECT_PATHS', content[1])
        self.assertIn('RESOURCE_PATHS', content[2])
        self.assertTrue(os.path.exists(os.path.join(self.execution_path, 'project')))
        self.assertTrue(os.path.exists(os.path.join(self.execution_path, 'resources')))

    def test_install_help(self):
        out = subprocess.check_output(['python', self.install_script, '--error'],
                                      cwd=self.execution_path, shell=False, universal_newlines=True)
        self.assertEqual(out, 'install.py -c <config_file> -p <project_path> -r <resource_dir> -u <url>\n')
        out = subprocess.check_output(['python', self.install_script, '-h'],
                                      cwd=self.execution_path, shell=False, universal_newlines=True)
        self.assertEqual(out, 'install.py -c <config_file> -p <project_path> -r <resource_dir> -u <url>\n')


if __name__ == '__main__':
    unittest.main()
