import os
import shutil
import subprocess
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
        subprocess.check_output(['python', self.install_script,
                                 '-c', os.path.join(self.execution_path, 'config'),
                                 '-r', os.path.join(self.execution_path, 'resources'),
                                 '-p', os.path.join(self.execution_path, 'project')],
                                cwd=self.execution_path, shell=False)
        with open(os.path.join(self.execution_path, 'config'), 'r') as f:
            content = f.readlines()
        self.assertEqual(content[0], '[DEFAULT]\n')
        self.assertIn('PROJECT_PATHS', content[1])
        self.assertIn('RESOURCE_PATHS', content[2])
        self.assertTrue(os.path.exists(os.path.join(self.execution_path, 'project')))
        self.assertTrue(os.path.exists(os.path.join(self.execution_path, 'resources')))


if __name__ == '__main__':
    unittest.main()
