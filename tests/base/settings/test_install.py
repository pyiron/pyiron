# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import shutil
import subprocess
from pyiron.base.settings.install import install_pyiron
import unittest


class TestInstall(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))

    @classmethod
    def tearDownClass(cls):
        execution_path = os.path.dirname(os.path.abspath(__file__))
        shutil.rmtree(os.path.join(execution_path, "resources"))
        shutil.rmtree(os.path.join(execution_path, "project"))
        os.remove(os.path.join(execution_path, "config"))

    def test_install(self):
        install_pyiron(
                config_file_name = os.path.join(self.execution_path, "config"),
                resource_directory = os.path.join(self.execution_path, "resources"),
                project_path = os.path.join(self.execution_path, "project"),
        )

        with open(os.path.join(self.execution_path, "config"), "r") as f:
            content = f.readlines()
        self.assertEqual(content[0], "[DEFAULT]\n")
        self.assertIn("PROJECT_PATHS", content[1])
        self.assertIn("RESOURCE_PATHS", content[2])
        self.assertTrue(os.path.exists(os.path.join(self.execution_path, "project")))
        self.assertTrue(os.path.exists(os.path.join(self.execution_path, "resources")))

if __name__ == "__main__":
    unittest.main()
