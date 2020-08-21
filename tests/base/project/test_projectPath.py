# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron.base.project.path import ProjectPath


class TestProjectPath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if os.name == "nt":
            cls.current_dir = os.path.dirname(os.path.abspath(__file__)).replace(
                "\\", "/"
            )
        else:
            cls.current_dir = os.path.dirname(os.path.abspath(__file__))
        cls.project_path = ProjectPath(path=cls.current_dir)
        cls.project_path = cls.project_path.open("test_project_path")

    def test_open(self):
        with self.project_path.open("test_open") as test_open:
            self.assertEqual(
                test_open.path, self.current_dir + "/test_project_path/test_open/"
            )
        self.project_path.removedirs("test_open")

    def test_close(self):
        with self.project_path.open("test_close") as test_close:
            self.assertEqual(
                test_close.path, self.current_dir + "/test_project_path/test_close/"
            )
        self.assertEqual(
            self.project_path.path, self.current_dir + "/test_project_path/"
        )
        self.project_path.removedirs("test_close")

    def test_copy(self):
        with self.project_path.open("test_copy") as test_copy:
            copied_path = test_copy.copy()
            self.assertEqual(copied_path.path, test_copy.path)
        self.project_path.removedirs("test_copy")

    def test_removedirs(self):
        self.project_path = self.project_path.open("test_removedirs")
        self.project_path = self.project_path.open("..")
        self.assertTrue("test_removedirs" in self.project_path.listdir())
        self.project_path.removedirs("test_removedirs")
        self.project_path.close()
        self.assertFalse("test_removedirs" in self.project_path.listdir())

if __name__ == "__main__":
    unittest.main()
