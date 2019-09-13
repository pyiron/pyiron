# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
from pyiron.base.generic.util import static_isinstance


class TestJobType(unittest.TestCase):
    def test_static_isinstance(self):
        self.assertTrue(
            static_isinstance(
                obj=list(), obj_type=["builtins.list", "__builtin__.list"]
            )
        )
        self.assertTrue(
            any(
                [
                    static_isinstance(obj=list(), obj_type="builtins.list"),
                    static_isinstance(obj=list(), obj_type="__builtin__.list"),
                ]
            )
        )
        self.assertRaises(TypeError, static_isinstance, list(), 1)


if __name__ == "__main__":
    unittest.main()
