import unittest
from pyiron.base.generic.util import static_isinstance


class TestJobType(unittest.TestCase):
    def test_static_isinstance(self):
        self.assertTrue(static_isinstance(obj=list(), obj_type=['builtins.list', '__builtin__.list']))
        self.assertTrue(any([static_isinstance(obj=list(), obj_type='builtins.list'),
                             static_isinstance(obj=list(), obj_type='__builtin__.list')]))
        self.assertRaises(TypeError, static_isinstance, list(), 1)


if __name__ == '__main__':
    unittest.main()
