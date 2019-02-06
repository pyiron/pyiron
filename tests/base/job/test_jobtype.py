import unittest
from pyiron.base.job.jobtype import static_isinstance


class TestJobType(unittest.TestCase):
    def test_static_isinstance(self):
        self.assertTrue(static_isinstance(obj=list(), obj_type=['builtins.list', 'builtins.object']))
        self.assertTrue(static_isinstance(obj=list(), obj_type='builtins.list'))
        self.assertTrue(static_isinstance(obj=list().__class__, obj_type='builtins.list'))
        with self.assertRaises(TypeError) as context:
            self.assertTrue(static_isinstance(obj=list(), obj_type=1))


if __name__ == '__main__':
    unittest.main()
