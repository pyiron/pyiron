import unittest
from base.core.project.readmode import ReadMode


class TestReadMode(unittest.TestCase):
    def setUp(self):
        self.readmode_default = ReadMode()
        self.readmode_inspect = ReadMode(mode='inspect')
        self.readmode_object = ReadMode(mode='object')
        self.readmode_modify = ReadMode()

    def test_inspect(self):
        self.assertFalse(self.readmode_default.inspect)

        self.assertTrue(self.readmode_inspect.inspect)
        self.assertFalse(self.readmode_inspect.object)
        self.readmode_inspect.inspect = False
        self.assertFalse(self.readmode_inspect.inspect)
        self.assertTrue(self.readmode_inspect.object)
        self.readmode_inspect.inspect = True
        self.assertTrue(self.readmode_inspect.inspect)
        self.assertFalse(self.readmode_inspect.object)

    def test_object(self):
        self.assertTrue(self.readmode_default.object)

        self.assertFalse(self.readmode_object.inspect)
        self.assertTrue(self.readmode_object.object)
        self.readmode_object.object = False
        self.assertTrue(self.readmode_object.inspect)
        self.assertFalse(self.readmode_object.object)
        self.readmode_object.object = True
        self.assertFalse(self.readmode_object.inspect)
        self.assertTrue(self.readmode_object.object)

    def test_mode(self):
        self.assertFalse(self.readmode_modify.inspect)
        self.assertTrue(self.readmode_modify.object)
        self.readmode_modify.mode = 'object'
        self.assertFalse(self.readmode_modify.inspect)
        self.assertTrue(self.readmode_modify.object)
        self.readmode_modify.mode = 'inspect'
        self.assertTrue(self.readmode_modify.inspect)
        self.assertFalse(self.readmode_modify.object)
        self.readmode_modify.mode = 'object'
        self.assertFalse(self.readmode_modify.inspect)
        self.assertTrue(self.readmode_modify.object)


if __name__ == '__main__':
    unittest.main()
