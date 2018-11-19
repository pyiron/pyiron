import unittest
from pyiron.base.objects.server.shelloption import ShellOption


class TestShellOption(unittest.TestCase):
    def setUp(self):
        self.shelloption0 = ShellOption(name='active_key_and_value', key='key', value='value', value_prefix='prefix=')
        self.shelloption1 = ShellOption(name='active_key_and_value', key='key', value='value')
        self.shelloption2 = ShellOption(name='active_key_only', key='key')
        self.shelloption3 = ShellOption(name='active_value_only', value='value')
        self.shelloption4 = ShellOption(name='passive', active=False)

    def test_active(self):
        self.assertTrue(self.shelloption0.active)
        self.assertTrue(self.shelloption1.active)
        self.assertTrue(self.shelloption2.active)
        self.assertTrue(self.shelloption3.active)
        self.assertFalse(self.shelloption4.active)

    def test_key(self):
        self.assertEqual(self.shelloption0.key, 'key')
        self.assertEqual(self.shelloption1.key, 'key')
        self.assertEqual(self.shelloption2.key, 'key')
        self.assertFalse(self.shelloption3.key)
        self.assertFalse(self.shelloption4.key)

    def test_value(self):
        self.assertEqual(self.shelloption0.value, 'value')
        self.assertEqual(self.shelloption1.value, 'value')
        self.assertFalse(self.shelloption2.value)
        self.assertEqual(self.shelloption3.value, 'value')
        self.assertFalse(self.shelloption4.value)

    def test_value_prefix(self):
        self.assertEqual(self.shelloption0.value_prefix, 'prefix=')
        self.assertFalse(self.shelloption1.value_prefix)
        self.assertFalse(self.shelloption2.value_prefix)
        self.assertFalse(self.shelloption3.value_prefix)
        self.assertFalse(self.shelloption4.value_prefix)

    def test_value_only(self):
        self.assertFalse(self.shelloption0.value_only)
        self.assertFalse(self.shelloption1.value_only)
        self.assertFalse(self.shelloption2.value_only)
        self.assertTrue(self.shelloption3.value_only)
        self.assertFalse(self.shelloption4.value_only)

    def test_key_only(self):
        self.assertFalse(self.shelloption0.key_only)
        self.assertFalse(self.shelloption1.key_only)
        self.assertTrue(self.shelloption2.key_only)
        self.assertFalse(self.shelloption3.key_only)
        self.assertFalse(self.shelloption4.key_only)

    def test_output_as_list(self):
        self.assertEqual(self.shelloption0.output_as_list(), ['key', 'prefix=value'])
        self.assertEqual(self.shelloption1.output_as_list(), ['key', 'value'])
        self.assertEqual(self.shelloption2.output_as_list(), ['key'])
        self.assertEqual(self.shelloption3.output_as_list(), ['value'])
        self.assertFalse(self.shelloption4.output_as_list())


if __name__ == '__main__':
    unittest.main()