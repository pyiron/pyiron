from copy import deepcopy
import pandas
from pyiron_base.objects.generic.parameters import GenericParameters
import unittest


class TestGenericParameters(unittest.TestCase):
    def setUp(self):
        self.generic_parameters_empty = GenericParameters()
        self.generic_parameters_str = GenericParameters()
        my_str = '''\
                par_1 1
                par_2 ab
                count 0
                write_restart True
                read_restart False'''
        self.generic_parameters_str.load_string(my_str)

    def test_load_string(self):
        self.assertEqual(self.generic_parameters_str.get("par_1"), 1)
        self.assertEqual(self.generic_parameters_str.get("par_2"), 'ab')
        self.assertEqual(self.generic_parameters_str.get("count"), 0)
        self.assertTrue(self.generic_parameters_str.get("write_restart"))
        self.assertFalse(self.generic_parameters_str.get("read_restart"))

    def test_get_pandas(self):
        self.assertEqual(str(self.generic_parameters_empty.get_pandas()),
                         str(pandas.DataFrame(columns=['Parameter', 'Value', 'Comment'])))

    def test_modify(self):
        self.assertEqual(self.generic_parameters_str.get("par_1"), 1)
        self.generic_parameters_str.modify(par_1=3)
        self.assertEqual(self.generic_parameters_str.get("par_1"), 3)

    def test_remove_keys(self):
        self.assertFalse(self.generic_parameters_str.get("read_restart"))
        data_frame_all_entries = deepcopy(self.generic_parameters_str)
        self.generic_parameters_str.remove_keys(["read_restart"])
        self.assertNotEqual(str(self.generic_parameters_str.get_pandas()), str(data_frame_all_entries.get_pandas()))
        self.generic_parameters_str.set(read_restart=False)
        self.assertFalse(self.generic_parameters_str.get("read_restart"))
        self.assertEqual(str(self.generic_parameters_str.get_pandas()), str(data_frame_all_entries.get_pandas()))


if __name__ == '__main__':
    unittest.main()