# import numpy as np
# from os import path as op
# from pandas import DataFrame, Series
# from pyiron.external.h5io import write_hdf5, read_hdf5, _TempDir, object_diff
# from scipy import sparse
# import unittest
#
#
# class TestHdf5(unittest.TestCase):
#     def test_hdf5(self):
#         """Test HDF5 IO
#         """
#         tempdir = _TempDir()
#         test_file = op.join(tempdir, 'test.hdf5')
#         sp = np.eye(3) if sparse is None else sparse.eye(3, 3, format='csc')
#         sp_csr = np.eye(3) if sparse is None else sparse.eye(3, 3, format='csr')
#         df = np.eye(3) if isinstance(DataFrame, type(None)) else DataFrame(
#             np.eye(3))
#         sr = np.eye(3) if isinstance(Series, type(None)) else Series(
#             np.random.randn(3))
#         sp[2, 2] = 2
#         sp_csr[2, 2] = 2
#         x = dict(a=dict(b=np.zeros(3)), c=np.zeros(2, np.complex128),
#                  d=[dict(e=(1, -2., 'hello', u'goodbyeu\u2764')), None], f=sp,
#                  g=dict(dfa=df, srb=sr), h=sp_csr)
#         write_hdf5(test_file, 1)
#         self.assertEqual(read_hdf5(test_file), 1)
#         self.assertRaises(IOError, write_hdf5, test_file, x)  # file exists
#         write_hdf5(test_file, x, overwrite=True)
#         self.assertRaises(IOError, read_hdf5, test_file + 'FOO')  # not found
#         xx = read_hdf5(test_file)
#         self.assertTrue(object_diff(x, xx) == '')  # no assert_equal, ugly output
#         write_hdf5(test_file, np.bool_(True), overwrite=True)
#         self.assertEqual(read_hdf5(test_file), np.bool_(True))
#
#         # bad title
#         self.assertRaises(ValueError, read_hdf5, test_file, title='nonexist')
#         self.assertRaises(ValueError, write_hdf5, test_file, x, overwrite=True,
#                       title=1)
#         self.assertRaises(ValueError, read_hdf5, test_file, title=1)
#         # unsupported objects
#         self.assertRaises(TypeError, write_hdf5, test_file, {1: 'foo'},
#                       overwrite=True)
#         self.assertRaises(TypeError, write_hdf5, test_file, object, overwrite=True)
#         # special_chars
#         spec_dict = {'first/second': 'third'}
#         self.assertRaises(ValueError, write_hdf5, test_file, spec_dict, overwrite=True)
#         self.assertRaises(ValueError, write_hdf5, test_file, spec_dict, overwrite=True,
#                       slash='brains')
#         write_hdf5(test_file, spec_dict, overwrite=True, slash='replace')
#         self.assertEqual(
#             read_hdf5(test_file, slash='replace').keys(), spec_dict.keys())
#         in_keys = list(read_hdf5(test_file, slash='ignore').keys())
#         self.assertTrue('{FWDSLASH}' in in_keys[0])
#         self.assertRaises(ValueError, read_hdf5, test_file, slash='brains')
#         # Testing that title slashes aren't replaced
#         write_hdf5(
#             test_file, spec_dict, title='one/two', overwrite=True, slash='replace')
#         self.assertEqual(read_hdf5(test_file, title='one/two', slash='replace').keys(),
#                          spec_dict.keys())
#
#         write_hdf5(test_file, 1, title='first', overwrite=True)
#         write_hdf5(test_file, 2, title='second', overwrite='update')
#         self.assertEqual(read_hdf5(test_file, title='first'), 1)
#         self.assertEqual(read_hdf5(test_file, title='second'), 2)
#         self.assertRaises(IOError, write_hdf5, test_file, 3, title='second')
#         write_hdf5(test_file, 3, title='second', overwrite='update')
#         self.assertEqual(read_hdf5(test_file, title='second'), 3)
#
#         write_hdf5(test_file, 5, title='second', overwrite='update', compression=5)
#         self.assertEqual(read_hdf5(test_file, title='second'), 5)
#
#
#     def test_path_support(self):
#         tempdir = _TempDir()
#         test_file = op.join(tempdir, 'test.hdf5')
#         write_hdf5(test_file, 1, title='first')
#         write_hdf5(test_file, 2, title='second/third', overwrite='update')
#         self.assertRaises(ValueError, read_hdf5, test_file, title='second')
#         self.assertEqual(read_hdf5(test_file, 'first'), 1)
#         self.assertEqual(read_hdf5(test_file, 'second/third'), 2)
#
#
#     def test_object_diff(self):
#         """Test object diff calculation
#         """
#         self.assertTrue('type' in object_diff(1, 1.))
#         self.assertTrue('missing' in object_diff({1: 1}, {}))
#         self.assertTrue('missing' in object_diff({}, {1: 1}))
#         self.assertTrue('length' in object_diff([], [1]))
#         self.assertTrue('value' in object_diff('a', 'b'))
#         self.assertTrue('None' in object_diff(None, 'b'))
#         self.assertTrue('array mismatch' in object_diff(np.array([1]), np.array([2])))
#         if sparse is not None:
#             a = sparse.coo_matrix([[1]])
#             b = sparse.coo_matrix([[1, 2]])
#             self.assertTrue('shape mismatch' in object_diff(a, b))
#             c = sparse.coo_matrix([[1, 1]])
#             self.assertTrue('1 element' in object_diff(b, c))
#         if not isinstance(DataFrame, type(None)):
#             for ob_type in (DataFrame, Series):
#                 a = ob_type([1])
#                 b = ob_type([1, 2])
#                 self.assertTrue('shape mismatch' in object_diff(a, b))
#                 c = ob_type([1, 3])
#                 self.assertTrue('1 element' in object_diff(b, c))
#         self.assertRaises(RuntimeError, object_diff, object, object)
#
#
# if __name__ == '__main__':
#     unittest.main()