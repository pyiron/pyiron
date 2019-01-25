import os
from pyiron.base.generic.hdfio import FileHDFio
import unittest


class TestFileHDFio(unittest.TestCase):
    def setUp(self):
        self.current_dir = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
        self.empty_hdf5 = FileHDFio(file_name=self.current_dir + '/filehdfio_empty.h5')
        self.full_hdf5 = FileHDFio(file_name=self.current_dir + '/filehdfio_full.h5')
        self.es_hdf5 = FileHDFio(file_name=self.current_dir + "/../../static/dft/es_hdf.h5")

    def doCleanups(self):
        # os.remove('filehdfio_empty.h5')
        # os.remove('filehdfio_full.h5')
        pass

    def test_copy(self):
        # copy_empty = self.empty_hdf5.copy()
        # self.assertEqual(self.empty_hdf5, copy_empty)
        pass

    def test_file_name(self):
        print ('cdr: ', self.current_dir)
        print ('file: ', self.empty_hdf5.file_name)
        self.assertEqual(self.empty_hdf5.file_name, self.current_dir + '/filehdfio_empty.h5')
        self.assertEqual(self.full_hdf5.file_name, self.current_dir + '/filehdfio_full.h5')

    def test_h5_path(self):
        pass

    def test_open(self):
        pass

    def test_close(self):
        pass

    def test_remove_file(self):
        pass

    def test_get_from_table(self):
        pass

    def test_get_pandas(self):
        pass

    def test_get(self):
        pass

    def test_to_object(self):
        pass

    def test_put(self):
        pass

    def test_list_all(self):
        pass

    def test_list_nodes(self):
        pass

    def test_list_groups(self):
        pass

    def test_listdirs(self):
        pass

    def test_show_hdf(self):
        pass

    # def test_is_open(self):
    #     self.assertFalse(self.empty_hdf5.is_open)
    #     self.assertFalse(self.full_hdf5.is_open)
    #     # with self.empty_hdf5.open(h5_rel_path='test') as empty_test:
    #     #     self.assertTrue(empty_test.is_open)
    #     # with self.full_hdf5.open(h5_rel_path='.') as full_test:
    #     #     self.assertTrue(full_test.is_open)

    def test_is_empty(self):
        self.assertTrue(self.empty_hdf5.is_empty)
        self.assertTrue(self.full_hdf5.is_empty)

    def test_file_size(self):
        self.assertTrue(self.es_hdf5.file_size(self.es_hdf5)>0)

    def test_get_size(self):
        self.assertTrue(self.es_hdf5.get_size(self.es_hdf5)>0)

    def test_copy(self):
        self.assertIsInstance(self.es_hdf5.copy(), FileHDFio)

if __name__ == '__main__':
    unittest.main()

