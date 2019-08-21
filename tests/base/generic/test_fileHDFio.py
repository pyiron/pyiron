# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import numpy as np
from pyiron.base.generic.hdfio import FileHDFio
import unittest


class TestFileHDFio(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_dir = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
        cls.empty_hdf5 = FileHDFio(file_name=cls.current_dir + '/filehdfio_empty.h5')
        cls.full_hdf5 = FileHDFio(file_name=cls.current_dir + '/filehdfio_full.h5')
        cls.es_hdf5 = FileHDFio(file_name=cls.current_dir + "/../../static/dft/es_hdf.h5")
        with cls.full_hdf5.open('content') as hdf:
            hdf['array'] = np.array([1, 2, 3, 4, 5, 6])
            hdf['array_3d'] = np.array([[1, 2, 3], [4, 5, 6]])
            hdf['traj'] = np.array([[[1, 2, 3], [4, 5, 6]], [[7, 8, 9]]])
            hdf['dict'] = {'key_1': 1, 'key_2': 'hallo'}
            hdf['dict_numpy'] = {'key_1': 1, 'key_2': np.array([1, 2, 3, 4, 5, 6])}

    @classmethod
    def tearDownClass(cls):
        cls.current_dir = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
        os.remove(cls.current_dir + '/filehdfio_full.h5')

    def test_get_item(self):
        self.assertTrue(all(np.equal(self.full_hdf5['content/array'], np.array([1, 2, 3, 4, 5, 6]))))
        self.assertTrue(all(np.equal(self.full_hdf5['content']['array_3d'],
                                     np.array([[1, 2, 3], [4, 5, 6]])).flatten()))
        self.assertTrue(all(np.equal(self.full_hdf5['content/traj'][0],
                                     np.array([[1, 2, 3], [4, 5, 6]])).flatten()))
        self.assertTrue(all(np.equal(self.full_hdf5['content/traj'][1],
                                     np.array([[7, 8, 9]])).flatten()))
        self.assertEqual(self.full_hdf5['content/dict']['key_1'], 1)
        self.assertEqual(self.full_hdf5['content/dict']['key_2'], 'hallo')
        self.assertEqual(self.full_hdf5['content/dict_numpy']['key_1'], 1)
        self.assertTrue(all(np.equal(self.full_hdf5['content/dict_numpy']['key_2'], np.array([1, 2, 3, 4, 5, 6]))))

    def test_file_name(self):
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
        empty_file_dict = self.empty_hdf5.list_all()
        self.assertEqual(empty_file_dict['groups'], [])
        self.assertEqual(empty_file_dict['nodes'], [])
        es_file_dict = self.es_hdf5.list_all()
        self.assertEqual(es_file_dict['groups'], ['es_new', 'es_old'])
        self.assertEqual(es_file_dict['nodes'], [])
        es_group_dict = self.es_hdf5['es_new'].list_all()
        self.assertEqual(es_group_dict['groups'], ['dos'])
        self.assertEqual(es_group_dict['nodes'],
                         ['TYPE', 'efermi', 'eig_matrix', 'k_points', 'k_weights', 'occ_matrix'])

    def test_list_nodes(self):
        self.assertEqual(self.empty_hdf5.list_nodes(), [])
        self.assertEqual(self.es_hdf5['es_new'].list_nodes(),
                         ['TYPE', 'efermi', 'eig_matrix', 'k_points', 'k_weights', 'occ_matrix'])

    def test_list_groups(self):
        self.assertEqual(self.empty_hdf5.list_groups(), [])
        self.assertEqual(self.es_hdf5.list_groups(), ['es_new', 'es_old'])

    def test_listdirs(self):
        self.assertEqual(self.empty_hdf5.listdirs(), [])
        self.assertEqual(self.es_hdf5.listdirs(), ['es_new', 'es_old'])

    def test_show_hdf(self):
        pass

    def test_is_empty(self):
        self.assertTrue(self.empty_hdf5.is_empty)
        self.assertFalse(self.full_hdf5.is_empty)

    def test_file_size(self):
        self.assertTrue(self.es_hdf5.file_size(self.es_hdf5) > 0)

    def test_get_size(self):
        self.assertTrue(self.es_hdf5.get_size(self.es_hdf5) > 0)

    def test_copy(self):
        self.assertIsInstance(self.es_hdf5.copy(), FileHDFio)


if __name__ == '__main__':
    unittest.main()
