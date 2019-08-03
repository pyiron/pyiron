import unittest
import numpy as np
import os
from pyiron.atomistics.volumetric.generic import VolumetricData
from pyiron.atomistics.structure.atoms import Atoms
from pyiron.vasp.volumetric_data import VaspVolumetricData


"""
@author: surendralal

Unittests for the pyiron.objects.volumetric module
"""


class TestVolumetricData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        cls.data_dict = dict()
        cls.data_dict["cubic"] = np.ones((100, 100, 100))
        cls.data_dict["non_cubic"] = np.zeros((200, 50, 100))

    @classmethod
    def tearDownClass(cls):
        cls.execution_path = os.path.dirname(os.path.abspath(__file__))
        os.remove(os.path.join(cls.execution_path, "chgcar.cube"))
        os.remove(os.path.join(cls.execution_path, "random_CHGCAR"))

    def test_total_data_assertion(self):
        vd = VolumetricData()
        self.assertIsNone(vd._atoms)
        self.assertIsNone(vd.atoms)
        with self.assertRaises(ValueError):
            vd.total_data = list()
        with self.assertRaises(TypeError):
            vd.total_data = 1947.0
        for val in self.data_dict.values():
            vd.total_data = val
            self.assertIsInstance(vd.total_data, np.ndarray)

    def test_get_average_along_axis(self):
        vd = VolumetricData()
        for key, val in self.data_dict.items():
            vd.total_data = val
            for i in range(3):
                nz = vd.total_data.shape[i]
                if key == "cubic":
                    answer = np.ones(nz)
                    self.assertTrue(all(np.equal(answer, vd.get_average_along_axis(ind=i))))
                if key == "non-cubic":
                    answer = np.zeros(nz)
                    self.assertTrue(all(np.equal(answer, vd.get_average_along_axis(ind=i))))

    def test_write_cube(self):
        cd_obj = VaspVolumetricData()
        # file_name = os.path.abspath(os.path.join(self.execution_path, "..", "..", "static", "vasp_test_files",
        #                                          "CHGCAR_samples", "CHGCAR_no_spin"))
        # cd_obj.from_file(filename=file_name)
        # data_before = cd_obj.total_data.copy()
        # cd_obj.write_cube_file(filename=os.path.join(self.execution_path, "chgcar.cube"))
        # cd_obj.read_cube_file(filename=os.path.join(self.execution_path, "chgcar.cube"))
        # data_after = cd_obj.total_data.copy()
        # self.assertTrue(np.allclose(data_before, data_after))
        n_x, n_y, n_z = (3, 4, 2)
        random_array = np.random.rand(n_x, n_y, n_z)
        rd_obj = VolumetricData()
        rd_obj.atoms = Atoms("H2O", cell=np.eye(3)*10, positions=np.eye(3))
        rd_obj.total_data = random_array
        rd_obj.write_vasp_volumetric(filename=os.path.join(self.execution_path, "random_CHGCAR"))
        cd_obj.from_file(filename=os.path.join(self.execution_path, "random_CHGCAR"))
        self.assertTrue(np.allclose(cd_obj.total_data * cd_obj.atoms.get_volume(), rd_obj.total_data))
        data_before = cd_obj.total_data.copy()
        cd_obj.write_cube_file(filename=os.path.join(self.execution_path, "chgcar.cube"))
        cd_obj.read_cube_file(filename=os.path.join(self.execution_path, "chgcar.cube"))
        data_after = cd_obj.total_data.copy()
        self.assertTrue(np.allclose(data_before, data_after))
        file_name = os.path.abspath(os.path.join(self.execution_path, "..", "..", "static", "vasp_test_files",
                                                 "CHGCAR_samples", "CHGCAR_water"))
        cd_obj = VaspVolumetricData()
        cd_obj.from_file(file_name)
        data_before = cd_obj.total_data.copy()
        cd_obj.write_cube_file(filename=os.path.join(self.execution_path, "chgcar.cube"))
        cd_obj.read_cube_file(filename=os.path.join(self.execution_path, "chgcar.cube"))
        self.assertIsNotNone(cd_obj.atoms)
        data_after = cd_obj.total_data.copy()
        self.assertTrue(np.allclose(data_before, data_after))
        data_before = cd_obj.total_data.copy()
        cd_obj.write_vasp_volumetric(filename=os.path.join(self.execution_path, "random_CHGCAR"))
        cd_obj.from_file(filename=os.path.join(self.execution_path, "random_CHGCAR"), normalize=False)
        data_after = cd_obj.total_data.copy()
        self.assertTrue(np.allclose(data_before, data_after))


if __name__ == '__main__':
    unittest.main()
