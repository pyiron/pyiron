# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
import posixpath
import numpy as np
from pyiron.vasp.outcar import Outcar


class TestOutcar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.file_list = list()
        cls.outcar_parser = Outcar()
        # file_list = ["OUTCAR_1", "OUTCAR_2", "OUTCAR_3", "OUTCAR_4", "OUTCAR_5", "OUTCAR_6", "OUTCAR_7"]
        file_list = os.listdir(
            os.path.abspath(
                os.path.join(
                    cls.file_location, "../static/vasp_test_files/outcar_samples"
                )
            )
        )
        file_list = sorted(file_list)
        for f in file_list:
            direc = os.path.abspath(
                os.path.join(
                    cls.file_location, "../static/vasp_test_files/outcar_samples"
                )
            )
            filename = posixpath.join(direc, f)
            cls.file_list.append(filename)

    def test_from_file(self):
        for filename in self.file_list:
            self.outcar_parser.from_file(filename=filename)
            type_dict = dict()
            type_dict["energies"] = np.ndarray
            type_dict["scf_energies"] = list
            type_dict["forces"] = np.ndarray
            type_dict["positions"] = np.ndarray
            type_dict["cells"] = np.ndarray
            type_dict["steps"] = np.ndarray
            type_dict["temperatures"] = np.ndarray
            type_dict["time"] = np.ndarray
            type_dict["fermi_level"] = float
            type_dict["scf_dipole_moments"] = list
            type_dict["kin_energy_error"] = float
            type_dict["stresses"] = np.ndarray
            type_dict["irreducible_kpoints"] = tuple
            type_dict["pressures"] = np.ndarray
            type_dict["energies_int"] = np.ndarray
            type_dict["energies_zero"] = np.ndarray
            parse_keys = self.outcar_parser.parse_dict.keys()
            for key, value in type_dict.items():
                self.assertTrue(key in parse_keys)
                try:
                    self.assertIsInstance(self.outcar_parser.parse_dict[key], value)
                except AssertionError:
                    if int(filename.split("/OUTCAR_")[-1]) == 8:
                        self.assertEqual(key, "fermi_level")
                    else:
                        print(key, self.outcar_parser.parse_dict[key])
                        raise AssertionError("{} has the wrong type".format(key))

    def test_get_positions_and_forces(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_positions_and_forces(filename)
            self.assertIsInstance(output[0], np.ndarray)
            self.assertIsInstance(output[1], np.ndarray)
            self.assertEqual(len(output[0]), len(output[1]))
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                positions = np.array(
                    [
                        [
                            [0.0, 0.0, 0.0],
                            [1.85, 1.85, 0.0],
                            [1.85, 0.0, 1.85],
                            [0.0, 1.85, 1.85],
                        ]
                    ]
                )
                forces = np.array(
                    [
                        [
                            [-0.0, 0.0, 0.0],
                            [-0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                            [0.0, -0.0, -0.0],
                        ]
                    ]
                )
                self.assertEqual(positions.__str__(), output[0].__str__())
                self.assertEqual(forces.__str__(), output[1].__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 2:
                positions = np.array([[[0.0, 0.0, 0.0], [1.4, 1.4, 1.4]]])
                forces = np.array([[[-0.0, 0.0, 0.0], [0.0, -0.0, -0.0]]])
                self.assertEqual(positions.__str__(), output[0].__str__())
                self.assertEqual(forces.__str__(), output[1].__str__())
            if (
                int(filename.split("/OUTCAR_")[-1]) == 3
                or int(filename.split("/OUTCAR_")[-1]) == 5
                or int(filename.split("/OUTCAR_")[-1]) == 6
            ):
                positions = np.array([[[0.0, 0.0, 0.0], [1.4, 1.4, 1.4]]])
                forces = np.array([[[0.0, -0.0, 0.0], [-0.0, 0.0, -0.0]]])
                self.assertEqual(positions.__str__(), output[0].__str__())
                self.assertEqual(forces.__str__(), output[1].__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 4:
                positions = np.array([[[0.0, 0.0, 0.0], [1.4, 1.4, 1.4]]])
                forces = np.array([[[-0.0, -0.0, 0.0], [0.0, 0.0, 0.0]]])
                self.assertEqual(positions.__str__(), output[0].__str__())
                self.assertEqual(forces.__str__(), output[1].__str__())

    def test_get_positions(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_positions(filename)
            self.assertIsInstance(output, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                positions = np.array(
                    [
                        [
                            [0.0, 0.0, 0.0],
                            [1.85, 1.85, 0.0],
                            [1.85, 0.0, 1.85],
                            [0.0, 1.85, 1.85],
                        ]
                    ]
                )
                self.assertEqual(positions.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [2, 3, 4, 5, 6]:
                positions = np.array([[[0.0, 0.0, 0.0], [1.4, 1.4, 1.4]]])
                self.assertEqual(positions.__str__(), output.__str__())

    def test_get_forces(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_forces(filename)
            self.assertIsInstance(output, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                forces = np.array(
                    [
                        [
                            [-0.0, 0.0, 0.0],
                            [-0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                            [0.0, -0.0, -0.0],
                        ]
                    ]
                )
                self.assertEqual(forces.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 2:
                forces = np.array([[[-0.0, 0.0, 0.0], [0.0, -0.0, -0.0]]])
                self.assertEqual(forces.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [3, 5, 6]:
                forces = np.array([[[0.0, -0.0, 0.0], [-0.0, 0.0, -0.0]]])
                self.assertEqual(forces.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 4:
                forces = np.array([[[-0.0, -0.0, 0.0], [0.0, 0.0, 0.0]]])
                self.assertEqual(forces.__str__(), output.__str__())

    def test_get_cells(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_cells(filename)
            self.assertIsInstance(output, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                cells = np.array([[[3.7, 0.0, 0.0], [0.0, 3.7, 0.0], [0.0, 0.0, 3.7]]])
                self.assertEqual(cells.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [2, 4, 5, 6]:
                cells = np.array([[[2.8, 0.0, 0.0], [0.0, 2.8, 0.0], [0.0, 0.0, 2.8]]])
                self.assertEqual(cells.__str__(), output.__str__())

    def test_get_total_energies(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_total_energies(filename)
            self.assertIsInstance(output, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                energies = np.array([-24.59765979])
                self.assertEqual(energies.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 2:
                energies = np.array([-17.73798679])
                self.assertEqual(energies.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [3, 4]:
                energies = np.array([-18.40140804])
                self.assertEqual(energies.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [5, 6]:
                energies = np.array([-18.40218607])
                self.assertEqual(energies.__str__(), output.__str__())

    def test_get_all_total_energies(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_all_total_energies(filename)
            self.assertIsInstance(output, list)
            self.assertIsInstance(output[-1], np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                energies = np.array(
                    [
                        164.50346113,
                        -22.34826105,
                        -29.05044151,
                        -29.09407035,
                        -29.0944954,
                        -27.08800718,
                        -24.63733705,
                        -24.60427656,
                        -24.59776189,
                        -24.59765131,
                        -24.59765979,
                    ]
                )
                self.assertEqual(energies.__str__(), output[0].__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 2:
                energies = np.array(
                    [
                        49.31292902,
                        -18.6411937,
                        -20.8983855,
                        -20.92056153,
                        -20.92073722,
                        -17.95357529,
                        -17.76756585,
                        -17.739524,
                        -17.7379502,
                        -17.73798679,
                    ]
                )
                self.assertEqual(energies.__str__(), output[0].__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [3, 4]:
                energies = np.array(
                    [
                        44.19769573,
                        -19.71777011,
                        -20.86598567,
                        -20.87149728,
                        -20.87156762,
                        -18.51206993,
                        -18.43794441,
                        -18.41311659,
                        -18.40310749,
                        -18.4014402,
                        -18.40140804,
                    ]
                )
                self.assertEqual(energies.__str__(), output[0].__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [5, 6]:
                energies = np.array(
                    [
                        69.15405871,
                        -18.00103756,
                        -20.85032272,
                        -20.87198829,
                        -20.87222361,
                        -18.51305208,
                        -18.43873509,
                        -18.41390999,
                        -18.40387841,
                        -18.40222137,
                        -18.40218607,
                    ]
                )
                self.assertEqual(energies.__str__(), output[0].__str__())

    def test_get_energy_without_entropy(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_energy_without_entropy(filename)
            energy_tot = self.outcar_parser.get_total_energies(filename)
            self.assertIsInstance(output, np.ndarray)
            self.assertEqual(len(energy_tot), len(output))

    def test_get_energy_sigma_0(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_energy_sigma_0(filename)
            energy_tot = self.outcar_parser.get_total_energies(filename)
            self.assertIsInstance(output, np.ndarray)
            self.assertEqual(len(energy_tot), len(output))

    def test_get_temperatures(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_temperatures(filename)
            self.assertIsInstance(output, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) in [1, 2, 3, 4, 5, 6]:
                temperatures = np.array([0.0])
                self.assertEqual(temperatures.__str__(), output.__str__())

    def test_get_steps(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_steps(filename)
            self.assertIsInstance(output, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) in [1, 2, 3, 4, 5, 6]:
                steps = np.array(
                    [
                        0.0,
                        0.02040816,
                        0.04081633,
                        0.06122449,
                        0.08163265,
                        0.10204082,
                        0.12244898,
                        0.14285714,
                        0.16326531,
                        0.18367347,
                        0.20408163,
                        0.2244898,
                        0.24489796,
                        0.26530612,
                        0.28571429,
                        0.30612245,
                        0.32653061,
                        0.34693878,
                        0.36734694,
                        0.3877551,
                        0.40816327,
                        0.42857143,
                        0.44897959,
                        0.46938776,
                        0.48979592,
                        0.51020408,
                        0.53061224,
                        0.55102041,
                        0.57142857,
                        0.59183673,
                        0.6122449,
                        0.63265306,
                        0.65306122,
                        0.67346939,
                        0.69387755,
                        0.71428571,
                        0.73469388,
                        0.75510204,
                        0.7755102,
                        0.79591837,
                        0.81632653,
                        0.83673469,
                        0.85714286,
                        0.87755102,
                        0.89795918,
                        0.91836735,
                        0.93877551,
                        0.95918367,
                        0.97959184,
                        1.0,
                    ]
                )
                self.assertEqual(steps.__str__(), output.__str__())

    def test_get_time(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_time(filename)
            self.assertIsInstance(output, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) in [1, 2, 3, 4, 5, 6]:
                time = np.array(
                    [
                        0.0,
                        0.02040816,
                        0.04081633,
                        0.06122449,
                        0.08163265,
                        0.10204082,
                        0.12244898,
                        0.14285714,
                        0.16326531,
                        0.18367347,
                        0.20408163,
                        0.2244898,
                        0.24489796,
                        0.26530612,
                        0.28571429,
                        0.30612245,
                        0.32653061,
                        0.34693878,
                        0.36734694,
                        0.3877551,
                        0.40816327,
                        0.42857143,
                        0.44897959,
                        0.46938776,
                        0.48979592,
                        0.51020408,
                        0.53061224,
                        0.55102041,
                        0.57142857,
                        0.59183673,
                        0.6122449,
                        0.63265306,
                        0.65306122,
                        0.67346939,
                        0.69387755,
                        0.71428571,
                        0.73469388,
                        0.75510204,
                        0.7755102,
                        0.79591837,
                        0.81632653,
                        0.83673469,
                        0.85714286,
                        0.87755102,
                        0.89795918,
                        0.91836735,
                        0.93877551,
                        0.95918367,
                        0.97959184,
                        1.0,
                    ]
                )
                self.assertEqual(time.__str__(), output.__str__())

    def test_get_fermi_level(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_fermi_level(filename)
            try:
                self.assertIsInstance(output, float)
            except AssertionError:
                self.assertEqual(int(filename.split("/OUTCAR_")[-1]), 8)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                fermi_level = 2.9738
                self.assertEqual(fermi_level, output)
            if int(filename.split("/OUTCAR_")[-1]) == 2:
                fermi_level = 5.9788
                self.assertEqual(fermi_level, output)
            if int(filename.split("/OUTCAR_")[-1]) in [3, 4]:
                fermi_level = 5.9613
                self.assertEqual(fermi_level, output)
            if int(filename.split("/OUTCAR_")[-1]) in [5, 6]:
                fermi_level = 5.9614
                self.assertEqual(fermi_level, output)

    def test_get_magnetization(self):
        for filename in self.file_list:
            output, final_magmoms = self.outcar_parser.get_magnetization(filename)
            positions = self.outcar_parser.get_positions(filename)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                magnetization = [np.array([])]
                self.assertEqual(magnetization.__str__(), output.__str__())
                self.assertFalse(final_magmoms)
            if int(filename.split("/OUTCAR_")[-1]) == 2:
                magnetization = [np.array([])]
                self.assertEqual(magnetization.__str__(), output.__str__())
                self.assertFalse(final_magmoms)
            if int(filename.split("/OUTCAR_")[-1]) == 3:
                magnetization = [
                    np.array(
                        [
                            4.0,
                            4.0,
                            4.0,
                            4.0,
                            3.1349091,
                            4.1157043,
                            4.2255076,
                            4.2514329,
                            4.2041297,
                            4.2125079,
                            4.2125079,
                        ]
                    )
                ]
                self.assertEqual(magnetization.__str__(), output.__str__())
                self.assertFalse(final_magmoms)
            if int(filename.split("/OUTCAR_")[-1]) == 4:
                magnetization = [
                    np.array(
                        [
                            4.0,
                            4.0,
                            4.0,
                            4.0,
                            3.1349091,
                            4.1157043,
                            4.2255076,
                            4.2514329,
                            4.2041297,
                            4.2125079,
                            4.2125079,
                        ]
                    )
                ]
                final_mag_lst = [[2.111, 2.111]]
                self.assertEqual(magnetization.__str__(), output.__str__())
                self.assertEqual(final_magmoms, final_mag_lst)
            if int(filename.split("/OUTCAR_")[-1]) == 5:
                magnetization = [
                    np.array(
                        [
                            [0.0, 4.0, 0.0],
                            [0.0, 4.0, 0.0],
                            [0.0, 4.0, 0.0],
                            [0.0, 4.0, 0.0],
                            [-0.0, 3.1360688, -0.0],
                            [-0.0, 4.1168471, 0.0],
                            [0.0, 4.2261074, -0.0],
                            [-0.0, 4.2517361, 0.0],
                            [-0.0, 4.2043676, 0.0],
                            [0.0, 4.2128061, 0.0],
                            [0.0, 4.2128061, 0.0],
                        ]
                    )
                ]
                self.assertEqual(magnetization.__str__(), output.__str__())
                self.assertFalse(final_magmoms)

            if int(filename.split("/OUTCAR_")[-1]) == 6:
                magnetization = [
                    np.array(
                        [
                            [0.0, 4.0, 0.0],
                            [0.0, 4.0, 0.0],
                            [0.0, 4.0, 0.0],
                            [0.0, 4.0, 0.0],
                            [-0.0, 3.1360688, -0.0],
                            [-0.0, 4.1168471, 0.0],
                            [0.0, 4.2261074, -0.0],
                            [-0.0, 4.2517361, 0.0],
                            [-0.0, 4.2043676, 0.0],
                            [0.0, 4.2128061, 0.0],
                            [0.0, 4.2128061, 0.0],
                        ]
                    )
                ]
                final_mag_lst = [[[0.0, 2.111, -0.0], [0.0, 2.111, 0.0]]]
                self.assertEqual(magnetization.__str__(), output.__str__())
                self.assertEqual(final_magmoms, final_mag_lst)
                final_magmoms = np.array(final_magmoms)
                final_magmoms[
                    :, np.arange(len(positions[0]), dtype=int), :
                ] = final_magmoms.copy()
                self.assertEqual(np.array(final_magmoms).shape[1], positions.shape[1])

            if int(filename.split("/OUTCAR_")[-1]) == 7:
                self.assertEqual((1, 49, 3), np.array(output).shape)
                self.assertEqual((11, 32, 3), np.array(final_magmoms).shape)
                final_magmoms = np.array(final_magmoms)
                final_magmoms[
                    :, np.arange(len(positions[0]), dtype=int), :
                ] = final_magmoms.copy()
                self.assertEqual(np.array(final_magmoms).shape[1], positions.shape[1])

    def test_get_broyden_mixing_mesh(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_broyden_mixing_mesh(filename)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                mixing = 729
                self.assertEqual(mixing, output)
            if int(filename.split("/OUTCAR_")[-1]) in [2, 3, 4, 5, 6]:
                mixing = 343
                self.assertEqual(mixing, output)
            if int(filename.split("/OUTCAR_")[-1]) == 8:
                mixing = 21 ** 3
                self.assertEqual(mixing, output)

    def test_get_dipole_moments(self):
        for filename in self.file_list:
            output = self.outcar_parser.get_dipole_moments(filename)
            self.assertTrue(output)
            self.assertIsInstance(output, list)
            if len(output[0]) > 1:
                self.assertIsInstance(output[-1], np.ndarray)
                self.assertIsInstance(output[-1][-1], np.ndarray)
                self.assertEqual(len(output[-1][-1]), 3)

    def test_get_stresses(self):
        for filename in self.file_list:
            output_si = self.outcar_parser.get_stresses(filename, si_unit=True)
            output_kb = self.outcar_parser.get_stresses(filename, si_unit=False)
            self.assertIsInstance(output_si, np.ndarray)
            self.assertIsInstance(output_kb, np.ndarray)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                pullay_si = np.array([[-14.41433, -14.41433, -14.41433, 0.0, 0.0, 0.0]])
                pullay_kb = np.array(
                    [[-455.93181, -455.93181, -455.93181, 0.0, 0.0, 0.0]]
                )
                self.assertEqual(output_si.__str__(), pullay_si.__str__())
                self.assertEqual(output_kb.__str__(), pullay_kb.__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 2:
                pullay_si = np.array([[-5.38507, -5.38507, -5.38507, -0.0, 0.0, -0.0]])
                pullay_kb = np.array([[-393.032, -393.032, -393.032, -0.0, 0.0, -0.0]])
                self.assertEqual(output_si.__str__(), pullay_si.__str__())
                self.assertEqual(output_kb.__str__(), pullay_kb.__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 3:
                pullay_si = np.array([[-3.29441, -3.29441, -3.29441, 0.0, -0.0, 0.0]])
                pullay_kb = np.array(
                    [[-240.44384, -240.44384, -240.44384, 0.0, -0.0, 0.0]]
                )
                self.assertEqual(output_si.__str__(), pullay_si.__str__())
                self.assertEqual(output_kb.__str__(), pullay_kb.__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 4:
                pullay_si = np.array([[-3.29441, -3.29441, -3.29441, -0.0, -0.0, -0.0]])
                pullay_kb = np.array(
                    [[-240.44384, -240.44384, -240.44384, -0.0, -0.0, -0.0]]
                )
                self.assertEqual(output_si.__str__(), pullay_si.__str__())
                self.assertEqual(output_kb.__str__(), pullay_kb.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [5, 6]:
                pullay_si = np.array([[-3.3066, -3.3066, -3.3066, 0.0, -0.0, 0.0]])
                pullay_kb = np.array(
                    [[-241.33405, -241.33409, -241.33405, 0.0, -0.0, 0.0]]
                )
                self.assertEqual(output_si.__str__(), pullay_si.__str__())
                self.assertEqual(output_kb.__str__(), pullay_kb.__str__())

    def test_get_kinetic_energy_error(self):
        for filename in self.file_list:
            output_total = self.outcar_parser.get_kinetic_energy_error(filename)
            self.assertIsInstance(output_total, float)
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                kinetic_energy_error = 0.0774
                self.assertEqual(4 * kinetic_energy_error, output_total)
            if int(filename.split("/OUTCAR_")[-1]) in [2, 3, 4, 5, 6]:
                kinetic_energy_error = 0.0664
                self.assertEqual(2 * kinetic_energy_error, output_total)
            if int(filename.split("/OUTCAR_")[-1]) == 7:
                kinetic_energy_error = 16.0 * 0.0088 + 16.0 * 0.002
                self.assertEqual(kinetic_energy_error, output_total)
            if int(filename.split("/OUTCAR_")[-1]) == 8:
                kinetic_energy_error = 1.0 * 0.0236 + 32.0 * 0.0091
                self.assertEqual(kinetic_energy_error, output_total)

    def test_get_kpoints_irreducible_reciprocal(self):
        for filename in self.file_list:
            output_all = self.outcar_parser.get_irreducible_kpoints(
                filename, reciprocal=True, weight=True, planewaves=True
            )
            output_w = self.outcar_parser.get_irreducible_kpoints(
                filename, reciprocal=True, weight=True, planewaves=False
            )
            output_p = self.outcar_parser.get_irreducible_kpoints(
                filename, reciprocal=True, weight=False, planewaves=True
            )
            output_k = self.outcar_parser.get_irreducible_kpoints(
                filename, reciprocal=True, weight=False, planewaves=False
            )
            self.assertIsInstance(output_all[0], np.ndarray)
            self.assertIsInstance(output_all[1], np.ndarray)
            self.assertIsInstance(output_all[2], np.ndarray)
            self.assertIsInstance(output_w[0], np.ndarray)
            self.assertIsInstance(output_w[1], np.ndarray)
            self.assertIsInstance(output_p[0], np.ndarray)
            self.assertIsInstance(output_p[1], np.ndarray)
            self.assertIsInstance(output_k, np.ndarray)
            self.assertEqual(output_all[0].__str__(), output_w[0].__str__())
            self.assertEqual(output_all[1].__str__(), output_w[1].__str__())
            self.assertEqual(output_all[0].__str__(), output_p[0].__str__())
            self.assertEqual(output_all[2].__str__(), output_p[1].__str__())
            self.assertEqual(output_all[0].__str__(), output_k.__str__())
            if int(filename.split("/OUTCAR_")[-1]) == 1:
                output = (
                    np.array(
                        [
                            [0.125, 0.125, 0.125],
                            [0.375, 0.125, 0.125],
                            [0.375, 0.375, 0.125],
                            [0.375, 0.375, 0.375],
                        ]
                    ),
                    np.array([8.0, 24.0, 24.0, 8.0]),
                    np.array([452.0, 457.0, 451.0, 459.0]),
                )
                self.assertEqual(output_all.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [2, 3, 4]:
                output = (
                    np.array(
                        [
                            [0.125, 0.125, 0.125],
                            [0.375, 0.125, 0.125],
                            [0.375, 0.375, 0.125],
                            [0.375, 0.375, 0.375],
                        ]
                    ),
                    np.array([8.0, 24.0, 24.0, 8.0]),
                    np.array([196.0, 199.0, 196.0, 190.0]),
                )
                self.assertEqual(output_all.__str__(), output.__str__())
            if int(filename.split("/OUTCAR_")[-1]) in [5, 6]:
                output = (
                    np.array(
                        [
                            [0.125, 0.125, 0.125],
                            [0.375, 0.125, 0.125],
                            [-0.375, 0.125, 0.125],
                            [0.125, 0.375, 0.125],
                            [0.375, 0.375, 0.125],
                            [-0.375, 0.375, 0.125],
                            [0.375, 0.125, 0.375],
                            [0.375, 0.375, 0.375],
                        ]
                    ),
                    np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]),
                    np.array([392.0, 398.0, 398.0, 398.0, 392.0, 392.0, 392.0, 380.0]),
                )
                self.assertEqual(output_all.__str__(), output.__str__())

    def test_get_nelect(self):
        n_elect_list = [40.0, 16.0, 16.0, 16.0, 16.0, 16.0, 224.0, 358.0, 8]
        for filename in self.file_list:
            i = int(filename.split("_")[-1]) - 1
            self.assertEqual(n_elect_list[i], self.outcar_parser.get_nelect(filename))

    def test_get_band_properties(self):
        fermi_level_list = [2.9738, 5.9788, 5.9613, 5.9613, 5.9614, 5.9614, 3.9092]
        vbm_level_list = [3.076, 6.5823, 6.6175, 6.6175, 6.618, 6.618]
        cbm_level_list = [4.3112, 6.7396, 6.8568, 6.8568, 6.8569, 6.8569]
        for i, filename in enumerate(self.file_list):
            fermi_list, vbm_list, cbm_list = \
                self.outcar_parser.get_band_properties(filename=filename)
            if i <= 5:
                self.assertEqual(fermi_list[-1], fermi_level_list[i])
                self.assertEqual(vbm_list[-1], vbm_level_list[i])
                self.assertEqual(cbm_list[-1], cbm_level_list[i])
            elif i == 8:
                self.assertTrue(np.array_equal(fermi_list, [-7.3532, -5.7586]))
                self.assertTrue(np.array_equal(vbm_list, [-9.7778, -7.9449]))
                self.assertTrue(np.array_equal(cbm_list, [-1.8327, -3.5781]))


if __name__ == "__main__":
    unittest.main()
