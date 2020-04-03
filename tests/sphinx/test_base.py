# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import numpy as np
import unittest
import warnings
import scipy
import scipy.constants
from pyiron.project import Project
from pyiron.atomistics.structure.periodic_table import PeriodicTable
from pyiron.atomistics.structure.atoms import Atoms
from collections import OrderedDict as odict

BOHR_TO_ANGSTROM = (
        scipy.constants.physical_constants["Bohr radius"][0] / scipy.constants.angstrom
)
HARTREE_TO_EV = scipy.constants.physical_constants["Hartree energy in eV"][0]
HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM = HARTREE_TO_EV / BOHR_TO_ANGSTROM


class TestSphinx(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "../static/sphinx"))
        pt = PeriodicTable()
        pt.add_element(parent_element="Fe", new_element="Fe_up", spin="0.5")
        Fe_up = pt.element("Fe_up")
        cls.basis = Atoms(
            elements=[Fe_up, Fe_up],
            scaled_positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            cell=2.6 * np.eye(3),
        )
        cls.sphinx = cls.project.create_job("Sphinx", "job_sphinx")
        cls.sphinx_bs = cls.project.create_job("Sphinx", "sphinx_test_bs")
        cls.sphinx_2_3 = cls.project.create_job("Sphinx", "sphinx_test_2_3")
        cls.sphinx_2_5 = cls.project.create_job("Sphinx", "sphinx_test_2_5")
        cls.sphinx_aborted = cls.project.create_job("Sphinx", "sphinx_test_aborted")
        cls.sphinx.structure = cls.basis
        cls.sphinx_bs.structure = cls.project.create_structure("Fe", "bcc", 2.81)
        cls.sphinx_bs.structure = cls.sphinx_bs.structure.create_line_mode_structure()
        cls.sphinx_2_3.structure = Atoms(
            elements=["Fe", "Fe"],
            scaled_positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            cell=2.6 * np.eye(3),
        )
        cls.sphinx_2_5.structure = Atoms(
            elements=["Fe", "Ni"],
            scaled_positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            cell=2.83 * np.eye(3),
        )
        cls.sphinx_aborted.structure = Atoms(
            elements=32 * ["Fe"],
            scaled_positions=np.arange(32 * 3).reshape(-1, 3) / (32 * 3),
            cell=3.5 * np.eye(3),
        )
        cls.sphinx_aborted.status.aborted = True
        cls.current_dir = os.path.abspath(os.getcwd())
        cls.sphinx._create_working_directory()
        cls.sphinx_2_3._create_working_directory()
        cls.sphinx.input["VaspPot"] = False
        cls.sphinx.write_input()
        cls.sphinx.version = "2.6"
        cls.sphinx_2_3.to_hdf()
        cls.sphinx_2_3.decompress()
        cls.sphinx_2_5.decompress()

    @classmethod
    def tearDownClass(cls):
        cls.sphinx_2_3.decompress()
        cls.sphinx_2_5.decompress()
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/basis.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/control.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/guess.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/hamilton.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/input.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/potentials.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/structure.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/userparameters.sx",
            )
        )
        os.remove(
            os.path.join(
                cls.file_location,
                "../static/sphinx/job_sphinx_hdf5/job_sphinx/Fe_GGA.atomicdata",
            )
        )
        os.rmdir(
            os.path.join(
                cls.file_location, "../static/sphinx/job_sphinx_hdf5/job_sphinx"
            )
        )
        os.rmdir(os.path.join(cls.file_location, "../static/sphinx/job_sphinx_hdf5"))
        os.remove(
            os.path.join(cls.file_location, "../static/sphinx/sphinx_test_2_3.h5")
        )

    def test_write_basis(self):
        file_content = [
            "eCut = EnCut/13.606;\n",
            "kPoint {\n",
            "\tcoords = KpointCoords;\n",
            "\tweight = 1;\n",
            "\trelative;\n",
            "}\n",
            "folding = KpointFolding;\n",
            "saveMemory;\n",
        ]
        with open(
                os.path.join(
                    self.file_location,
                    "../static/sphinx/job_sphinx_hdf5/job_sphinx/basis.sx",
                )
        ) as basis_sx:
            lines = basis_sx.readlines()
        self.assertEqual(file_content, lines)

    def test_id_pyi_to_spx(self):
        self.assertEqual(len(self.sphinx.id_pyi_to_spx), len(self.sphinx.structure))
        self.assertEqual(len(self.sphinx.id_spx_to_pyi), len(self.sphinx.structure))

    def test_write_control(self):
        file_content = [
            "scfDiag {\n",
            "\trhoMixing = 1.0;\n",
            "\tspinMixing = 1.0;\n",
            "\tdEnergy = Ediff/" + str(scipy.constants.physical_constants["Hartree energy in eV"][0]) + ";\n",
            "\tmaxSteps = 400;\n",
            "\tblockCCG {}\n",
            "}\n",
            "evalForces {\n",
            '\tfile = "relaxHist.sx";\n',
            "}\n",
        ]
        file_name = os.path.join(
            self.file_location, "../static/sphinx/job_sphinx_hdf5/job_sphinx/control.sx"
        )
        with open(file_name) as control_sx:
            lines = control_sx.readlines()
        self.assertEqual(file_content, lines)

    def test_write_input(self):
        file_content = [
            "//job_sphinx;\n",
            "//SPHinX input file generated by pyiron;\n",
            "format paw;\n",
            "include <parameters.sx>;\n",
            "include <userparameters.sx>;\n",
            "pawPot {\n",
            "\tinclude <potentials.sx>;\n",
            "}\n",
            "structure {\n",
            "\tinclude <structure.sx>;\n",
            "}\n",
            "basis {\n",
            "\tinclude <basis.sx>;\n",
            "}\n",
            "PAWHamiltonian {\n",
            "\tinclude <hamilton.sx>;\n",
            "}\n",
            "initialGuess {\n",
            "\tinclude <guess.sx>;\n",
            "}\n",
            "main {\n",
            "\tinclude <control.sx>;\n",
            "}\n",
        ]
        file_name = os.path.join(
            self.file_location, "../static/sphinx/job_sphinx_hdf5/job_sphinx/input.sx"
        )
        with open(file_name) as input_sx:
            lines = input_sx.readlines()
        self.assertEqual(file_content, lines)

    def test_write_userparameters(self):
        file_content = [
            "EnCut=340;\n",
            "KpointCoords=[0.5, 0.5, 0.5];\n",
            "KpointFolding=[4, 4, 4];\n",
            "EmptyStates=6;\n",
            "Sigma=0.2;\n",
            "Xcorr=PBE;\n",
            "VaspPot=false;\n",
            "Estep=400;\n",
            "Ediff=0.0001;\n",
            "WriteWaves=true;\n",
            "KJxc=false;\n",
            "SaveMemory=true;\n",
            "CoarseRun=false;\n",
            "rhoMixing=1.0;\n",
            "spinMixing=1.0;\n",
            "CheckOverlap=true;\n",
            "THREADS=1;\n",
        ]
        file_name = os.path.join(
            self.file_location,
            "../static/sphinx/job_sphinx_hdf5/job_sphinx/userparameters.sx",
        )
        with open(file_name) as userparameters_sx:
            lines = userparameters_sx.readlines()
        self.assertEqual(file_content, lines)

    def test_plane_wave_cutoff(self):
        with self.assertRaises(ValueError):
            self.sphinx.plane_wave_cutoff = -1
        self.sphinx.plane_wave_cutoff = 340
        self.assertEqual(self.sphinx.plane_wave_cutoff, 340)

    def test_set_kpoints(self):

        mesh = [2, 3, 4]
        center_shift = [0.1, 0.1, 0.1]

        trace = {"my_path": [("GAMMA", "H"), ("H", "N"), ("P", "H")]}
        kpoints_dict = odict([('kPoints',
                               odict([('relative', None),
                                      ('from',
                                       odict([('coords', '[0.0, 0.0, 0.0]'),
                                              ('label', '"GAMMA"')])),
                                      ('to___0',
                                       odict([('coords', '[0.5, -0.5, 0.5]'),
                                              ('nPoints', 20),
                                              ('label', '"H"')])),
                                      ('to___1',
                                       odict([('coords', '[0.0, 0.0, 0.5]'),
                                              ('nPoints', 20),
                                              ('label', '"N"')])),
                                      ('to___1___1',
                                       odict([('coords', '[0.25, 0.25, 0.25]'),
                                              ('nPoints', 0),
                                              ('label', '"P"')])),
                                      ('to___2',
                                       odict([('coords', '[0.5, -0.5, 0.5]'),
                                              ('nPoints', 20),
                                              ('label', '"H"')]))]))])

        with self.assertRaises(ValueError):
            self.sphinx_bs.set_kpoints(symmetry_reduction="pyiron rules!")
        with self.assertRaises(ValueError):
            self.sphinx_bs.set_kpoints(scheme="no valid scheme")
        with self.assertRaises(ValueError):
            self.sphinx_bs.set_kpoints(scheme="Line", path_name="my_path")

        self.sphinx_bs.structure.add_high_symmetry_path(trace)
        with self.assertRaises(ValueError):
            self.sphinx_bs.set_kpoints(scheme="Line", n_path=20)
        with self.assertRaises(AssertionError):
            self.sphinx_bs.set_kpoints(scheme="Line", path_name="worng name", n_path=20)

        self.sphinx_bs.set_kpoints(scheme="Line", path_name="my_path", n_path=20)
        self.assertEqual(kpoints_dict, self.sphinx_bs._kpoints_odict)

        self.sphinx_bs.set_kpoints(scheme="MP", mesh=mesh, center_shift=center_shift)
        self.assertIsNone(self.sphinx_bs._kpoints_odict)
        self.assertEqual(self.sphinx_bs.input["KpointFolding"], mesh)
        self.assertEqual(self.sphinx_bs.input["KpointCoords"], center_shift)

    def test_set_empty_states(self):
        with self.assertRaises(ValueError):
            self.sphinx.set_empty_states(-1)
        self.sphinx.set_empty_states(666)
        self.assertEqual(self.sphinx.input["EmptyStates"], 666)
        self.sphinx.set_empty_states()
        self.assertEqual(self.sphinx.input["EmptyStates"], "auto")

    def test_write_guess(self):
        file_content = [
            "waves {\n",
            "\tlcao {}\n",
            "\tpawBasis;\n",
            "}\n",
            "rho {\n",
            "\tatomicOrbitals;\n",
            "\tatomicSpin {\n",
            '\t\tlabel = "spin_0.5";\n',
            "\t\tspin = 0.5;\n",
            "\t}\n",
            "\tatomicSpin {\n",
            '\t\tlabel = "spin_0.5";\n',
            "\t\tspin = 0.5;\n",
            "\t}\n",
            "}\n",
        ]
        file_name = os.path.join(
            self.file_location, "../static/sphinx/job_sphinx_hdf5/job_sphinx/guess.sx"
        )
        with open(file_name) as guess_sx:
            lines = guess_sx.readlines()
        self.assertEqual(file_content, lines)

    def test_write_potentials(self):
        file_content = [
            "species {\n",
            '\tname = "Fe";\n',
            '\tpotType = "AtomPAW";\n',
            '\telement = "Fe";\n',
            '\tpotential = "Fe_GGA.atomicdata";\n',
            "}\n",
        ]
        file_name = os.path.join(
            self.file_location,
            "../static/sphinx/job_sphinx_hdf5/job_sphinx/potentials.sx",
        )
        with open(file_name) as potentials_sx:
            lines = potentials_sx.readlines()
        self.assertEqual(file_content, lines)

    def test_fix_spin_constraint(self):
        self.assertIsNone(self.sphinx.fix_spin_constraint)
        with self.assertRaises(ValueError):
            self.sphinx.fix_spin_constraint = 3
        self.sphinx.fix_spin_constraint = False
        self.assertIsInstance(self.sphinx.fix_spin_constraint, bool)

    def test_calc_static(self):
        self.sphinx.calc_static(algorithm="wrong_algorithm")
        self.assertFalse(
            "keepRho"
            in self.sphinx.input_writer._odict_to_spx_input(self.sphinx._control_str)
        )
        self.assertTrue(
            "blockCCG"
            in self.sphinx.input_writer._odict_to_spx_input(self.sphinx._control_str)
        )
        self.sphinx.restart_file_list.append("randomfile")
        self.sphinx.calc_static(algorithm="ccg")
        self.assertTrue(
            "keepRho"
            in self.sphinx.input_writer._odict_to_spx_input(self.sphinx._control_str)
        )
        self.assertEqual(self.sphinx.input["Estep"], 400)
        self.assertTrue(
            "CCG"
            in self.sphinx.input_writer._odict_to_spx_input(self.sphinx._control_str)
        )

    def test_calc_minimize(self):
        self.sphinx.calc_minimize(electronic_steps=100, ionic_steps=50)
        self.assertEqual(self.sphinx.input["Estep"], 100)
        self.assertEqual(self.sphinx.input["Istep"], 50)
        self.assertEqual(self.sphinx._control_str['linQN']['maxSteps'], '50')

    def test_check_setup(self):
        self.assertFalse(self.sphinx.check_setup())

    def test_validate_ready_to_run(self):
        self.sphinx.validate_ready_to_run()

    def test_set_mixing_parameters(self):
        self.assertRaises(
            AssertionError, self.sphinx.set_mixing_parameters, "LDA", 7, 1.0, 1.0
        )
        self.assertRaises(
            AssertionError, self.sphinx.set_mixing_parameters, "PULAY", 1.2, 1.0, 1.0
        )
        self.assertRaises(
            ValueError, self.sphinx.set_mixing_parameters, "PULAY", 7, -0.1, 1.0
        )
        self.assertRaises(
            ValueError, self.sphinx.set_mixing_parameters, "PULAY", 7, 1.0, 2.0
        )
        self.sphinx.set_mixing_parameters("PULAY", 7, 0.5, 0.2)
        self.assertEqual(self.sphinx.input["rhoMixing"], 0.5)
        self.assertEqual(self.sphinx.input["spinMixing"], 0.2)

    def test_exchange_correlation_functional(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.sphinx.exchange_correlation_functional = "llda"
            self.assertEqual(len(w), 1)
            self.assertIsInstance(w[-1].message, SyntaxWarning)
        self.sphinx.exchange_correlation_functional = "pbe"
        self.assertEqual(self.sphinx.exchange_correlation_functional, "PBE")

    def test_write_structure(self):
        cell = (self.sphinx.structure.cell / BOHR_TO_ANGSTROM).tolist()
        pos_2 = (self.sphinx.structure.positions[1] / BOHR_TO_ANGSTROM).tolist()
        file_content = [
            "cell = " + str(cell) + ";\n",
            "species {\n",
            '\telement = "Fe";\n',
            "\tatom {\n",
            '\t\tlabel = "spin_0.5";\n',
            "\t\tcoords = [0.0, 0.0, 0.0];\n",
            "\t}\n",
            "\tatom {\n",
            '\t\tlabel = "spin_0.5";\n',
            "\t\tcoords = " + str(pos_2) + ";\n",
            "\t}\n",
            "}\n",
        ]
        file_name = os.path.join(
            self.file_location,
            "../static/sphinx/job_sphinx_hdf5/job_sphinx/structure.sx",
        )
        with open(file_name) as structure_sx:
            lines = structure_sx.readlines()
        self.assertEqual(file_content, lines)

    def test_collect_aborted(self):
        with self.assertRaises(AssertionError):
            self.sphinx_aborted.collect_output()

    def test_collect_2_5(self):
        output = self.sphinx_2_5._output_parser
        output.collect(directory=self.sphinx_2_5.working_directory)
        self.assertTrue(
            all(
                (
                        output._parse_dict["scf_computation_time"][0]
                        - np.roll(output._parse_dict["scf_computation_time"][0], 1)
                )[1:]
                > 0
            )
        )
        self.assertTrue(
            all(
                np.array(output._parse_dict["scf_energy_free"][0])
                - np.array(output._parse_dict["scf_energy_int"][0])
                < 0
            )
        )
        self.assertTrue(
            all(
                np.array(output._parse_dict["scf_energy_free"][0])
                - np.array(output._parse_dict["scf_energy_zero"][0])
                < 0
            )
        )
        list_values = [
            "scf_energy_int",
            "scf_energy_zero",
            "scf_energy_free",
            "scf_convergence",
            "scf_electronic_entropy",
            "atom_scf_spins",
        ]
        for list_one in list_values:
            for list_two in list_values:
                self.assertEqual(
                    len(output._parse_dict[list_one]), len(output._parse_dict[list_two])
                )

    def test_collect_2_3(self):
        file_location = os.path.join(
            self.file_location, "../static/sphinx/sphinx_test_2_3_hdf5/sphinx_test_2_3/"
        )
        residue_lst = np.loadtxt(file_location + "residue.dat")[:, 1].reshape(1, -1)
        residue_lst = (residue_lst * HARTREE_TO_EV).tolist()
        energy_int_lst = np.loadtxt(file_location + "energy.dat")[:, 2].reshape(1, -1)
        energy_int_lst = (energy_int_lst * HARTREE_TO_EV).tolist()
        with open(file_location + "sphinx.log") as ffile:
            energy_free_lst = [[float(line.split('=')[-1]) * HARTREE_TO_EV for line in ffile if line.startswith('F(')]]
        energy_zero_lst = [(0.5 * (np.array(ff) + np.array(uu))).tolist() for ff, uu in
                           zip(energy_free_lst, energy_int_lst)]
        eig_lst = [np.loadtxt(file_location + "eps.dat")[:, 1:].tolist()]
        self.sphinx_2_3.collect_output()
        self.assertEqual(
            residue_lst, self.sphinx_2_3._output_parser._parse_dict["scf_residue"]
        )
        self.assertEqual(
            energy_int_lst, self.sphinx_2_3._output_parser._parse_dict["scf_energy_int"]
        )
        self.assertEqual(
            energy_zero_lst,
            self.sphinx_2_3._output_parser._parse_dict["scf_energy_zero"],
        )
        self.assertEqual(
            eig_lst,
            self.sphinx_2_3._output_parser._parse_dict["bands_eigen_values"].tolist(),
        )
        self.assertEqual(
            energy_free_lst,
            self.sphinx_2_3._output_parser._parse_dict["scf_energy_free"],
        )
        self.assertEqual(
            21.952 * BOHR_TO_ANGSTROM ** 3, self.sphinx_2_3._output_parser._parse_dict["volume"]
        )

    def test_structure_parsing(self):
        self.sphinx_2_3._output_parser.collect_relaxed_hist(
            file_name="relaxedHist_2.sx", cwd=self.sphinx_2_3.working_directory
        )


if __name__ == "__main__":
    unittest.main()
