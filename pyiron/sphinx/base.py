# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, division

import numpy as np
import os
import posixpath
import re
import stat
from shutil import copyfile
import scipy.constants
import subprocess
from ase import io
from pyiron.atomistics.structure.atoms import ase_to_pyiron
import warnings
import json
from collections import OrderedDict as odict
from collections import defaultdict

from pyiron.dft.job.generic import GenericDFTJob
from pyiron.vasp.potential import VaspPotentialFile
from pyiron.vasp.potential import find_potential_file as find_potential_file_vasp
from pyiron.sphinx.potential import SphinxJTHPotentialFile
from pyiron.sphinx.potential import find_potential_file as find_potential_file_jth
from pyiron.base.settings.generic import Settings
from pyiron.base.generic.parameters import GenericParameters

__author__ = "Osamu Waseda, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"

s = Settings()

BOHR_TO_ANGSTROM = (
    scipy.constants.physical_constants["Bohr radius"][0] / scipy.constants.angstrom
)
HARTREE_TO_EV = scipy.constants.physical_constants["Hartree energy in eV"][0]
HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM = HARTREE_TO_EV / BOHR_TO_ANGSTROM


class SphinxBase(GenericDFTJob):
    """
    Class to setup and run Sphinx simulations which is a derivative of pyiron_atomistics.job.generic.GenericJob.
    The functions in these modules are written in such the function names and attributes are very generic
    (get_structure(), molecular_dynamics(), version) but the functions are written to handle Sphinx specific input and
    output.

    Args:
        project: Project object (defines path where job will be created and stored)
        job_name (str): name of the job (must be unique within this project path)
    """

    def __init__(self, project, job_name):
        super(SphinxBase, self).__init__(project, job_name)
        self.input = Input()
        self._main_str = None
        self._species_str = None
        self._structure_str = None
        self._basis_str = None
        self._hamilston_str = None
        self._guess_str = None
        self._spins_str = None
        self._save_memory = False
        self._output_parser = Output(self)
        self.input_writer = InputWriter()
        if self.check_vasp_potentials():
            self.input["PotType"] = "VASP"  # use VASP potentials if available
        self._kpoints_odict = None
        self._generic_input["restart_for_band_structure"] = False
        self._generic_input["path_name"] = None
        self._generic_input["n_path"] = None

    @property
    def id_pyi_to_spx(self):
        if self.input_writer.id_pyi_to_spx is None:
            self.input_writer.structure = self.structure
        return self.input_writer.id_pyi_to_spx

    @property
    def id_spx_to_pyi(self):
        if self.input_writer.id_spx_to_pyi is None:
            self.input_writer.structure = self.structure
        return self.input_writer.id_spx_to_pyi

    @property
    def plane_wave_cutoff(self):
        return self.input["EnCut"]

    @property
    def fix_spin_constraint(self):
        return self._generic_input["fix_spin_constraint"]

    @fix_spin_constraint.setter
    def fix_spin_constraint(self, boolean):
        if not isinstance(boolean, bool):
            raise ValueError("fix_spin_constraint has to be a boolean")
        self._generic_input["fix_spin_constraint"] = boolean
        self.structure.add_tag(spin_constraint=boolean)

    @plane_wave_cutoff.setter
    def plane_wave_cutoff(self, val):
        """
        Function to setup the energy cut-off for the Sphinx job in eV.

        Args:
            encut (int): energy cut-off in eV
        """
        if val <= 0:
            raise ValueError("Cutoff radius value not valid")
        self.input["EnCut"] = val

    @property
    def exchange_correlation_functional(self):
        return self.input["Xcorr"]

    @exchange_correlation_functional.setter
    def exchange_correlation_functional(self, val):
        """

        Args:
            exchange_correlation_functional:

        Returns:

        """
        if val.upper() in ["PBE", "LDA"]:
            self.input["Xcorr"] = val.upper()
        else:
            warnings.warn(
                "Exchange correlation function not recognized (recommended: PBE or LDA)",
                SyntaxWarning,
            )
            self.input["Xcorr"] = val

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(SphinxBase, self).set_input_to_read_only()
        self.input.read_only = True

    def _input_control_scf_string(
        self, maxSteps=None, keepRhoFixed=False, dEnergy=None, algorithm="blockCCG"
    ):
        """
            scf control string setting for SPHInX
            for all args refer to calc_static or calc_minimize
        """
        if algorithm.upper() == "CCG":
            algorithm = "CCG"
        else:
            if algorithm.upper() != "BLOCKCCG":
                warnings.warn(
                    "Algorithm not recognized -> setting to blockCCG. Alternatively, choose algorithm=CCG",
                    SyntaxWarning,
                )
            algorithm = "blockCCG"
        control_str = odict()
        if keepRhoFixed:
            control_str["keepRhoFixed"] = None
        else:
            control_str["rhoMixing"] = str(self.input["rhoMixing"])
            control_str["spinMixing"] = str(self.input["spinMixing"])
            if self.input["nPulaySteps"] is not None:
                control_str["nPulaySteps"] = str(self.input["nPulaySteps"])
        if dEnergy is None:
            control_str["dEnergy"] = "Ediff/" + str(HARTREE_TO_EV)
        else:
            control_str["dEnergy"] = str(dEnergy)
        if maxSteps is None:
            control_str["maxSteps"] = str(self.input["Estep"])
        else:
            control_str["maxSteps"] = str(maxSteps)
        if self.input["preconditioner"] is not None:
            control_str["preconditioner"] = odict(
                [("type", self.input["preconditioner"])]
            )
        control_str[algorithm] = odict()
        if self.input["maxStepsCCG"] is not None:
            control_str[algorithm]["maxStepsCCG"] = self.input["maxStepsCCG"]
        if self.input["blockSize"] is not None and algorithm == "blockCCG":
            control_str[algorithm]["blockSize"] = self.input["blockSize"]
        if self.input["nSloppy"] is not None and algorithm == "blockCCG":
            control_str[algorithm]["nSloppy"] = self.input["nSloppy"]
        if self.input["WriteWaves"] is False:
            control_str["noWavesStorage"] = None
        return control_str

    @property
    def _control_str(self):
        control_str = odict()
        control_str.setdefault("scfDiag", [])
        if len(self.restart_file_list) != 0 and not self._generic_input["restart_for_band_structure"]:
            control_str["scfDiag"].append(
                self._input_control_scf_string(
                    maxSteps=10, keepRhoFixed=True, dEnergy=1.0e-4
                )
            )
        if self.input["Istep"] is not None:
            control_str["linQN"] = odict()
            control_str["linQN"]["maxSteps"] = str(self.input["Istep"])
            if self.input["dE"] is None and self.input["dF"] is None:
                self.input["dE"] = 1e-3
            if self.input["dE"] is not None:
                control_str["linQN"]["dEnergy"] = str(self.input["dE"] / HARTREE_TO_EV)
            if self.input["dF"] is not None:
                control_str["linQN"]["dF"] = str(
                    self.input["dF"] / HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM
                )
            control_str["linQN"]["bornOppenheimer"] = odict(
                [("scfDiag", self._input_control_scf_string())]
            )
        else:
            if self._generic_input["restart_for_band_structure"]:
                control_str["scfDiag"].append(self._input_control_scf_string(keepRhoFixed=True))
            else:
                control_str["scfDiag"].append(self._input_control_scf_string())
            if self.executable.version is not None:
                vers_num = [
                    int(vv) for vv in self.executable.version.split("_")[0].split(".")
                ]
                if vers_num[0] > 2 or (vers_num[0] == 2 and vers_num[1] > 5):
                    control_str["evalForces"] = odict([("file", '"relaxHist.sx"')])
            else:
                warnings.warn("executable version could not be identified")
        return control_str

    def calc_static(
        self,
        electronic_steps=None,
        algorithm=None,
        retain_charge_density=False,
        retain_electrostatic_potential=False,
    ):
        """
        Function to setup the hamiltonian to perform static SCF DFT runs

        Args:
            retain_electrostatic_potential:
            retain_charge_density:
            algorithm (str): CCG or blockCCG (not implemented)
            electronic_steps (int): maximum number of electronic steps, which can be used
                                    to achieve convergence
        """
        if electronic_steps is not None:
            self.input["Estep"] = electronic_steps
        for arg in ["Istep", "dF", "dE"]:
            if self.input[arg] is not None:
                del self.input[arg]
        super(SphinxBase, self).calc_static(
            electronic_steps=electronic_steps,
            algorithm=algorithm,
            retain_charge_density=retain_charge_density,
            retain_electrostatic_potential=retain_electrostatic_potential,
        )

    def calc_minimize(
        self,
        electronic_steps=None,
        ionic_steps=None,
        max_iter=None,
        pressure=None,
        algorithm=None,
        retain_charge_density=False,
        retain_electrostatic_potential=False,
        ionic_energy=None,
        ionic_forces=None,
        volume_only=False,
    ):
        """
        Function to setup the hamiltonian to perform ionic relaxations using DFT. The convergence goal can be set using
        either the iconic_energy as an limit for fluctuations in energy or the iconic_forces.
        Args:
            retain_electrostatic_potential:
            retain_charge_density:
            algorithm:
            pressure:
            max_iter:
            electronic_steps (int): maximum number of electronic steps per electronic convergence
            ionic_steps (int): maximum number of ionic steps
            ionic_energy (float): convergence goal in terms of energy (optional)
            ionic_forces (float): convergence goal in terms of forces (optional)
        """
        if pressure is not None:
            raise NotImplementedError(
                "pressure minimization is not implemented in SPHInX"
            )
        if electronic_steps is not None:
            self.input["Estep"] = electronic_steps
        if ionic_steps is not None:
            self.input["Istep"] = ionic_steps
        elif self.input["Istep"] is None:
            self.input["Istep"] = 100
        if ionic_forces is not None:
            if ionic_forces < 0:
                raise ValueError("ionic_forces must be a positive integer")
            self.input["dF"] = float(ionic_forces)
        if ionic_energy is not None:
            if ionic_energy < 0:
                raise ValueError("ionic_forces must be a positive integer")
            self.input["dE"] = float(ionic_energy)
        super(SphinxBase, self).calc_minimize(
            electronic_steps=electronic_steps,
            ionic_steps=ionic_steps,
            max_iter=max_iter,
            pressure=pressure,
            algorithm=algorithm,
            retain_charge_density=retain_charge_density,
            retain_electrostatic_potential=retain_electrostatic_potential,
            ionic_energy=ionic_energy,
            ionic_forces=ionic_forces,
            volume_only=volume_only,
        )

    def calc_md(
        self,
        temperature=None,
        n_ionic_steps=1000,
        n_print=1,
        time_step=1.0,
        retain_charge_density=False,
        retain_electrostatic_potential=False,
        **kwargs
    ):
        raise NotImplementedError("calc_md() not implemented in SPHInX.")

    def restart_from_charge_density(
            self,
            job_name=None,
            job_type=None,
            band_structure_calc=False
    ):
        """
        Restart a new job created from an existing Vasp calculation by reading the charge density.

        Args:
            job_name (str): Job name
            job_type (str): Job type. If not specified a Vasp job type is assumed
            band_structure_calc (bool): has to be True for band structure calculations.

        Returns:
            pyiron.sphinx.sphinx.sphinx: new job instance
        """
        ham_new = self.restart(
            job_name=job_name,
            job_type=job_type,
            from_wave_functions=False,
            from_charge_density=True
        )
        if band_structure_calc:
            ham_new._generic_input["restart_for_band_structure"] = True
        return ham_new

    def restart_from_wave_functions(
            self,
            job_name=None,
            job_type=None,
    ):
        """
        Restart a new job created from an existing Vasp calculation by reading the wave functions.

        Args:
            job_name (str): Job name
            job_type (str): Job type. If not specified a Vasp job type is assumed

        Returns:
            pyiron.sphinx.sphinx.sphinx: new job instance
        """
        return self.restart(
            job_name=job_name,
            job_type=job_type,
            from_wave_functions=True,
            from_charge_density=False
        )

    def restart(
        self,
        job_name=None,
        job_type=None,
        from_charge_density=True,
        from_wave_functions=True,
    ):
        if self.status!='finished' and not self.is_compressed():
            # self.decompress()
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                try:
                    self.collect_output()
                except AssertionError:
                    from_charge_density=False
                    from_wave_functions=False
                if len(w) > 0:
                    self.status.not_converged = True
        new_job = super(SphinxBase, self).restart(
            job_name=job_name, job_type=job_type
        )
        if from_charge_density and os.path.isfile(
            posixpath.join(self.working_directory, "rho.sxb")
        ):
            new_job.restart_file_list.append(
                posixpath.join(self.working_directory, "rho.sxb")
            )
        elif from_charge_density:
            self._logger.warning(
                msg="A charge density from job: {} is not generated and therefore it can't be read.".format(
                    self.job_name
                )
            )
        if from_wave_functions and os.path.isfile(
            posixpath.join(self.working_directory, "waves.sxb")
        ):
            new_job.restart_file_list.append(
                posixpath.join(self.working_directory, "waves.sxb")
            )
        elif from_wave_functions:
            self._logger.warning(
                msg="A WAVECAR from job: {} is not generated and therefore it can't be read.".format(
                    self.job_name
                )
            )
        return new_job

    def to_hdf(self, hdf=None, group_name=None):
        """
        Stores the instance attributes into the hdf5 file

        Args:
            hdf (str): Path to the hdf5 file
            group_name (str): Name of the group which contains the object
        """
        super(SphinxBase, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()
        self.input.to_hdf(self._hdf5)
        self._output_parser.to_hdf(self._hdf5)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Recreates instance from the hdf5 file

        Args:
            hdf (str): Path to the hdf5 file
            group_name (str): Name of the group which contains the object
        """
        super(SphinxBase, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()
        self.input.from_hdf(self._hdf5)
        if self.status.finished:
            self._output_parser.from_hdf(self._hdf5)

    def from_directory(self, directory):
        try:
            if not self.status.finished:
                subprocess.call(
                    "module load sphinx && sx2aims", cwd=directory, shell=True
                )
                # self._output_parser.to_hdf(self._hdf5)
                if directory[-1] == "/":
                    directory = directory[:-1]
                if os.path.isfile(directory + "/geometry.in"):
                    self.structure = ase_to_pyiron(
                        io.read(filename=directory + "/geometry.in")
                    )
                else:
                    print(
                        "WARNING: input structure not found: "
                        + directory
                        + "/geometry.in"
                    )
                subprocess.call(
                    "rm " + directory + "/geometry.in", cwd=directory, shell=True
                )
                self._output_parser.collect(directory=directory)
                self.to_hdf(self._hdf5)
            else:
                self._output_parser.from_hdf(self._hdf5)
            self.status.finished = True
        except Exception as err:
            print(err)
            self.status.aborted = True

    def set_check_overlap(self, check_overlap=True):
        """
        Args:
            check_overlap (bool): Whether to check overlap

        Comments:
            Certain PAW-pseudo-potentials have an intrinsic pathology: their PAW overlap
            operator is not generally positive definite (i.e., the PAW-corrected
            norm of a wavefunction could become negative). SPHInX usually refuses to
            use such problematic potentials. This behavior can be overridden by setting
            check_overlap to False.
        """
        if not isinstance(check_overlap, bool):
            raise ValueError("check_overlap has to be a boolean")
        if self.executable.version != "2.5.1" and not check_overlap:
            vers_num = [
                int(vv) for vv in self.executable.version.split("_")[0].split(".")
            ]
            if (
                vers_num[0] < 2
                or vers_num[1] < 5
                or (vers_num[0] <= 2 and sum(vers_num[1:]) <= 5)
            ):
                warnings.warn(
                    "SPHInX executable version has to be 2.5.1 or above in order for the overlap to be considered. "
                    + "Change it via job.executable.version"
                )
        self.input["CheckOverlap"] = check_overlap

    def set_mixing_parameters(
        self,
        method=None,
        n_pulay_steps=None,
        density_mixing_parameter=None,
        spin_mixing_parameter=None,
    ):
        """
        args:
            method ('PULAY' or 'LINEAR'): mixing method (default: PULAY)
            n_pulay_steps (int): number of previous densities to use for the Pulay mixing (default: 7)
            density_mixing_parameter (float): mixing proportion m defined by rho^n = (m-1)*rho^(n-1)+m*preconditioner*rho_(opt) (default: 1)
            spin_mixing_parameter (float): linear mixing parameter for spin densities (default: 1)

        comments:
            A low value of density mixing parameter may lead to a more stable convergence,
            but will slow down the calculation if set too low.

            Further information can be found on the website: https://sxrepo.mpie.de
        """
        method_list = ["PULAY", "LINEAR"]
        assert (
            method is None or method.upper() in method_list
        ), "Mixing method has to be PULAY or LINEAR"
        assert n_pulay_steps is None or isinstance(
            n_pulay_steps, int
        ), "n_pulay_steps has to be an integer"
        if density_mixing_parameter is not None and (
            density_mixing_parameter > 1.0 or density_mixing_parameter < 0
        ):
            raise ValueError(
                "density_mixing_parameter has to be between 0 and 1 (default value is 1)"
            )
        if spin_mixing_parameter is not None and (
            spin_mixing_parameter > 1.0 or spin_mixing_parameter < 0
        ):
            raise ValueError(
                "spin_mixing_parameter has to be between 0 and 1 (default value is 1)"
            )

        if method is not None:
            self.input["mixingMethod"] = method.upper()
        if n_pulay_steps is not None:
            self.input["nPulaySteps"] = n_pulay_steps
        if density_mixing_parameter is not None:
            self.input["rhoMixing"] = density_mixing_parameter
        if spin_mixing_parameter is not None:
            self.input["spinMixing"] = spin_mixing_parameter

    def set_occupancy_smearing(self, smearing=None, width=None):
        """
        Set how the finite temperature smearing is applied in determining partial occupancies

        Args:
            smearing (str): Type of smearing (only fermi si implemented anything else will be ignored)
            width (float): Smearing width (eV) (default: 0.2)
        """
        if smearing is not None and not isinstance(smearing, str):
            raise ValueError(
                "Smearing must be a string (only fermi is supported in SPHInX)"
            )
        if width is not None and width < 0:
            raise ValueError("Smearing value must be a float >= 0")
        if width is not None:
            self.input["Sigma"] = width

    def set_convergence_precision(
        self, ionic_energy=None, electronic_energy=None, ionic_forces=None
    ):
        """
        Sets the electronic and ionic convergence precision. For ionic convergence either the energy or the force
        precision is required

        Args:
            ionic_energy (float): Ionic energy convergence precision (eV)
            electronic_energy (float): Electronic energy convergence precision (eV)
            ionic_forces (float): Ionic force convergence precision (eV/A)
        """
        assert (
            ionic_energy is None or ionic_energy > 0
        ), "ionic_energy must be a positive float"
        assert (
            ionic_forces is None or ionic_forces > 0
        ), "ionic_forces must be a positive float"
        assert (
            electronic_energy is None or electronic_energy > 0
        ), "electronic_energy must be a positive float"
        if ionic_energy is not None or ionic_forces is not None:
            print("Setting calc_minimize")
            self.calc_minimize(ionic_energy=ionic_energy, ionic_forces=ionic_forces)
        if electronic_energy is not None:
            self.input["Ediff"] = electronic_energy

    def set_empty_states(self, n_empty_states=None):
        """
        Function to set the number of empty states.

        Args:
            n_empty_states (int/None): Number of empty states. If None, sets it to 'auto'.

        Comments:
            If this number is too low, the algorithm will not be able to able to swap wave functions
            near the chemical potential. If the number is too high, computation time will be wasted
            for the higher energy states and potentially lead to a memory problem.

            In contrast to VASP, this function sets only the number of empty states and not the number of
            total states.

            The default value is 0.5*NIONS+3 for non-magnetic systems and 1.5*NIONS+3 for magnetic systems
        """
        if n_empty_states is None:
            self.input["EmptyStates"] = "auto"
        else:
            if n_empty_states < 0:
                raise ValueError("Number of empty states must be greater than 0")
            self.input["EmptyStates"] = n_empty_states

    def _set_kpoints(
        self,
        mesh=None,
        scheme="MP",
        center_shift=None,
        symmetry_reduction=True,
        manual_kpoints=None,
        weights=None,
        reciprocal=True,
        n_path=None,
        path_name=None,
    ):
        """
        Function to setup the k-points for the Sphinx job

        Args:
            reciprocal (bool): Tells if the supplied values are in reciprocal (direct)
                               or cartesian coordinates (in reciprocal space) (not implemented)
            weights (list): Manually supplied weights to each k-point in case of the manual mode (not implemented)
            manual_kpoints (list): Manual list of k-points (not implemented)
            symmetry_reduction (bool): Tells if the symmetry reduction is to be applied to the k-points
            scheme (str): Type of k-point generation scheme ('MP' or 'Line')
            mesh (list): Size of the mesh (in the MP scheme)
            center_shift (list): Shifts the center of the mesh from the gamma point by the given vector
            n_trace (int): Number of points per trace part for line mode
            path_name (str): Name of high symmetry path used for band structure calculations.
        """
        if not isinstance(symmetry_reduction, bool):
            raise ValueError("symmetry_reduction has to be a boolean")
        if manual_kpoints is not None:
            raise ValueError("manual_kpoints is not implemented in SPHInX yet")
        if weights is not None:
            raise ValueError("manual weights are not implmented in SPHInX yet")
        if scheme == "MP":
            self._kpoints_odict = None
            if mesh is not None:
                self.input["KpointFolding"] = str(list(mesh))
            if center_shift is not None:
                self.input["KpointCoords"] = str(list(center_shift))
        elif scheme == "Line":
            if n_path is None and self._generic_input["n_path"] is None:
                raise ValueError("'n_path' has to be defined")
            if n_path is None:
                n_path = self._generic_input["n_path"]
            else:
                self._generic_input["n_path"] = n_path

            if self.structure.get_high_symmetry_points() is None:
                raise ValueError("no 'high_symmetry_points' defined for 'structure'.")

            if path_name is None and self._generic_input["path_name"] is None:
                raise ValueError("'path_name' has to be defined")
            if path_name is None:
                path_name = self._generic_input["path_name"]
            else:
                self._generic_input["path_name"] = path_name

            try:
                path = self.structure.get_high_symmetry_path()[path_name]
            except KeyError:
                raise AssertionError("'{}' is not a valid path!".format(path_name))

            kpoints = odict([("relative", None)])

            kpoints["from"] = odict(
                [
                    ("coords", str(self.structure.get_high_symmetry_points()[path[0][0]])),
                    ("label", '"' + path[0][0].replace("'", "p") + '"'),
                ]
            )
            kpoints["to___0"] = odict(
                [
                    ("coords", str(self.structure.get_high_symmetry_points()[path[0][1]])),
                    ("nPoints", n_path),
                    ("label", '"' + path[0][1].replace("'", "p") + '"'),
                ]
            )

            for i, path in enumerate(zip(path[:-1], np.roll(path, -1, 0)[:-1])):
                if not path[0][1] == path[1][0]:
                    name = "to___{}___1".format(i)
                    kpoints[name] = odict(
                        [
                            ("coords", str(self.structure.get_high_symmetry_points()[path[1][0]])),
                            ("nPoints", 0),
                            ("label", '"' + path[1][0].replace("'", "p") + '"'),
                        ]
                    )

                name = "to___{}".format(i + 1)
                kpoints[name] = odict(
                    [
                        ("coords", str(self.structure.get_high_symmetry_points()[path[1][1]])),
                        ("nPoints", n_path),
                        ("label", '"' + path[1][1].replace("'", "p") + '"'),
                    ]
                )

            self._kpoints_odict = odict([("kPoints", kpoints)])
        else:
            raise ValueError("only Monkhorst-Pack mesh and Line mode are implemented in SPHInX")

    def write_input(self):
        """
        The write_input function is called when the job is executed to generate all the required input files for the
        calculation of the Sphinx job.
        """
        self._coarse_run = self.input["CoarseRun"]
        if self.input["EmptyStates"] == "auto":
            self.input["EmptyStates"] = int(len(self.structure) + 3)
            if np.any(self.structure.get_initial_magnetic_moments() != None):
                self.input["EmptyStates"] = int(1.5 * len(self.structure) + 3)
        self.input_writer.structure = self.structure
        write_waves = self.input["WriteWaves"]
        save_memory = self.input["SaveMemory"]
        check_overlap = self.input["CheckOverlap"]
        enable_kjxc = self.input["KJxc"]
        if self._main_str is None:
            self.input_writer.write_potentials(
                file_name="potentials.sx",
                cwd=self.working_directory,
                species_str=self._species_str,
                check_overlap=check_overlap,
                xc=self.input["Xcorr"],
                potformat=self.input["PotType"]
            )
            self.input_writer.write_guess(
                file_name="guess.sx",
                cwd=self.working_directory,
                guess_str=self._guess_str,
                restart_file_str=self.restart_file_list,
                write_waves=write_waves,
            )
            self.input_writer.write_structure(
                file_name="structure.sx",
                cwd=self.working_directory,
                structure_str=self._structure_str,
                symmetry_enabled=self._generic_input["fix_symmetry"],
            )
            self.input_writer.write_control(
                file_name="control.sx",
                cwd=self.working_directory,
                control_str=self._control_str,
            )
            self.input_writer.write_basis(
                file_name="basis.sx",
                cwd=self.working_directory,
                basis_str=self._basis_str,
                save_memory=save_memory,
                kpoints_odict=self._kpoints_odict
            )
            self.input_writer.write_hamilton(
                file_name="hamilton.sx",
                cwd=self.working_directory,
                hamilston_str=self._hamilston_str,
            )
            if self._generic_input["fix_spin_constraint"]:
                self.input_writer.write_spins_constraints(
                    file_name="spins.in",
                    cwd=self.working_directory,
                    spins_str=self._spins_str,
                )
        self.input.write_file(file_name="userparameters.sx", cwd=self.working_directory)
        self.input_writer.write_main(
            file_name="input.sx",
            cwd=self.working_directory,
            main_str=self._main_str,
            spin_constraint_enabled=self._generic_input["fix_spin_constraint"],
            enable_kjxc=enable_kjxc,
            job_name=self.job_name,
        )

    def collect_output(self, force_update=False):
        """
        Collects the outputs and stores them to the hdf file
        """
        self._output_parser.collect(directory=self.working_directory)
        self._output_parser.to_hdf(self._hdf5, force_update=force_update)

    def convergence_check(self):
        if (
            self._generic_input["calc_mode"] == "minimize"
            and self._output_parser._parse_dict["scf_convergence"][-1]
        ):
            return True
        elif self._generic_input["calc_mode"] == "static" and np.all(
            self._output_parser._parse_dict["scf_convergence"]
        ):
            return True
        else:
            return False

    def collect_logfiles(self):
        """
        Collect errors and warnings.
        """
        self.collect_errors()
        self.collect_warnings()

    def collect_warnings(self):
        """
        Collects warnings from the Sphinx run
        """
        # TODO: implement for Sphinx
        self._logger.info("collect_warnings() is not yet implemented for Sphinx")

    def collect_errors(self):
        """
        Collects errors from the Sphinx run
        """
        # TODO: implement for Sphinx
        self._logger.info("collect_errors() is not yet implemented for Sphinx")

    def get_n_ir_reciprocal_points(
        self, is_time_reversal=True, symprec=1e-5, ignore_magmoms=False
    ):
        from phonopy.structure import spglib

        lattice = self.structure.cell
        positions = self.structure.get_scaled_positions()
        numbers = self.structure.get_atomic_numbers()
        magmoms = self.structure.get_initial_magnetic_moments()
        if np.all(magmoms == None) or ignore_magmoms:
            magmoms = np.zeros(len(magmoms))
        mag_num = np.array(list(zip(magmoms, numbers)))
        satz = np.unique(mag_num, axis=0)
        numbers = []
        for nn in np.all(satz == mag_num[:, np.newaxis], axis=-1):
            numbers.append(np.arange(len(satz))[nn][0])
        mapping, _ = spglib.get_ir_reciprocal_mesh(
            mesh=[int(k) for k in self.input["KpointFolding"]],
            cell=(lattice, positions, numbers),
            is_shift=np.dot(self.structure.cell, np.array(self.input["KpointCoords"])),
            is_time_reversal=is_time_reversal,
            symprec=symprec,
        )
        return len(np.unique(mapping))

    def check_setup(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            if (
                not (
                    isinstance(self.input["EnCut"], int)
                    or isinstance(self.input["EnCut"], float)
                )
                or self.input["EnCut"] == 340
            ):
                warnings.warn(
                    "Energy cut-off value wrong or not modified from default 340 eV; change it via job.set_encut()"
                )
            if not (
                isinstance(self.input["KpointCoords"], list)
                or len(self.input["KpointCoords"]) != 3
            ):
                warnings.warn("K point coordinates seem to be inappropriate")
            if (
                not (
                    isinstance(self.input["Sigma"], int)
                    or isinstance(self.input["Sigma"], float)
                )
                or self.input["Sigma"] == 0.2
            ):
                warnings.warn(
                    "Fermi smearing value wrong or not modified from default 0.2 eV; change it via job.set_occupancy_smearing()"
                )
            if not (
                isinstance(self.input["KpointFolding"], list)
                or len(self.input["KpointFolding"]) != 3
            ) or self.input["KpointFolding"] == [4, 4, 4]:
                warnings.warn(
                    "K point folding wrong or not modified from default [4,4,4]; change it via job.set_kpoints()"
                )
            if self.get_n_ir_reciprocal_points() < self.server.cores:
                warnings.warn(
                    "Number of cores exceed number of irreducible reciprocal points: "
                    + str(self.get_n_ir_reciprocal_points())
                )
            if self.input["EmptyStates"] == "auto":
                if any(self.structure.get_initial_magnetic_moments() != None):
                    warnings.warn(
                        "Number of empty states was not specified. Default: NIONS*1.5 for magnetic systems. "
                    )
                else:
                    warnings.warn(
                        "Number of empty states was not specified. Default: NIONS for non-magnetic systems"
                    )
            if len(w) > 0:
                print("WARNING:")
                for ww in w:
                    print(ww.message)
                return False
            else:
                return True

    def validate_ready_to_run(self):
        """
            Checks whether parameters are set appropriately. It does not mean the simulation won't run even if it returns False
        """
        if self._control_str is None:
            self.calc_static()
        if self.structure is None:
            raise AssertionError(
                "Structure not set; set it via job.structure = Project().create_structure()"
            )
        if self._control_str is None:
            raise AssertionError(
                "Control string not set; set it e.g. via job.calc_static()"
            )
        if self.input["THREADS"] > self.server.cores:
            raise AssertionError(
                "Number of cores cannot be smaller than the number of OpenMP threads"
            )
        if self.input["EmptyStates"] != "auto" and self.input["EmptyStates"] < 0:
            raise AssertionError("Number of empty states not valid")

    def compress(self, files_to_compress=None):
        """
        Compress the output files of a job object.

        Args:
            files_to_compress (list): A list of files to compress (optional)
        """
        # delete empty files
        if files_to_compress is None:
            files_to_compress = [
                f for f in list(self.list_files()) if (f not in ["rho.sxb", "waves.sxb"]
                                                       and not stat.S_ISFIFO(os.stat(os.path.join(self.working_directory, f)).st_mode))
            ]
        for f in list(self.list_files()):
            filename = os.path.join(self.working_directory, f)
            if (
                f not in files_to_compress
                and os.path.exists(filename)
                and os.stat(filename).st_size == 0
            ):
                os.remove(filename)
        super(SphinxBase, self).compress(files_to_compress=files_to_compress)

    @staticmethod
    def check_vasp_potentials():
        return any([os.path.exists(os.path.join(p, 'vasp')) for p in s.resource_paths])


class InputWriter(object):
    """
    The Sphinx Input writer is called to write the Sphinx specific input files.
    """

    def __init__(self):
        self.structure = None
        self._spin_enabled = False
        self._id_pyi_to_spx = []
        self._id_spx_to_pyi = []
        self.file_dict = {}

    def _odict_to_spx_input(self, element, level=0):
        """
        Convert collections.OrderedDict containing SPHInX input
        hierarchy to string.

        If the item contains:
            - no value -> considered as flag -> output format: flag;
            - value is odict
                -> considered as a group
                -> output format: group { ...recursive... }
            - else
                -> considered as parameter and value
                -> output format: parameter = value;
        """
        line = ""
        for k, v in element.items():
            k = k.split("___")[0]
            if type(v) != list:
                v = [v]
            for vv in v:
                if vv is None:
                    line += level * "\t" + str(k) + ";\n"
                elif type(vv) == odict:
                    if len(vv) == 0:
                        line += level * "\t" + k + " {}\n"
                    else:
                        line += (
                            level * "\t"
                            + k
                            + " {\n"
                            + self._odict_to_spx_input(vv, level + 1)
                            + level * "\t"
                            + "}\n"
                        )
                else:
                    line += level * "\t" + k + " = " + str(vv) + ";\n"
        return line

    def write_main(
        self,
        file_name="input.sx",
        cwd=None,
        main_str=None,
        spin_constraint_enabled=False,
        enable_kjxc=False,
        job_name=None,
    ):
        """
        Write the main Sphinx script named input.sx.

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            main_str (str): the input to write;
                            if no input is given,
                            the default input will be written. (optinal)
        """
        if main_str is None:
            self.file_dict['input'] = self.get_main(enable_kjxc=enable_kjxc,
                                                    spin_constraint_enabled=spin_constraint_enabled,
                                                    job_name=job_name)
        else:
            self.file_dict['input'] = main_str
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(self._odict_to_spx_input(self.file_dict['input']))

    @staticmethod
    def get_main(enable_kjxc=False, spin_constraint_enabled=False, job_name=None):
        line = odict(
            [
                ("//" + str(job_name), None),
                ("//SPHinX input file generated by pyiron", None),
                ("format paw", None),
                ("include <parameters.sx>", None),
                ("include <userparameters.sx>", None),
            ]
        )
        if enable_kjxc:
            line["pawPot"] = odict(
                [("include <potentials.sx>", None), ("kjxc", None)]
            )
        else:
            line["pawPot"] = odict([("include <potentials.sx>", None)])
        line["structure"] = odict([("include <structure.sx>", None)])
        line["basis"] = odict([("include <basis.sx>", None)])
        line["PAWHamiltonian"] = odict([("include <hamilton.sx>", None)])
        line["initialGuess"] = odict([("include <guess.sx>", None)])
        if spin_constraint_enabled:
            line["spinConstraint"] = odict([("file", '"spins.in"')])
        line["main"] = odict([("include <control.sx>", None)])
        return line

    @property
    def id_spx_to_pyi(self):
        if self.structure is None:
            return None
        if len(self._id_spx_to_pyi) == 0:
            self._initialize_order()
        return self._id_spx_to_pyi

    @property
    def id_pyi_to_spx(self):
        if self.structure is None:
            return None
        if len(self._id_pyi_to_spx) == 0:
            self._initialize_order()
        return self._id_pyi_to_spx

    def _initialize_order(self):
        for elm_species in self.structure.get_species_objects():
            self._id_pyi_to_spx.append(
                np.arange(len(self.structure))[
                    self.structure.get_chemical_symbols() == elm_species.Abbreviation
                ]
            )
        self._id_pyi_to_spx = np.array(
            [ooo for oo in self._id_pyi_to_spx for ooo in oo]
        )
        self._id_spx_to_pyi = np.array([0] * len(self._id_pyi_to_spx))
        for i, p in enumerate(self._id_pyi_to_spx):
            self._id_spx_to_pyi[p] = i

    def write_potentials(
        self,
        file_name="potentials.sx",
        cwd=None,
        species_str=None,
        check_overlap=True,
        xc=None,
        potformat='JTH',
    ):
        """
        Write the Sphinx potential configuration named potentials.sx.

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optional)
            species_str (str): the input to write, if no input is given the default input will be written. (optional)
        """
        if potformat == 'JTH':
            potentials = SphinxJTHPotentialFile(xc=xc)
            find_potential_file = find_potential_file_jth
            pot_path_dict = {"PBE": "jth-gga-pbe"}
        elif potformat == 'VASP':
            potentials = VaspPotentialFile(xc=xc)
            find_potential_file = find_potential_file_vasp
            pot_path_dict = {"PBE": "paw-gga-pbe", "LDA": "paw-lda"}
        else:
            raise ValueError('Only JTH and VASP potentials are supported!')
        if species_str is None:
            species_str = odict()
            for species_obj in self.structure.get_species_objects():
                if species_obj.Parent is not None:
                    elem = species_obj.Parent
                else:
                    elem = species_obj.Abbreviation

                if "pseudo_potcar_file" in species_obj.tags.keys():
                    new_element = species_obj.tags["pseudo_potcar_file"]
                    potentials.add_new_element(
                        parent_element=elem, new_element=new_element
                    )
                    potential_path = find_potential_file(
                        path=potentials.find_default(new_element)["Filename"].values[0][
                            0
                        ],
                        pot_path_dict=self.pot_path_dict,
                    )
                    assert os.path.isfile(
                        potential_path
                    ), "such a file does not exist in the pp directory"
                else:
                    potential_path = find_potential_file(
                        path=potentials.find_default(elem)["Filename"].values[0][0],
                        pot_path_dict=self.pot_path_dict,
                    )
                if potformat == "JTH":
                    copyfile(potential_path, posixpath.join(cwd, elem + "_GGA.atomicdata"))
                else: 
                    copyfile(potential_path, posixpath.join(cwd, elem + "_POTCAR"))
                check_overlap_str = ""
                species_str.setdefault("species", [])
                if potformat == "JTH":
                    species_str["species"].append(
                        odict(
                            [
                                ("name", '"' + elem + '"'),
                                ("potType", '"AtomPAW"'),
                                ("element", '"' + elem + '"'),
                                ("potential", '"' + elem + "_GGA.atomicdata" + '"'),
                            ]
                        )
                    )
                elif potformat == "VASP":
                    species_str["species"].append(
                        odict(
                            [
                                ("name", '"' + elem + '"'),
                                ("potType", '"VASP"'),
                                ("element", '"' + elem + '"'),
                                ("potential", '"' + elem + "_POTCAR" + '"'),
                            ]
                        )
                    )
                else:
                    raise ValueError()
                if not check_overlap:
                    species_str["species"][-1]["checkOverlap"] = "false"
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(self._odict_to_spx_input(species_str))

    def write_control(self, file_name="control.sx", cwd=None, control_str=None):
        """
        Write the Sphinx control script named control.sx

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            control_str (OrderedDict): the input to write, if no input is given the default input will be written.
        """

        if type(control_str) == odict:
            control_str = self._odict_to_spx_input(control_str)
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(control_str)

    def write_structure(
        self,
        file_name="structure.sx",
        cwd=None,
        structure_str=None,
        symmetry_enabled=True,
        keep_angstrom=False,
    ):
        """
        Write the Sphinx structure file named structure.sx

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            structure_str (str): the input to write, if no input is given the default input will be written. (optinal)
        """
        if structure_str is None:
            if keep_angstrom:
                structure_str = odict(
                    [("cell", str(np.array(self.structure.cell).tolist()))]
                )
            else:
                structure_str = odict(
                    [
                        (
                            "cell",
                            str(
                                np.array(
                                    self.structure.cell * 1 / BOHR_TO_ANGSTROM
                                ).tolist()
                            ),
                        )
                    ]
                )
            if "selective_dynamics" in self.structure._tag_list.keys():
                selective_dynamics_list = self.structure.selective_dynamics.list()
            else:
                selective_dynamics_list = [3 * [False]] * len(self.structure.positions)
            for i_elem, elm_species in enumerate(self.structure.get_species_objects()):
                if elm_species.Parent:
                    element = elm_species.Parent
                else:
                    element = elm_species.Abbreviation
                structure_str.setdefault("species", [])
                structure_str["species"].append(
                    odict([("element", '"' + str(element) + '"')])
                )
                elm_list = np.array(
                    self.structure.get_chemical_symbols() == elm_species.Abbreviation
                )
                for elm_pos, elm_magmon, selective in zip(
                    self.structure.positions[elm_list],
                    np.array(self.structure.get_initial_magnetic_moments())[elm_list],
                    np.array(selective_dynamics_list)[elm_list],
                ):
                    structure_str["species"][-1].setdefault("atom", [])
                    structure_str["species"][-1]["atom"].append(odict())
                    if self._spin_enabled:
                        structure_str["species"][-1]["atom"][-1]["label"] = (
                            '"spin_' + str(elm_magmon) + '"'
                        )
                    if keep_angstrom:
                        structure_str["species"][-1]["atom"][-1]["coords"] = str(
                            np.array(elm_pos).tolist()
                        )
                    else:
                        structure_str["species"][-1]["atom"][-1]["coords"] = str(
                            np.array(elm_pos * 1 / BOHR_TO_ANGSTROM).tolist()
                        )
                    if all(selective):
                        structure_str["species"][-1]["atom"][-1]["movable"] = None
                    elif any(selective):
                        for ss, xx in zip(selective, ["X", "Y", "Z"]):
                            if ss:
                                structure_str["species"][-1]["atom"][-1][
                                    "movable" + xx
                                ] = None
            if not symmetry_enabled:
                structure_str["symmetry"] = odict(
                    [("operator", odict([("S", "[[1,0,0],[0,1,0],[0,0,1]]")]))]
                )
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(self._odict_to_spx_input(structure_str))

    def write_basis(
        self, file_name="basis.sx", cwd=None, basis_str=None, save_memory=False, kpoints_odict=None
    ):
        """
        Write the Sphinx bases set configuration file named basis.sx

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            basis_str (str): the input to write, if no input is given the default input will be written. (optinal)
            kpoints_odict (collection.OrderedDict):
        """

        if kpoints_odict is None:
            kpoints = odict(
                [
                    (
                        "kPoint",
                        odict(
                            [
                                ("coords", "KpointCoords"),
                                ("weight", 1),
                                ("relative", None),
                            ]
                        ),
                    ),
                    ("folding", "KpointFolding")
                ]
            )
        else:
            kpoints = kpoints_odict

        if basis_str is None:
            basis_str = odict([("eCut", "EnCut/13.606")])
            basis_str.update(kpoints)
            if save_memory:
                basis_str["saveMemory"] = None
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(self._odict_to_spx_input(basis_str))

    def write_hamilton(self, file_name="hamilton.sx", cwd=None, hamilston_str=None):
        """
        Write the Sphinx hamilton configuration file named hamilton.sx

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            hamilston_str (list): the input to write, if no input is given the default input will be written. (optinal)
        """
        if hamilston_str is None:
            hamilston_str = odict(
                [("nEmptyStates", "EmptyStates"), ("ekt", "Sigma"), ("xc", "Xcorr")]
            )
        if self._spin_enabled:
            hamilston_str["spinPolarized"] = None
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(self._odict_to_spx_input(hamilston_str))

    def write_guess(
        self,
        file_name="guess.sx",
        cwd=None,
        guess_str=None,
        restart_file_str=[],
        write_waves=False,
    ):
        """
        Write the Sphinx initial guess configuration file named guess.sx

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            guess_str (str): the input to write, if no input is given the default input will be written. (optinal)
        """
        charge_density_file = None
        for ff in restart_file_str:
            if "rho.sxb" in ff.split("/")[-1]:
                charge_density_file = ff
        wave_function_file = None
        for ff in restart_file_str:
            if "waves.sxb" in ff.split("/")[-1]:
                wave_function_file = ff
        if guess_str is None:
            guess_str = odict(
                [("waves", odict([("lcao", odict()), ("pawBasis", None)]))]
            )
            if wave_function_file is not None:
                guess_str["exchange"] = odict(
                    [("file", '"' + wave_function_file + '"')]
                )
            if charge_density_file is None:
                guess_str["rho"] = odict([("atomicOrbitals", None)])
            else:
                guess_str["rho"] = odict([("file", '"' + charge_density_file + '"')])
            if np.any(self.structure.get_initial_magnetic_moments().flatten() != None):
                if any(
                    [
                        True
                        if isinstance(spin, list) or isinstance(spin, np.ndarray)
                        else False
                        for spin in self.structure.get_initial_magnetic_moments()
                    ]
                ):
                    raise ValueError("Sphinx only supports collinear spins.")
                else:
                    for spin in self.structure.get_initial_magnetic_moments()[
                        self.id_pyi_to_spx
                    ]:
                        guess_str["rho"].setdefault("atomicSpin", [])
                        guess_str["rho"]["atomicSpin"].append(
                            odict(
                                [
                                    ("label", '"spin_' + str(spin) + '"'),
                                    ("spin", str(spin)),
                                ]
                            )
                        )
                self._spin_enabled = True
            if not write_waves:
                guess_str["noWavesStorage"] = None
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(self._odict_to_spx_input(guess_str))

    def write_spins_constraints(self, file_name="spins.in", cwd=None, spins_str=None):
        """
        Write a text file containing a list of all spins named spin.in - which is used for the external control scripts.

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            spins_str (str): the input to write, if no input is given the default input will be written. (optinal)
        """
        s.logger.debug("Writing spins.in")
        if spins_str is None:
            s.logger.debug("Getting magnetic moments via get_initial_magnetic_moments")
            if any(self.structure.get_initial_magnetic_moments().flatten() != None):
                if any(
                    [
                        True
                        if isinstance(spin, list) or isinstance(spin, np.ndarray)
                        else False
                        for spin in self.structure.get_initial_magnetic_moments()
                    ]
                ):
                    raise ValueError(
                        "Sphinx only supports collinear spins at the moment."
                    )
                else:
                    spins_str = []
                    for spin, value in zip(
                        self.structure.spin_constraint[self.id_pyi_to_spx],
                        self.structure.get_initial_magnetic_moments()[
                            self.id_pyi_to_spx
                        ],
                    ):
                        if spin:
                            spins_str.append(str(value))
                        else:
                            spins_str.append("X")
                    spins_str = "\n".join(spins_str)
        if spins_str is not None:
            if cwd is not None:
                file_name = posixpath.join(cwd, file_name)
            with open(file_name, "w") as f:
                f.write(spins_str)
        else:
            s.logger.debug("No magnetic moments")


class Input(GenericParameters):
    """
    class to control the generic input for a Sphinx calculation.

    Args:
        input_file_name (str): name of the input file
        table_name (str): name of the GenericParameters table
    """

    def __init__(self, input_file_name=None, table_name="input"):
        super(Input, self).__init__(
            input_file_name=input_file_name,
            table_name=table_name,
            comment_char="//",
            separator_char="=",
            end_value_char=";",
        )
        self._bool_dict = {True: "true", False: "false"}

    def load_default(self):
        """
        Loads the default file content
        """
        file_content = (
            "EnCut = 340\n"
            "KpointCoords = [0.5, 0.5, 0.5]\n"
            "KpointFolding = [4,4,4]\n"
            "EmptyStates = auto\n"
            "Sigma = 0.2\n"
            "Xcorr = PBE\n"
            "PotType = JTH\n"
            "Estep = 400\n"
            "Ediff = 1.0e-4\n"
            "WriteWaves = True\n"
            "KJxc = False\n"
            "SaveMemory = True\n"
            "CoarseRun = False\n"
            "rhoMixing = 1.0\n"
            "spinMixing = 1.0\n"
            "CheckOverlap = True\n"
            "THREADS = 1\n"
        )
        self.load_string(file_content)


class Output(object):
    """
    Handles the output from a Sphinx simulation.
    """

    def __init__(self, job):
        self._job = job
        self._parse_dict = defaultdict(list)

    @staticmethod
    def splitter(arr, counter):
        if len(arr) == 0 or len(counter) == 0:
            return []
        arr_new = []
        spl_loc = list(np.where(np.array(counter) == min(counter))[0])
        spl_loc.append(None)
        for ii, ll in enumerate(spl_loc[:-1]):
            arr_new.append(np.array(arr[ll : spl_loc[ii + 1]]).tolist())
        return arr_new

    def collect_spins_dat(self, file_name="spins.dat", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        if not os.path.isfile(posixpath.join(cwd, file_name)):
            return None
        spins = np.loadtxt(posixpath.join(cwd, file_name))
        self._parse_dict["atom_scf_spins"] = self.splitter(
            np.array([ss[self._job.id_spx_to_pyi] for ss in spins[:, 1:]]), spins[:, 0]
        )

    def collect_energy_dat(self, file_name="energy.dat", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        if not os.path.isfile(posixpath.join(cwd, file_name)):
            return None
        energies = np.loadtxt(posixpath.join(cwd, file_name))
        self._parse_dict["scf_computation_time"] = self.splitter(
            energies[:, 1], energies[:, 0]
        )
        self._parse_dict["scf_energy_int"] = self.splitter(
            energies[:, 2] * HARTREE_TO_EV, energies[:, 0]
        )
        if len(energies[0]) == 7:
            self._parse_dict["scf_energy_free"] = self.splitter(
                energies[:, 3] * HARTREE_TO_EV, energies[:, 0]
            )
            self._parse_dict["scf_energy_zero"] = self.splitter(
                energies[:, 4] * HARTREE_TO_EV, energies[:, 0]
            )
            self._parse_dict["scf_energy_band"] = self.splitter(
                energies[:, 5] * HARTREE_TO_EV, energies[:, 0]
            )
            self._parse_dict["scf_electronic_entropy"] = self.splitter(
                energies[:, 6] * HARTREE_TO_EV, energies[:, 0]
            )
        else:
            self._parse_dict["scf_energy_band"] = self.splitter(
                energies[:, 3] * HARTREE_TO_EV, energies[:, 0]
            )

    def collect_residue_dat(self, file_name="residue.dat", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        if not os.path.isfile(posixpath.join(cwd, file_name)):
            return None
        residue = np.loadtxt(posixpath.join(cwd, file_name))
        if len(residue) == 0:
            return None
        if len(residue[0]) == 2:
            self._parse_dict["scf_residue"] = self.splitter(
                residue[:, 1] * HARTREE_TO_EV, residue[:, 0]
            )
        else:
            self._parse_dict["scf_residue"] = self.splitter(
                residue[:, 1:] * HARTREE_TO_EV, residue[:, 0]
            )

    def collect_eps_dat(self, file_name="eps.dat", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = posixpath.join(cwd, file_name)
        if len(self._parse_dict["bands_eigen_values"]) != 0:
            return None
        if os.path.isfile(file_name):
            try:
                self._parse_dict["bands_eigen_values"] = np.loadtxt(file_name)[:, 1:]
            except:
                self._parse_dict["bands_eigen_values"] = np.loadtxt(file_name)[1:]
        else:
            if os.path.isfile(posixpath.join(cwd, "eps.0.dat")) and os.path.isfile(
                posixpath.join(cwd, "eps.1.dat")
            ):
                eps_up = np.loadtxt(posixpath.join(cwd, "eps.0.dat"))
                eps_down = np.loadtxt(posixpath.join(cwd, "eps.1.dat"))
                if len(eps_up.shape) == 2:
                    eps_up = eps_up[:, 1:]
                    eps_down = eps_down[:, 1:]
                else:
                    eps_up = eps_up[1:]
                    eps_down = eps_down[1:]
                self._parse_dict["bands_eigen_values"] = np.array(
                    list(zip(eps_up.tolist(), eps_down.tolist()))
                )
        return None

    def collect_energy_struct(self, file_name="energy-structOpt.dat", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        energy_free_lst = []
        file_name = posixpath.join(cwd, file_name)
        if os.path.isfile(file_name):
            with open(file_name, "r") as f:
                for line in f.readlines():
                    line = line.split()
                    energy_free_lst.append(float(line[1]) * HARTREE_TO_EV)
        self._energy_free_struct_lst = energy_free_lst

    def collect_sphinx_log(
        self, file_name="sphinx.log", cwd=None, check_consistency=True
    ):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """

        if not os.path.isfile(posixpath.join(cwd, file_name)):
            return None

        def check_conv(line):
            if line.startswith("WARNING: Maximum number of steps exceeded"):
                return False
            elif line.startswith("Convergence reached"):
                return True
            else:
                return None

        with open(posixpath.join(cwd, file_name), "r") as sphinx_log_file:
            log_file = sphinx_log_file.readlines()
            if not np.any(["Enter Main Loop" in line for line in log_file]):
                self._job.status.aborted = True
                raise AssertionError("SPHInX did not enter the main loop; output not collected")
            if not np.any(["Program exited normally." in line for line in log_file]):
                self._job.status.aborted = True
                warnings.warn("SPHInX parsing failed; most likely SPHInX crashed.")
            main_start = np.where(["Enter Main Loop" in line for line in log_file])[0][
                0
            ]
            log_main = log_file[main_start:]

            def get_partial_log(file_content, start_line, end_line):
                start_line = np.where([line == start_line for line in file_content])[0][
                    0
                ]
                end_line = np.where(
                    [line == end_line for line in file_content[start_line:]]
                )[0][0]
                return file_content[start_line : start_line + end_line]
            k_points = get_partial_log(
                log_file,
                "| Symmetrized k-points:               in cartesian coordinates\n",
                "\n",
            )[2:-1]
            self._parse_dict["bands_k_weights"] = np.array([float(kk.split()[6]) for kk in k_points])
            k_points = (
                np.array(
                    [[float(kk.split()[i]) for i in range(2, 5)] for kk in k_points]
                )
                / BOHR_TO_ANGSTROM
            )
            counter = [
                int(line.replace("F(", "").replace(")", " ").split()[0])
                for line in log_main
                if line.startswith("F(")
            ]
            energy_free = [
                float(line.replace("=", " ").replace(",", " ").split()[1])
                * HARTREE_TO_EV
                for line in log_main
                if line.startswith("F(")
            ]
            energy_int = [
                float(line.replace("=", " ").replace(",", " ").split()[1])
                * HARTREE_TO_EV
                for line in log_main
                if line.startswith("eTot(") and not line.startswith("eTot(Val)")
            ]
            energy_zero = 0.5 * (np.array(energy_free) + np.array(energy_int))
            energy_band = [
                float(line.split()[2]) * HARTREE_TO_EV
                for line in log_main
                if line.startswith("eBand")
            ]

            forces = [
                float(re.split("{|}", line)[1].split(",")[i])
                * HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM
                for line in log_main
                for i in range(3)
                if line.startswith("Species: ")
            ]
            magnetic_forces = [
                HARTREE_TO_EV * float(line.split()[-1])
                for line in log_main
                if line.startswith("nu(")
            ]
            convergence = [
                check_conv(line) for line in log_main if check_conv(line) is not None
            ]
            self._parse_dict["bands_e_fermi"] = np.array(
                [
                    float(line.split()[3])
                    for line in log_main
                    if line.startswith("| Fermi energy:")
                ]
            )
            line_vol = np.where(["Omega:" in line for line in log_file])[0][0]
            volume = float(log_file[line_vol].split()[2]) * BOHR_TO_ANGSTROM ** 3
            self._parse_dict["bands_occ"] = [
                line.split()[3:]
                for line in log_main
                if line.startswith("| final focc:")
            ]
            self._parse_dict["bands_occ_initial"] = [
                line.split()[3:] for line in log_main if line.startswith("| focc:")
            ]
            self._parse_dict["bands_eigen_values"] = [
                line.split()[4:]
                for line in log_main
                if line.startswith("| final eig [eV]:")
            ]
            self._parse_dict["bands_eigen_values_initial"] = [
                line.split()[4:] for line in log_main if line.startswith("| eig [eV]:")
            ]

            def eig_converter(
                arr,
                magnetic=np.any(
                    self._job.structure.get_initial_magnetic_moments() != None
                ),
                len_k_points=len(k_points),
            ):
                if len(arr) == 0:
                    return np.array([])
                elif magnetic:
                    return np.array([float(ff) for f in arr for ff in f]).reshape(
                        -1, 2, len_k_points, len(arr[0])
                    )
                else:
                    return np.array([float(ff) for f in arr for ff in f]).reshape(
                        -1, len_k_points, len(arr[0])
                    )

            self._parse_dict["bands_occ"] = eig_converter(self._parse_dict["bands_occ"])
            self._parse_dict["bands_occ_initial"] = eig_converter(
                self._parse_dict["bands_occ_initial"]
            )
            self._parse_dict["bands_eigen_values"] = eig_converter(
                self._parse_dict["bands_eigen_values"]
            )
            self._parse_dict["bands_eigen_values_initial"] = eig_converter(
                self._parse_dict["bands_eigen_values_initial"]
            )
            energy_free_lst = self.splitter(energy_free, counter)
            energy_int_lst = self.splitter(energy_int, counter)
            energy_zero_lst = self.splitter(energy_zero, counter)
            energy_band_lst = self.splitter(energy_band, counter)
            if len(forces) != 0:
                forces = np.array(forces).reshape(-1, len(self._job.structure), 3)
                for ii, ff in enumerate(forces):
                    forces[ii] = ff[self._job.id_spx_to_pyi]
            if len(magnetic_forces) != 0:
                magnetic_forces = np.array(magnetic_forces).reshape(
                    -1, len(self._job.structure)
                )
                for ii, mm in enumerate(magnetic_forces):
                    magnetic_forces[ii] = mm[self._job.id_spx_to_pyi]
                magnetic_forces = self.splitter(magnetic_forces, counter)
        if len(convergence) == len(energy_free_lst) - 1:
            convergence.append(False)
        self._parse_dict["scf_convergence"] = convergence
        self._parse_dict["volume"] = np.array(len(convergence) * [volume])
        if len(self._parse_dict["scf_energy_int"]) == 0 and len(energy_int_lst) != 0:
            self._parse_dict["scf_energy_int"] = energy_int_lst
        if len(self._parse_dict["scf_energy_free"]) == 0 and len(energy_free_lst) != 0:
            self._parse_dict["scf_energy_free"] = energy_free_lst
        if len(self._parse_dict["forces"]) == 0 and len(forces) != 0:
            self._parse_dict["forces"] = forces
        if len(self._parse_dict["scf_magnetic_forces"]) == 0 and len(magnetic_forces) != 0:
            self._parse_dict["scf_magnetic_forces"] = magnetic_forces

    def collect_relaxed_hist(self, file_name="relaxHist.sx", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = posixpath.join(cwd, file_name)
        if not os.path.isfile(file_name):
            return None
        with open(file_name, "r") as file_content:
            file_content = file_content.readlines()
            natoms = len(self._job.id_spx_to_pyi)
            coords = np.array(
                [
                    json.loads(line.split("=")[1].split(";")[0])
                    for line in file_content
                    if "coords" in line
                ]
            )
            self._parse_dict["positions"] = (
                coords.reshape(-1, natoms, 3) * BOHR_TO_ANGSTROM
            )
            self._parse_dict["positions"] = np.array(
                [ff[self._job.id_spx_to_pyi] for ff in self._parse_dict["positions"]]
            )
            force = np.array(
                [
                    json.loads(line.split("=")[1].split(";")[0])
                    for line in file_content
                    if "force" in line
                ]
            )
            self._parse_dict["forces"] = (
                force.reshape(-1, natoms, 3) * HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM
            )
            self._parse_dict["forces"] = np.array(
                [ff[self._job.id_spx_to_pyi] for ff in self._parse_dict["forces"]]
            )
            self._parse_dict["cell"] = (
                np.array(
                    [
                        json.loads(
                            "".join(file_content[i_line : i_line + 3])
                            .split("=")[1]
                            .split(";")[0]
                        )
                        for i_line, line in enumerate(file_content)
                        if "cell" in line
                    ]
                )
                * BOHR_TO_ANGSTROM
            )

    def collect(self, directory=os.getcwd()):
        """
        The collect function, collects all the output from a Sphinx simulation.

        Args:
            directory (str): the directory to collect the output from.
        """
        self.collect_energy_dat(file_name="energy.dat", cwd=directory)
        self.collect_residue_dat(file_name="residue.dat", cwd=directory)
        self.collect_eps_dat(file_name="eps.dat", cwd=directory)
        self.collect_spins_dat(file_name="spins.dat", cwd=directory)
        self.collect_energy_struct(file_name="energy-structOpt.dat", cwd=directory)
        self.collect_sphinx_log(file_name="sphinx.log", cwd=directory)
        self.collect_relaxed_hist(file_name="relaxHist.sx", cwd=directory)
        self._job.compress()

    def to_hdf(self, hdf, force_update=False):
        """
        Store output in an HDF5 file

        Args:
            hdf: HDF5 group
            force_update(bool):
        """

        if len(self._parse_dict["scf_energy_zero"]) == 0:
            self._parse_dict["scf_energy_zero"] = [
                (0.5 * (np.array(fr) + np.array(en))).tolist()
                for fr, en in zip(
                    self._parse_dict["scf_energy_free"],
                    self._parse_dict["scf_energy_int"],
                )
            ]
        with hdf.open("input") as hdf5_input:
            with hdf5_input.open("generic") as hdf5_generic:
                if "dft" not in hdf5_generic.list_groups():
                    hdf5_generic.create_group("dft")
                with hdf5_generic.open("dft") as hdf5_dft:
                    if (
                        len(self._parse_dict["atom_spin_constrains"]) > 0
                        and "atom_spin_constraints" not in hdf5_dft.list_nodes()
                    ):
                        hdf5_dft["atom_spin_constraints"] = [
                            self._parse_dict["atom_spin_constrains"]
                        ]

        with hdf.open("output") as hdf5_output:
            if "sphinx" not in hdf5_output.list_groups():
                hdf5_output.create_group("sphinx")
            with hdf5_output.open("sphinx") as hdf5_sphinx:
                hdf5_sphinx["bands_occ_initial"] = self._parse_dict["bands_occ_initial"]
                hdf5_sphinx["bands_eigen_values_initial"] = self._parse_dict[
                    "bands_eigen_values_initial"
                ]
            with hdf5_output.open("generic") as hdf5_generic:
                if "dft" not in hdf5_generic.list_groups():
                    hdf5_generic.create_group("dft")
                with hdf5_generic.open("dft") as hdf5_dft:
                    hdf5_dft["scf_convergence"] = self._parse_dict["scf_convergence"]
                    for k in [
                        "scf_residue",
                        "scf_energy_free",
                        "scf_energy_zero",
                        "scf_energy_int",
                        "scf_electronic_entropy",
                        "scf_energy_band",
                        "scf_magnetic_forces",
                        "scf_computation_time",
                        "bands_occ",
                        "bands_e_fermi",
                        "bands_k_weights",
                        "bands_eigen_values",
                        "atom_scf_spins",
                    ]:
                        if len(self._parse_dict[k]) > 0:
                            hdf5_dft[k] = self._parse_dict[k]
                            if "scf_" in k:
                                hdf5_dft[k.replace("scf_", "")] = np.array(
                                    [vv[-1] for vv in self._parse_dict[k]]
                                )
                if len(self._parse_dict["scf_computation_time"]) > 0:
                    hdf5_generic["computation_time"] = np.array(
                        [tt[-1] for tt in self._parse_dict["scf_computation_time"]]
                    )
                if len([en[-1] for en in self._parse_dict["scf_energy_free"]]) > 0:
                    hdf5_generic["energy_tot"] = np.array(
                        [en[-1] for en in self._parse_dict["scf_energy_free"]]
                    )
                    hdf5_generic["energy_pot"] = np.array(
                        [en[-1] for en in self._parse_dict["scf_energy_free"]]
                    )
                hdf5_generic["volume"] = self._parse_dict["volume"]
                if "positions" not in hdf5_generic.list_nodes() or force_update:
                    if len(self._parse_dict["positions"]) > 0:
                        hdf5_generic["positions"] = np.array(
                            self._parse_dict["positions"]
                        )
                    elif len(self._parse_dict["scf_convergence"]) == 1:
                        hdf5_generic["positions"] = np.array(
                            [self._job.structure.positions]
                        )
                if ("forces" not in hdf5_generic.list_nodes() or force_update) and len(
                    self._parse_dict["forces"]
                ) > 0:
                    hdf5_generic["forces"] = np.array(self._parse_dict["forces"])
                if "cells" not in hdf5_generic.list_nodes() or force_update:
                    if len(self._parse_dict["cell"]) > 0:
                        hdf5_generic["cells"] = np.array(self._parse_dict["cell"])
                    elif len(self._parse_dict["scf_convergence"]) == 1:
                        hdf5_generic["cells"] = np.array([self._job.structure.cell])

    def from_hdf(self, hdf):
        """
        Load output from an HDF5 file
        """
