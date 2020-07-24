# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH -Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, division

import numpy as np
import os
import posixpath
import re
import stat
from shutil import copyfile
import scipy.constants
import warnings
import json
from collections import defaultdict

from pyiron.dft.job.generic import GenericDFTJob
from pyiron.vasp.potential import VaspPotentialFile
from pyiron.vasp.potential import find_potential_file \
    as find_potential_file_vasp
from pyiron.sphinx.structure import read_atoms
from pyiron.sphinx.potential import SphinxJTHPotentialFile
from pyiron.sphinx.potential import find_potential_file \
    as find_potential_file_jth
from pyiron.sphinx.volumetric_data import SphinxVolumetricData
from pyiron.base.settings.generic import Settings
from pyiron.base.generic.inputlist import InputList

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
    scipy.constants.physical_constants["Bohr radius"][0] /
    scipy.constants.angstrom
)
HARTREE_TO_EV = scipy.constants.physical_constants["Hartree energy in eV"][0]
RYDBERG_TO_EV = HARTREE_TO_EV / 2
HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM = HARTREE_TO_EV / BOHR_TO_ANGSTROM


class SphinxBase(GenericDFTJob):
    """
    Class to setup and run Sphinx simulations.

    Inherits pyiron_atomistics.job.generic.GenericJob. The functions in
    these modules are written such that the function names and attributes
    are very Pyiron-generic (get_structure(), molecular_dynamics(),
    version) but internally handle Sphinx specific input and output.

    Alternatively, because SPHInX inputs are built on a group-based
    format, users have the option to set specific groups and parameters
    directly, e.g.

    ```python
    # Modify/add a new parameter via
    job.input.PAWHamiltonian.nEmptyStates = 15
    job.input.PAWHamiltonian.dipoleCorrection = True
    # or
    job.input.PAWHamiltonian.set("nEmptyStates", 15)
    job.input.PAWHamiltonian.set("dipoleCorrection", True)
    # Modify/add a sub-group via
    job.input.initialGuess.rho.charged = {"charge": 2, "z": 25}
    # or
    job.input.initialGuess.rho.set("charged", {"charge": 2, "z": 25})
    ```

    Args:
        project: Project object (defines path where job will be
                 created and stored)
        job_name (str): name of the job (must be unique within
                        this project path)
    """

    """Version of the data format in hdf5"""
    __hdf_version__ = "0.1.0"

    def __init__(self, project, job_name):
        super(SphinxBase, self).__init__(project, job_name)

        # keeps both the generic parameters as well as the sphinx specific
        # input groups
        self.input = Group(table_name = "parameters")
        self.load_default_input()
        self._save_memory = False
        self._output_parser = Output(self)
        self.input_writer = InputWriter()
        if self.check_vasp_potentials():
            self.input["VaspPot"] = True  # use VASP potentials if available
        self._generic_input["restart_for_band_structure"] = False
        self._generic_input["path_name"] = None
        self._generic_input["n_path"] = None

    def get_version_float(self):
        version_str = self.executable.version.split("_")[0]
        version_float = float(
            version_str.split(".")[0]
        )
        if len(version_str.split(".")) > 1:
            version_float += float(
                "0." + "".join(version_str.split(".")[1:])
                )
        return version_float

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
        if "eCut" in self.input.sphinx.basis.keys():
            return self.input.sphinx.basis["eCut"] * RYDBERG_TO_EV
        else:
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
        elif val < 50:
            warnings.warn(
                "The given cutoff is either very small (probably "
                + "too small) or was accidentally given in Ry. "
                + "Please make sure it is in eV (1eV = 13.606 Ry)."
            )
        self.input["EnCut"] = val
        self.input.sphinx.basis.eCut = self.input["EnCut"] / RYDBERG_TO_EV

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
                "Exchange correlation function not recognized (\
                    recommended: PBE or LDA)",
                SyntaxWarning,
            )
            self.input["Xcorr"] = val
        if "xc" in self.input.sphinx.PAWHamiltonian.keys():
            self.input.sphinx.PAWHamiltonian.xc = self.input["Xcorr"]

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes,
        but it has to be implemented in the individual classes.
        """
        super(SphinxBase, self).set_input_to_read_only()
        self.input.read_only = True

    def get_scf_group(
        self, maxSteps=None, keepRhoFixed=False, dEnergy=None,
        algorithm="blockCCG"
    ):
        """
        SCF group setting for SPHInX
        for all args refer to calc_static or calc_minimize
        """

        scf_group = Group()
        if algorithm.upper() == "CCG":
            algorithm = "CCG"
        elif algorithm.upper() != "BLOCKCCG":
            warnings.warn(
                "Algorithm not recognized -> setting to blockCCG. \
                    Alternatively, choose algorithm=CCG",
                SyntaxWarning,
            )
            algorithm = "blockCCG"

        if keepRhoFixed:
            scf_group["keepRhoFixed"] = True
        else:
            scf_group["rhoMixing"] = str(self.input["rhoMixing"])
            scf_group["spinMixing"] = str(self.input["spinMixing"])
            if "nPulaySteps" in self.input:
                scf_group["nPulaySteps"] = str(self.input["nPulaySteps"])
        if dEnergy is None:
            scf_group["dEnergy"] = self.input["Ediff"] / HARTREE_TO_EV
        else:
            scf_group["dEnergy"] = str(dEnergy)
        if maxSteps is None:
            scf_group["maxSteps"] = str(self.input["Estep"])
        else:
            scf_group["maxSteps"] = str(maxSteps)
        if "preconditioner" in self.input:
            scf_group.create_group("preconditioner")["type"] = \
                    self.input["preconditioner"]
        scf_group.create_group(algorithm)
        if "maxStepsCCG" in self.input:
            scf_group[algorithm]["maxStepsCCG"] = self.input["maxStepsCCG"]
        if "blockSize" in self.input and algorithm == "blockCCG":
            scf_group[algorithm]["blockSize"] = self.input["blockSize"]
        if "nSloppy" in self.input and algorithm == "blockCCG":
            scf_group[algorithm]["nSloppy"] = self.input["nSloppy"]
        if self.input["WriteWaves"] is False:
            scf_group["noWavesStorage"] = True
        return scf_group

    def get_structure_group(self, keep_angstrom=False):
        """
        create a Sphinx Group object based on self.structure

        Args:
            keep_angstrom (bool): Store distances in Angstroms or Bohr
        """
        if keep_angstrom:
            structure_group = Group( {"cell": np.array(self.structure.cell)} )
        else:
            structure_group = Group( {
                    "cell": np.array(self.structure.cell * 1 / BOHR_TO_ANGSTROM),
            })
        if "selective_dynamics" in self.structure._tag_list.keys():
            selective_dynamics_list = \
                self.structure.selective_dynamics.list()
        else:
            selective_dynamics_list = [3 * [False]] * len(
                self.structure.positions)
        species = structure_group.create_group("species")
        for elm_species in self.structure.get_species_objects():
            if elm_species.Parent:
                element = elm_species.Parent
            else:
                element = elm_species.Abbreviation
            species.append(
                Group({"element": '"' + str(element) + '"'})
            )
            elm_list = np.array(
                self.structure.get_chemical_symbols() == \
                    elm_species.Abbreviation
            )
            atom_group = species[-1].create_group("atom")
            for elm_pos, elm_magmon, selective in zip(
                self.structure.positions[elm_list],
                np.array(self.structure.get_initial_magnetic_moments())[
                    elm_list],
                np.array(selective_dynamics_list)[elm_list],
            ):
                atom_group.append(Group())
                if self._spin_enabled:
                    atom_group[-1]["label"] \
                        = '"spin_' + str(elm_magmon) + '"'
                if keep_angstrom:
                    atom_group[-1]["coords"] = np.array(elm_pos)
                else:
                    atom_group[-1]["coords"] = \
                        np.array(elm_pos * 1 / BOHR_TO_ANGSTROM)
                if all(selective):
                    atom_group[-1]["movable"] = True
                elif any(selective):
                    for ss, xx in zip(selective, ["X", "Y", "Z"]):
                        if ss:
                            atom_group[-1]["movable" + xx] = True
        if not self.fix_symmetry:
            structure_group.symmetry = Group({
                "operator": {
                    "S": "[[1,0,0],[0,1,0],[0,0,1]]"
            } })
        return structure_group

    def load_default_input(self):
        """
        Set defaults for generic parameters and create sphinx input groups.
        """

        sph = self.input.create_group('sphinx')
        sph.create_group('pawPot')
        sph.create_group('structure')
        sph.create_group('basis')
        sph.create_group('PAWHamiltonian')
        sph.create_group('initialGuess')
        sph.create_group('main')

        self.input.EnCut = 340
        self.input.KpointCoords = [0.5, 0.5, 0.5]
        self.input.KpointFolding = [4,4,4]
        self.input.EmptyStates = "auto"
        self.input.Sigma = 0.2
        self.input.Xcorr = "PBE"
        self.input.VaspPot = False
        self.input.Estep = 400
        self.input.Ediff = 1.0e-4
        self.input.WriteWaves = True
        self.input.KJxc = False
        self.input.SaveMemory = True
        self.input.CoarseRun = False
        self.input.rhoMixing = 1.0
        self.input.spinMixing = 1.0
        self.input.CheckOverlap = True
        self.input.THREADS = 1

    def load_structure_group(self, keep_angstrom=False):
        """
        Build + load the structure group based on self.structure

        Args:
            keep_angstrom (bool): Store distances in Angstroms or Bohr
        """
        self.input.sphinx.structure = self.get_structure_group(
            keep_angstrom=keep_angstrom
            )

    def load_species_group(self, check_overlap=True, potformat='VASP'):
        """
        Build the species Group object based on self.structure

        Args:
            check_overlap (bool): Whether to check overlap
                (see set_check_overlap)
            potformat (str): type of pseudopotentials that will be
                read. Options are JTH or VASP.
        """

        self.input.sphinx.pawPot = Group({"species": []})
        for species_obj in self.structure.get_species_objects():
            if species_obj.Parent is not None:
                elem = species_obj.Parent
            else:
                elem = species_obj.Abbreviation
            if potformat == "JTH":
                self.input.sphinx.pawPot["species"].append({
                            "name": '"' + elem + '"',
                            "potType": '"AtomPAW"',
                            "element": '"' + elem + '"',
                            "potential": f'"{elem}_GGA.atomicdata"',
                })
            elif potformat == "VASP":
                self.input.sphinx.pawPot["species"].append({
                            "name": '"' + elem + '"',
                            "potType": '"VASP"',
                            "element": '"' + elem + '"',
                            "potential": '"' + elem + "_POTCAR" + '"',
                })
            else:
                raise ValueError()
        if not check_overlap:
            self.input.sphinx.pawPot["species"][-1]["checkOverlap"] = "false"
        if self.input["KJxc"]:
            self.input.sphinx.pawPot["kjxc"] = True

    def load_main_group(self):
        """
        Load the main Group.

        The group is populated based on the type of calculation and settings in
        the self.input.
        """

        if len(self.restart_file_list) != 0 \
        and not self._generic_input["restart_for_band_structure"]:
            self.input.sphinx.main.get("scfDiag", create = True).append(
                self.get_scf_group(
                    maxSteps=10, keepRhoFixed=True, dEnergy=1.0e-4
                )
            )
        if "Istep" in self.input:
            self.input.sphinx.main["ricQN"] = Group(table_name = "input")
            self.input.sphinx.main["ricQN"]["maxSteps"] = str(self.input["Istep"])
            if "dE" in self.input and "dF" in self.input:
                self.input["dE"] = 1e-3
            if "dE" in self.input:
                self.input.sphinx.main["ricQN"]["dEnergy"] = str(
                    self.input["dE"] / HARTREE_TO_EV
                    )
            if "dF" in self.input:
                self.input.sphinx.main["ricQN"]["dF"] = str(
                    self.input["dF"] / HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM
                )
            self.input.sphinx.main.ricQN.create_group("bornOppenheimer")
            self.input.sphinx.main.ricQN.bornOppenheimer["scfDiag"] = \
                self.get_scf_group()
        else:
            scf = self.input.sphinx.main.get("scfDiag", create = True)
            if self._generic_input["restart_for_band_structure"]:
                scf.append(
                    self.get_scf_group(keepRhoFixed=True)
                    )
            else:
                scf.append(self.get_scf_group())
            if self.executable.version is not None:
                vers_num = [
                    int(vv)
                    for vv in self.executable.version.split("_")[0].split(".")
                ]
                if self.get_version_float() > 2.5:
                    self.input.sphinx.main.create_group("evalForces")["file"] = \
                            '"relaxHist.sx"'
            else:
                warnings.warn("executable version could not be identified")

    def load_basis_group(self):
        """
        Load the basis Group.

        The group is populated using setdefault to avoid
        overwriting values that were previously (intentionally)
        modified.
        """
        self.input.sphinx.basis.setdefault("eCut", self.input["EnCut"]/RYDBERG_TO_EV)
        self.input.sphinx.basis.get("kPoint", create = True)
        if "KpointCoords" in self.input:
            self.input.sphinx.basis.kPoint.setdefault("coords",
                    np.array(self.input["KpointCoords"]))
        self.input.sphinx.basis.kPoint.setdefault("weight", 1)
        self.input.sphinx.basis.kPoint.setdefault("relative", True)
        if "KpointFolding" in self.input:
            self.input.sphinx.basis.setdefault("folding",
                    np.array(self.input["KpointFolding"]))
        self.input.sphinx.basis.setdefault("saveMemory",
                self.input["SaveMemory"])

    def load_hamilton_group(self):
        """
        Load the PAWHamiltonian Group.

        The group is populated using setdefault to avoid
        overwriting values that were previously (intentionally)
        modified.
        """
        self.input.sphinx.PAWHamiltonian.setdefault(
            "nEmptyStates", self.input["EmptyStates"]
            )
        self.input.sphinx.PAWHamiltonian.setdefault(
            "ekt", self.input["Sigma"]/HARTREE_TO_EV
            )
        self.input.sphinx.PAWHamiltonian.setdefault("xc", self.input["Xcorr"])
        self.input.sphinx.PAWHamiltonian["spinPolarized"] = self._spin_enabled

    def load_guess_group(self, update_spins=True):
        """
        Load the initialGuess Group.

        The group is populated using setdefault to avoid
        overwriting values that were previously (intentionally)
        modified.

        Args:
            update_spins (bool): whether or not to reload the
                atomicSpin groups based on the latest structure.
                Defaults to True.
        """

        charge_density_file = None
        for ff in self.restart_file_list:
            if "rho.sxb" in ff.split("/")[-1]:
                charge_density_file = ff
        wave_function_file = None
        for ff in self.restart_file_list:
            if "waves.sxb" in ff.split("/")[-1]:
                wave_function_file = ff

        self.input.sphinx.initialGuess.setdefault("waves", Group())
        self.input.sphinx.initialGuess.waves.setdefault("lcao", Group())
        self.input.sphinx.initialGuess.waves.setdefault("pawBasis", True)
        if wave_function_file is not None:
            self.input.sphinx.initialGuess.setdefault("exchange", Group())
            self.input.sphinx.initialGuess.exchange.setdefault(
                "file", '"' + wave_function_file + '"'
            )
        if charge_density_file is None:
            self.input.sphinx.initialGuess.setdefault("rho", Group({"atomicOrbitals": True}))
        else:
            self.input.sphinx.initialGuess.setdefault(
                "rho", Group({"file": '"' + charge_density_file + '"'})
                )
        if self._spin_enabled:
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
                rho = self.input.sphinx.initialGuess.rho
                rho.get("atomicSpin", create = True)
                if update_spins:
                    rho.atomicSpin.clear()
                if len(rho.atomicSpin) == 0:
                    for spin in self.structure.get_initial_magnetic_moments()[
                        self.id_pyi_to_spx
                    ]:
                        rho["atomicSpin"].append(
                            Group({
                                "label": '"spin_' + str(spin) + '"',
                                "spin": str(spin)
                            })
                        )

        if "noWavesStorage" not in self.input.sphinx.initialGuess:
            self.input.sphinx.initialGuess["noWavesStorage"] = \
                    not self.input["WriteWaves"]
    def calc_static(
        self,
        electronic_steps=None,
        algorithm=None,
        retain_charge_density=False,
        retain_electrostatic_potential=False,
    ):
        """
        Setup the hamiltonian to perform a static SCF run.

        Loads defaults for all Sphinx input groups, including a static
        main Group.

        Args:
            electronic_steps (float): max # of electronic steps
            retain_electrostatic_potential:
            retain_charge_density:
            algorithm (str): CCG or blockCCG (not implemented)
            electronic_steps (int): maximum number of electronic steps
                which can be used to achieve convergence
        """
        if electronic_steps is not None:
            self.input["Estep"] = electronic_steps
        for arg in ["Istep", "dF", "dE"]:
            if arg in self.input:
                del self.input[arg]
        super(SphinxBase, self).calc_static(
            electronic_steps=electronic_steps,
            algorithm=algorithm,
            retain_charge_density=retain_charge_density,
            retain_electrostatic_potential=retain_electrostatic_potential,
        )
        self.load_default_groups()


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
        ionic_energy_tolerance=0.0,
        ionic_force_tolerance=1.0e-2,
        volume_only=False,
    ):
        """
        Setup the hamiltonian to perform ionic relaxations.

        The convergence goal can be set using either the
        ionic_energy_tolerance as a limit for fluctuations in energy or the
        ionic_force_tolerance.

        Loads defaults for all Sphinx input groups, including a
        ricQN-based main Group.

        Args:
            retain_electrostatic_potential:
            retain_charge_density:
            algorithm:
            pressure:
            max_iter:
            electronic_steps (int): maximum number of electronic steps
                                    per electronic convergence
            ionic_steps (int): maximum number of ionic steps
            ionic_energy_tolerance (float): convergence goal in terms of
                                  energy (optional)
            ionic_force_tolerance (float): convergence goal in terms of
                                  forces (optional)
        """
        if pressure is not None:
            raise NotImplementedError(
                "pressure minimization is not implemented in SPHInX"
            )
        if electronic_steps is not None:
            self.input["Estep"] = electronic_steps
        if ionic_steps is not None:
            self.input["Istep"] = ionic_steps
        elif "Istep" not in self.input:
            self.input["Istep"] = 100
        if ionic_force_tolerance is not None:
            if ionic_force_tolerance < 0:
                raise ValueError("ionic_force_tolerance must be a positive integer")
            self.input["dF"] = float(ionic_force_tolerance)
        if ionic_energy_tolerance is not None:
            if ionic_energy_tolerance < 0:
                raise ValueError("ionic_force_tolerance must be a positive integer")
            self.input["dE"] = float(ionic_energy_tolerance)
        super(SphinxBase, self).calc_minimize(
            electronic_steps=electronic_steps,
            ionic_steps=ionic_steps,
            max_iter=max_iter,
            pressure=pressure,
            algorithm=algorithm,
            retain_charge_density=retain_charge_density,
            retain_electrostatic_potential=retain_electrostatic_potential,
            ionic_energy_tolerance=ionic_energy_tolerance,
            ionic_force_tolerance=ionic_force_tolerance,
            volume_only=volume_only,
        )
        self.load_default_groups()


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

    def restart_for_band_structure_calculations(self, job_name=None):
        """
        Restart a new job created from an existing calculation
        by reading the charge density for band structures.

        Args:
            job_name (str/None): Job name

        Returns:
            pyiron.sphinx.sphinx.sphinx: new job instance
        """
        return self.restart_from_charge_density(
            job_name=job_name,
            job_type=None,
            band_structure_calc=True
        )

    def restart_from_charge_density(
            self,
            job_name=None,
            job_type="Sphinx",
            band_structure_calc=False
    ):
        """
        Restart a new job created from an existing calculation
        by reading the charge density.

        Args:
            job_name (str/None): Job name
            job_type (str/None): Job type. If not specified a Sphinx
                                 job type is assumed (actually this is
                                 all that's currently supported)
            band_structure_calc (bool): has to be True for band
                                        structure calculations.

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
            job_type="Sphinx",
    ):
        """
        Restart a new job created from an existing calculation
        by reading the wave functions.

        Args:
            job_name (str): Job name
            job_type (str): Job type. If not specified a Sphinx
                            job type is assumed (actually this is
                            all that's currently supported.)

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

        new_job.input = self.input

        if from_charge_density and os.path.isfile(
            posixpath.join(self.working_directory, "rho.sxb")
        ):
            new_job.restart_file_list.append(posixpath.join(self.working_directory, "rho.sxb"))

        elif from_charge_density:
            self._logger.warning(
                msg=f"A charge density from job: {self.job_name} "
                + "is not generated and therefore it can't be read."
            )
        if from_wave_functions and os.path.isfile(
            posixpath.join(self.working_directory, "waves.sxb")
        ):
            new_job.restart_file_list.append(posixpath.join(self.working_directory, "waves.sxb"))
        elif from_wave_functions:
            self._logger.warning(
                msg="No wavefunction file (waves.sxb) was found for "
                + f"job {self.job_name} in {self.working_directory}."
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
        with self._hdf5.open("input") as hdf:
            self.input.to_hdf(hdf)
        self._output_parser.to_hdf(self._hdf5)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Recreates instance from the hdf5 file

        Args:
            hdf (str): Path to the hdf5 file
            group_name (str): Name of the group which contains the object
        """


        if "HDF_VERSION" not in self._hdf5.keys():
            from pyiron.base.generic.parameters import GenericParameters
            super(SphinxBase, self).from_hdf(hdf=hdf, group_name=group_name)
            self._structure_from_hdf()
            gp = GenericParameters(table_name = "input")
            gp.from_hdf(self._hdf5)
            for k in gp.keys():
                self.input[k] = gp[k]
            if self.status.finished:
                self._output_parser.from_hdf(self._hdf5)
        elif self._hdf5["HDF_VERSION"] == "0.1.0":
            super(SphinxBase, self).from_hdf(hdf=hdf, group_name=group_name)
            self._structure_from_hdf()
            with self._hdf5.open("input") as hdf:
                self.input.from_hdf(hdf, group_name = "parameters")
            if self.status.finished:
                self._output_parser.from_hdf(self._hdf5)


    def from_directory(self, directory, file_name="structure.sx"):
        try:
            if not self.status.finished:

                file_path = posixpath.join(directory, file_name)
                if os.path.isfile(file_path):
                    self.structure = read_atoms(file_path)
                else:
                    raise ValueError(f"File {file_path} not found. "
                    + "Please double check the directory and file name.")

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
            Certain PAW-pseudo-potentials have an intrinsic pathology:
            their PAW overlap operator is not generally positive definite
            (i.e., the PAW-corrected norm of a wavefunction could become
            negative). SPHInX usually refuses to use such problematic
            potentials. This behavior can be overridden by setting
            check_overlap to False.
        """
        if not isinstance(check_overlap, bool):
            raise ValueError("check_overlap has to be a boolean")

        if self.get_version_float() < 2.51 and not check_overlap:
            warnings.warn(
                "SPHInX executable version has to be 2.5.1 or above "
                + "in order for the overlap to be considered. "
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
            n_pulay_steps (int): number of previous densities to use for
                                 the Pulay mixing (default: 7)
            density_mixing_parameter (float): mixing proportion m defined by

                rho^n = (m-1)*rho^(n-1)+m*preconditioner*rho_(opt) (default: 1)

            spin_mixing_parameter (float): linear mixing parameter for
                                           spin densities (default: 1)

        comments:
            A low value of density mixing parameter may lead
            to a more stable convergence, but will slow down
            the calculation if set too low.

            Further information can be found on the website:
            https://sxrepo.mpie.de
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
                "density_mixing_parameter has to be between 0 and 1 "+
                "(default value is 1)"
            )
        if spin_mixing_parameter is not None and (
            spin_mixing_parameter > 1.0 or spin_mixing_parameter < 0
        ):
            raise ValueError(
                "spin_mixing_parameter has to be between 0 and 1 "+
                "(default value is 1)"
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
        Set how the finite temperature smearing is applied in
        determining partial occupancies

        Args:
            smearing (str): Type of smearing (only fermi is
                            implemented anything else will be ignored)
            width (float): Smearing width (eV) (default: 0.2)
        """
        if smearing is not None and not isinstance(smearing, str):
            raise ValueError(
                "Smearing must be a string"
            )
        if width is not None and width < 0:
            raise ValueError("Smearing value must be a float >= 0")
        if width is not None:
            self.input["Sigma"] = width

    def set_convergence_precision(
        self, ionic_energy_tolerance=None, ionic_force_tolerance=None, ionic_energy=None, electronic_energy=None, ionic_forces=None
    ):
        """
        Sets the electronic and ionic convergence precision.

        For ionic convergence either the energy or the force
        precision is required.

        Args:
            ionic_energy_tolerance (float): Ionic energy convergence precision
            electronic_energy (float): Electronic energy convergence
                                       precision
            ionic_force_tolerance (float): Ionic force convergence precision
        """
        if ionic_forces is not None:
            warnings.warn(
                "ionic_forces is deprecated as of vers. 0.3.0. It is not guaranteed to be in service in vers. 0.4.0. Use ionic_force_tolerance instead.",
                DeprecationWarning
            )
            ionic_force_tolerance = ionic_forces
        if ionic_energy is not None:
            warnings.warn(
                "ionic_energy is deprecated as of vers. 0.3.0. It is not guaranteed to be in service in vers. 0.4.0. Use ionic_energy_tolerance instead.",
                DeprecationWarning
            )
            ionic_energy_tolerance =ionic_energy
        assert (
            ionic_energy_tolerance is None or ionic_energy_tolerance > 0
        ), "ionic_energy_tolerance must be a positive float"
        assert (
            ionic_force_tolerance is None or ionic_force_tolerance > 0
        ), "ionic_force_tolerance must be a positive float"
        assert (
            electronic_energy is None or electronic_energy > 0
        ), "electronic_energy must be a positive float"
        if ionic_energy_tolerance is not None or ionic_force_tolerance is not None:
            print("Setting calc_minimize")
            self.calc_minimize(ionic_energy_tolerance=ionic_energy_tolerance,
                               ionic_force_tolerance=ionic_force_tolerance)
        if electronic_energy is not None:
            self.input["Ediff"] = electronic_energy

    def set_empty_states(self, n_empty_states=None):
        """
        Function to set the number of empty states.

        Args:
            n_empty_states (int/None): Number of empty states.
            If None, sets it to 'auto'.

        Comments:
            If this number is too low, the algorithm will not be
            able to able to swap wave functions near the chemical
            potential. If the number is too high, computation time
            will be wasted for the higher energy states and
            potentially lead to a memory problem.

            In contrast to VASP, this function sets only the number
            of empty states and not the number of total states.

            The default value is 0.5*NIONS+3 for non-magnetic systems
            and 1.5*NIONS+3 for magnetic systems
        """
        if n_empty_states is None:
            # will be converted later; see load_default_groups
            self.input["EmptyStates"] = "auto"
        else:
            if n_empty_states < 0:
                raise ValueError(
                    "Number of empty states must be greater than 0"
                    )
            self.input["EmptyStates"] = n_empty_states
        self.input.sphinx.PAWHamiltonian.nEmptyStates = self.input["EmptyStates"]

    def _set_kpoints(
        self,
        mesh=None,
        scheme="MP",
        center_shift=None,
        symmetry_reduction=True,
        manual_kpoints=None,
        weights=None,
        reciprocal=True,
        kpoints_per_reciprocal_angstrom=None,
        n_path=None,
        path_name=None,
    ):
        """
        Function to setup the k-points for the Sphinx job

        Args:
            reciprocal (bool): Tells if the supplied values are in
                               reciprocal (direct) or cartesian coordinates
                               (in reciprocal space) (not implemented)
            weights (list): Manually supplied weights to each k-point in
                            case of the manual mode (not implemented)
            manual_kpoints (list): Manual list of k-points (not implemented)
            symmetry_reduction (bool): Tells if the symmetry reduction is
                                       to be applied to the k-points
            scheme (str): Type of k-point generation scheme ('MP' or 'Line')
            mesh (list): Size of the mesh (in the MP scheme)
            center_shift (list): Shifts the center of the mesh from the
                                 gamma point by the given vector
            kpoints_per_reciprocal_angstrom (float): Number of kpoint per angstrom
                                          in each direction
            n_path (int): Number of points per trace part for line mode
            path_name (str): Name of high symmetry path used for band
                             structure calculations.
        """
        if not isinstance(symmetry_reduction, bool):
            raise ValueError("symmetry_reduction has to be a boolean")
        if manual_kpoints is not None:
            raise ValueError("manual_kpoints is not yet implemented in "
                + "Pyiron for SPHInX")
        if weights is not None:
            raise ValueError(
                "manual weights are not yet implmented in Pyiron for "
                + "SPHInX"
                )

        if scheme == "MP":
            # Remove kPoints and set kPoint
            if "kPoints" in self.input.sphinx.basis:
                del self.input.sphinx.basis.kPoints
            if kpoints_per_reciprocal_angstrom is not None:
                if mesh is not None:
                    warnings.warn("mesh value is overwritten "
                    + "by kpoints_per_reciprocal_angstrom")
                mesh = self.get_k_mesh_by_cell(
                    kpoints_per_reciprocal_angstrom=kpoints_per_reciprocal_angstrom
                    )
            self.input.sphinx.basis.get("kPoint", create = True)
            if mesh is not None:
                self.input["KpointFolding"] = list(mesh)
                self.input.sphinx.basis["folding"] = np.array(self.input["KpointFolding"])
            if center_shift is not None:
                self.input["KpointCoords"] = list(center_shift)
                self.input.sphinx.basis["kPoint"]["coords"] = \
                    np.array(self.input["KpointCoords"])
                self.input.sphinx.basis.kPoint["weight"] = 1
                self.input.sphinx.basis.kPoint["relative"] = True

        elif scheme == "Line":
            # Remove Kpoint and set Kpoints

            if "kPoint" in self.input.sphinx.basis:
                del self.input.sphinx.basis["kPoint"]
                del self.input["KpointFolding"]
                del self.input["KpointCoords"]
                if "folding" in self.input.sphinx.basis:
                    del self.input.sphinx.basis['folding']
            if n_path is None and self._generic_input["n_path"] is None:
                raise ValueError("'n_path' has to be defined")
            if n_path is None:
                n_path = self._generic_input["n_path"]
            else:
                self._generic_input["n_path"] = n_path

            if self.structure.get_high_symmetry_points() is None:
                raise ValueError(
                    "no 'high_symmetry_points' defined for 'structure'."
                    )

            if path_name is None and self._generic_input["path_name"] is None:
                raise ValueError("'path_name' has to be defined")
            if path_name is None:
                path_name = self._generic_input["path_name"]
            else:
                self._generic_input["path_name"] = path_name

            try:
                path = self.structure.get_high_symmetry_path()[path_name]
            except KeyError:
                raise AssertionError(
                    "'{}' is not a valid path!".format(path_name)
                    )

            def make_point(point, n_path):
                return Group({"coords":
                    np.array(self.structure.get_high_symmetry_points()[point]),
                    "nPoints": n_path,
                    "label": "\"{}\"".format(point.replace("'", "p"))})

            kpoints = Group({"relative": True})
            kpoints["from"] = make_point(path[0][0], None)
            # from nodes are not supposed to have a nPoints attribute
            del kpoints["from/nPoints"]

            kpoints.create_group("to").append(make_point(path[0][1], n_path))

            for segment in path[1:]:
                # if the last node on the so far is not the same as the first
                # node of this path segment, then we need to insert another
                # node into the path to alert sphinx that we want a cut in our
                # band structure (n_path = 0)
                if '"{}"'.format(segment[0]) != kpoints.to[-1].label:
                    kpoints["to"].append(
                            make_point(segment[0], 0)
                    )

                kpoints["to"].append(make_point(segment[1], n_path))

            self.input.sphinx.basis["kPoints"] = kpoints
        else:
            raise ValueError("only Monkhorst-Pack mesh and Line mode\
                are currently implemented in Pyiron for SPHInX")


    def load_default_groups(self):
        """
        Populates input groups with the default values.

        Nearly every default simply points to a variable stored in
        self.input.

        Does not load job.input.structure or job.input.pawPot.
        These groups should usually be modified via job.structure,
        in which case they will be set at the last minute when
        the job is run. These groups can be synced to job.structure
        at any time using job.load_structure_group() and
        job.load_species_group().
        """

        if self.structure is None:
            raise AssertionError(f"{self.job_name} has not been assigned "
                + "a structure. Please load one first (e.g. "
                + f"{self.job_name}.structure = ...)")

        self._coarse_run = self.input["CoarseRun"]

        if self.input["EmptyStates"] == "auto":
            if self._spin_enabled:
                self.input["EmptyStates"] = int(
                    1.5 * len(self.structure) + 3)
            else:
                self.input["EmptyStates"] = int(len(self.structure) + 3)

        if not self.input.sphinx.basis.locked:
            self.load_basis_group()
        if not self.input.sphinx.structure.locked:
            self.load_structure_group()
        if self.input["VaspPot"]:
            potformat = "VASP"
        else:
            potformat = "JTH"
        if not self.input.sphinx.pawPot.locked:
            self.load_species_group(potformat=potformat)
        if not self.input.sphinx.initialGuess.locked:
            self.load_guess_group()
        if not self.input.sphinx.PAWHamiltonian.locked:
            self.load_hamilton_group()
        if not self.input.sphinx.main.locked:
            self.load_main_group()


    def write_input(self):
        """
        Generate all the required input files for the Sphinx job.

        Creates:
        structure.sx: structure associated w/ job
        all pseudopotential files
        spins.in (if necessary): constrained spin moments
        input.sx: main input file with all sub-groups

        Automatically called by job.run()
        """

        # If the structure group was not modified directly by the
        # user, via job.input.structure (which is likely True),
        # load it based on job.structure.
        structure_sync = (str(self.input.sphinx.structure)
                          == str(self.get_structure_group()))
        if not structure_sync and not self.input.sphinx.structure.locked:
            self.load_structure_group()

        # copy potential files to working directory
        if self.input["VaspPot"]:
            potformat = "VASP"
        else:
            potformat = "JTH"

        # If the species group was not modified directly by the user,
        # via job.input.pawPot (which is likely True),
        # load it based on job.structure.
        if not structure_sync and not self.input.sphinx.pawPot.locked:
            self.load_species_group(potformat=potformat)

        self.input_writer.structure = self.structure
        self.input_writer.copy_potentials(
            potformat=potformat,
            xc=self.input["Xcorr"],
            cwd=self.working_directory
            )

        # Write spin constraints, if set via _generic_input.
        all_groups = [
            self.input.sphinx.pawPot,
            self.input.sphinx.structure,
            self.input.sphinx.basis,
            self.input.sphinx.PAWHamiltonian,
            self.input.sphinx.initialGuess,
            self.input.sphinx.main
        ]

        if self._generic_input["fix_spin_constraint"]:
            self.input.sphinx.spinConstraint = Group()
            all_groups.append(self.input.sphinx.spinConstraint)
            self.input_writer.write_spin_constraints(
                cwd=self.working_directory
                )
            self.input.sphinx.spinConstraint.setdefault("file", '"spins.in"')

        # In case the entire group was
        # set/overwritten as a normal dict.
        for group in all_groups:
            group = Group(group)

        # write input.sx
        file_name = posixpath.join(self.working_directory, "input.sx")
        with open(file_name, "w") as f:
            f.write(f"//{self.job_name}\n")
            f.write("//SPHInX input file generated by pyiron\n\n")
            f.write("format paw;\n")
            f.write("include <parameters.sx>;\n\n")
            f.write(self.input.sphinx.to_sphinx(indent=0))

    @property
    def _spin_enabled(self):
        if np.any(self.structure.get_initial_magnetic_moments().flatten() != None):
            return True
        return False


    def get_charge_density(self):
        """
        Gets the charge density from the hdf5 file. This value is normalized by the volume

        Returns:
            pyiron.atomistics.volumetric.generic.VolumetricData
        """
        if self.status not in [
            "finished", "warning", "not_converged"
        ]:
            return
        else:
            with self.project_hdf5.open("output") as ho:
                cd_obj = SphinxVolumetricData()
                cd_obj.from_hdf(ho, "charge_density")
            cd_obj.atoms = self.get_structure(-1)
            return cd_obj

    def get_electrostatic_potential(self):
        """
        Gets the electrostatic potential from the hdf5 file.

        Returns:
            pyiron.atomistics.volumetric.generic.VolumetricData
        """
        if self.status not in [
            "finished", "warning", "not_converged"
        ]:
            return
        else:
            with self.project_hdf5.open("output") as ho:
                es_obj = SphinxVolumetricData()
                es_obj.from_hdf(ho, "electrostatic_potential")
            es_obj.atoms = self.get_structure(-1)
            return es_obj

    def collect_output(self, force_update=False):
        """
        Collects the outputs and stores them to the hdf file
        """
        self._output_parser.collect(directory=self.working_directory)
        self._output_parser.to_hdf(self._hdf5, force_update=force_update)

    def convergence_check(self):
        """
        Checks if job has converged according to given cutoffs.
        """
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
        self._logger.info("collect_warnings() is not yet \
            implemented for Sphinx")

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
            is_shift=np.dot(self.structure.cell,
                np.array(self.input["KpointCoords"])),
            is_time_reversal=is_time_reversal,
            symprec=symprec,
        )
        return len(np.unique(mapping))

    def check_setup(self):

        with warnings.catch_warnings(record=True) as w:

            # Check for parameters that were not modified but
            # possibly should have (encut, kpoints, smearing, etc.),
            # or were set to nonsensical values.

            if (
                not (
                    isinstance(self.input.sphinx.basis["eCut"], int)
                    or isinstance(self.input.sphinx.basis["eCut"], float)
                )
                or round(self.input.sphinx.basis["eCut"]*RYDBERG_TO_EV, 0) == 340
            ):
                warnings.warn(
                    "Energy cut-off value wrong or not modified from default "+
                    "340 eV; change it via job.set_encut()"
                )
            if not (
                isinstance(self.input.sphinx.basis["kPoint"]["coords"], np.ndarray)
                or len(self.input.sphinx.basis["kPoint"]["coords"]) != 3
            ):
                warnings.warn("K point coordinates seem to be inappropriate")
            if (
                not (
                    isinstance(self.input.sphinx.PAWHamiltonian["ekt"], int)
                    or isinstance(self.input.sphinx.PAWHamiltonian["ekt"], float)
                )
                or round(self.input.sphinx.PAWHamiltonian["ekt"]*HARTREE_TO_EV, 1) == 0.2
            ):
                warnings.warn(
                    "Fermi smearing value wrong or not modified from default "+
                    "0.2 eV; change it via job.set_occupancy_smearing()"
                )
            if not (
                isinstance(self.input.sphinx.basis["folding"], np.ndarray)
                or len(self.input.sphinx.basis["folding"]) != 3
            ) or self.input.sphinx.basis["folding"].tolist() == [4,4,4]:
                warnings.warn(
                    "K point folding wrong or not modified from default "+
                    "[4,4,4]; change it via job.set_kpoints()"
                )
            if self.get_n_ir_reciprocal_points() < self.server.cores:
                warnings.warn(
                    "Number of cores exceed number of irreducible \
                        reciprocal points: "
                    + str(self.get_n_ir_reciprocal_points())
                )
            if self.input["EmptyStates"] == "auto":
                if any(self.structure.get_initial_magnetic_moments() != None):
                    warnings.warn(
                        "Number of empty states was not specified. Default: "
                        + "3+NIONS*1.5 for magnetic systems. "
                    )
                else:
                    warnings.warn(
                        "Number of empty states was not specified. Default: "
                        + "3+NIONS for non-magnetic systems"
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
        Checks whether parameters are set appropriately. It does not
        mean the simulation won't run if it returns False.
        """

        all_groups = {
            "job.input.pawPot": self.input.sphinx.pawPot,
            "job.input.structure": self.input.sphinx.structure,
            "job.input.basis": self.input.sphinx.basis,
            "job.input.PAWHamiltonian": self.input.sphinx.PAWHamiltonian,
            "job.input.initialGuess": self.input.sphinx.initialGuess,
            "job.input.main": self.input.sphinx.main
        }

        if np.any([len(all_groups[group]) == 0 for group in all_groups]):
            self.load_default_groups()

        if self.structure is None:
            raise AssertionError(
                "Structure not set; set it via job.structure = "
                + "Project().create_structure()"
            )
        if self.input["THREADS"] > self.server.cores:
            raise AssertionError(
                "Number of cores cannot be smaller than the number "
                + "of OpenMP threads"
            )
        with warnings.catch_warnings(record=True) as w:
            # Warn about discrepancies between values in
            # self.input and individual groups, in case
            # a user modified them directly
            if round(self.input["EnCut"], 0)\
                != round(self.input.sphinx.basis.eCut * RYDBERG_TO_EV, 0):
                warnings.warn("job.input.basis.eCut was modified directly. "
                + "It is recommended to set it via job.set_encut()")

            if round(self.input["Sigma"], 1)\
                != round(self.input.sphinx.PAWHamiltonian.ekt * HARTREE_TO_EV, 1):
                warnings.warn("job.input.PAWHamiltonian.ekt was modified directly. "
                + "It is recommended to set it via "
                + "job.set_occupancy_smearing()")

            if self.input["Xcorr"] != self.input.sphinx.PAWHamiltonian.xc:
                warnings.warn("job.input.PAWHamiltonian.xc was modified directly. "
                + "It is recommended to set it via "
                + "job.exchange_correlation_functional()")

            if self.input["EmptyStates"] != self.input.sphinx.PAWHamiltonian.nEmptyStates:
                warnings.warn("job.input.PAWHamiltonian.nEmptyStates was modified "
                + "directly. It is recommended to set it via "
                + "job.set_empty_states()")

            if (
                "KpointCoords" in self.input
                and np.array(self.input.KpointCoords).tolist()
                    != np.array(self.input.sphinx.basis.kPoint.coords).tolist()
                ) \
            or (
                "KpointFolding" in self.input
                and np.array(self.input.KpointFolding).tolist()
                    != np.array(self.input.sphinx.basis.folding).tolist()
                ):

                warnings.warn("job.input.basis.kPoint was modified directly. "
                + "It is recommended to set all k-point settings via "
                + "job.set_kpoints()")

            structure_sync = (str(self.input.sphinx.structure)
                              == str(self.get_structure_group()))
            if not structure_sync and not self.input.sphinx.structure.locked:
                warnings.warn(
                    "job.input.structure != job.structure. "
                    + "The current job.structure will overwrite "
                    + "any changes you may might have made to "
                    + "job.input.structure in the meantime. "
                    + "To disable this overwrite, "
                    + "set job.input.structure.locked = True. "
                    + "To disable this warning, call "
                    + "job.load_structure_group() after making changes "
                    + "to job.structure."
                    )

            if len(w) > 0:
                print("WARNING:")
                for ww in w:
                    print(ww.message)
                return False
            else:
                return True


    def compress(self, files_to_compress=None):
        """
        Compress the output files of a job object.

        Args:
            files_to_compress (list): A list of files to compress (optional)
        """
        # delete empty files
        if files_to_compress is None:
            files_to_compress = [
                f for f in list(self.list_files())
                if (f not in ["rho.sxb", "waves.sxb"]
                and not stat.S_ISFIFO(os.stat(os.path.join(
                    self.working_directory, f
                )).st_mode))
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
        return any(
            [os.path.exists(
                os.path.join(p, 'vasp', 'potentials', 'potpaw', 'Fe', 'POTCAR')
                )
            for p in s.resource_paths]
        )

class InputWriter(object):
    """
    The Sphinx Input writer is called to write the
    Sphinx specific input files.
    """

    def __init__(self):
        self.structure = None
        self._id_pyi_to_spx = []
        self._id_spx_to_pyi = []
        self.file_dict = {}

    def copy_potentials(self, potformat="JTH", xc=None, cwd=None,
                        pot_path_dict=None):

        if pot_path_dict is None:
            pot_path_dict = {}

        if potformat == 'JTH':
            potentials = SphinxJTHPotentialFile(xc=xc)
            find_potential_file = find_potential_file_jth
            pot_path_dict.setdefault("PBE", "jth-gga-pbe")
        elif potformat == 'VASP':
            potentials = VaspPotentialFile(xc=xc)
            find_potential_file = find_potential_file_vasp
            pot_path_dict.setdefault("PBE", "paw-gga-pbe")
            pot_path_dict.setdefault("LDA", "paw-lda")
        else:
            raise ValueError('Only JTH and VASP potentials are supported!')

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
                    path=potentials.find_default(new_element)[
                        "Filename"].values[0][0]
                )
                assert os.path.isfile(
                    potential_path
                ), "such a file does not exist in the pp directory"
            else:
                potential_path = find_potential_file(
                    path=potentials.find_default(elem)[
                        "Filename"].values[0][0]
                )
            if potformat == "JTH":
                copyfile(potential_path, posixpath.join(
                    cwd, elem + "_GGA.atomicdata"
                ))
            else:
                copyfile(potential_path, posixpath.join(
                    cwd, elem + "_POTCAR"
                ))

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
                    self.structure.get_chemical_symbols()
                    == elm_species.Abbreviation
                ]
            )
        self._id_pyi_to_spx = np.array(
            [ooo for oo in self._id_pyi_to_spx for ooo in oo]
        )
        self._id_spx_to_pyi = np.array([0] * len(self._id_pyi_to_spx))
        for i, p in enumerate(self._id_pyi_to_spx):
            self._id_spx_to_pyi[p] = i


    def write_spin_constraints(self, file_name="spins.in", cwd=None,
                               spins_list=None):
        """
        Write a text file containing a list of all spins named spins.in -
        which is used for the external control scripts.

        Args:
            file_name (str): name of the file to be written (optional)
            cwd (str): the current working directory (optinal)
            spins_list (list): the input to write, if no input is
                given the default input will be written. (optional)
        """
        s.logger.debug(f"Writing {file_name}")
        if spins_list is None or len(spins_list) == 0:
            spins_list = []
            s.logger.debug("Getting magnetic moments via \
                get_initial_magnetic_moments")
            if any(self.structure.get_initial_magnetic_moments().flatten()
                != None):
                if any([
                    True
                    if isinstance(spin, list) or isinstance(spin, np.ndarray)
                    else False
                    for spin in self.structure.get_initial_magnetic_moments()
                ]):
                    raise ValueError(
                        "Sphinx only supports collinear spins at the moment."
                    )
                else:
                    for spin, value in zip(
                        self.structure.spin_constraint[self.id_pyi_to_spx],
                        self.structure.get_initial_magnetic_moments()[
                            self.id_pyi_to_spx
                        ],
                    ):
                        if spin:
                            spins_list.append(str(value))
                        else:
                            spins_list.append("X")
                    spins_str = "\n".join(spins_list)
        if spins_str is not None:
            if cwd is not None:
                file_name = posixpath.join(cwd, file_name)
            with open(file_name, "w") as f:
                f.write(spins_str)
        else:
            s.logger.debug("No magnetic moments")

class Group(InputList):
    """
    Dictionary-like object to store SPHInX inputs.

    Attributes (sub-groups, parameters, & flags) can be set
    and accessed via dot notation, or as standard dictionary
    key/values.

    `to_{job_type}` converts the Group to the format
    expected by the given DFT code in its input files.
    """

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        object.__setattr__(self,'locked', False)

    def set(self, name, content):
        self[name] = content

    def set_group(self, name, content=None):
        if content is None:
            self.create_group(name)
        else:
            self.set(name, content)

    def set_flag(self, flag, val=True):
        self.set(flag, val)

    def set_parameter(self, parameter, val):
        self.set(parameter, val)

    def remove(self, name):
        if name in self.keys():
            del self[name]

    def to_sphinx(self, content="__self__", indent=0):

        if content == "__self__":
            content = self

        def format_value(v):
            if isinstance(v, bool):
                if v:
                    return ";"
                else:
                    return " = false;"
            elif isinstance(v, Group):
                if len(v) == 0:
                    return " {}"
                else:
                    return (
                        " {\n"
                        + self.to_sphinx(v, indent+1)
                        + indent * "\t"
                        + "}"
                    )
            else:
                if isinstance(v, np.ndarray):
                    v = v.tolist()
                return " = {!s};".format(v)

        line = ""
        for k, v in content.items():
            if isinstance(v, Group) and len(v) > 0 and not v.has_keys():
                for vv in v.values():
                    line += indent * "\t" + str(k) + format_value(vv) + "\n"
            else:
                line += indent * "\t" + str(k) + format_value(v) + "\n"

        return line

class Output(object):
    """
    Handles the output from a Sphinx simulation.
    """

    def __init__(self, job):
        self._job = job
        self._parse_dict = defaultdict(list)
        self.charge_density = SphinxVolumetricData()
        self.electrostatic_potential = SphinxVolumetricData()

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
            np.array([ss[self._job.id_spx_to_pyi] for ss in spins[:, 1:]]),
            spins[:, 0]
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
                self._parse_dict["bands_eigen_values"] = \
                    np.loadtxt(file_name)[:, 1:]
            except:
                self._parse_dict["bands_eigen_values"] = \
                    np.loadtxt(file_name)[1:]
        else:
            if os.path.isfile(posixpath.join(
                cwd, "eps.0.dat")) and os.path.isfile(
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

    def collect_energy_struct(self, file_name="energy-structOpt.dat",
                              cwd=None):
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
                raise AssertionError("SPHInX did not enter the main loop; \
                    output not collected")
            if not np.any(["Program exited normally." in line
                           for line in log_file]):
                self._job.status.aborted = True
                warnings.warn("SPHInX parsing failed; most likely \
                    SPHInX crashed.")
            main_start = np.where([
                "Enter Main Loop" in line
                for line in log_file]
            )[0][0]
            log_main = log_file[main_start:]

            self._parse_dict["n_valence"] = {
                log_file[ii-1].split()[1]:int(ll.split('=')[-1])
                for ii, ll in enumerate(log_file)
                if ll.startswith('| Z=')
            }

            def get_partial_log(file_content, start_line, end_line):
                start_line = np.where([
                    line == start_line
                    for line in file_content]
                )[0][0]
                end_line = np.where(
                    [line == end_line for line in file_content[start_line:]]
                )[0][0]
                return file_content[start_line : start_line + end_line]
            k_points = get_partial_log(
                log_file,
                "| Symmetrized k-points:               "
                + "in cartesian coordinates\n",
                "\n",
            )[2:-1]
            self._parse_dict["bands_k_weights"] = np.array(
                [float(kk.split()[6]) for kk in k_points]
            )
            k_points = (
                np.array(
                    [[float(kk.split()[i]) for i in range(2, 5)]
                     for kk in k_points]
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
                if line.startswith("eTot(") and not line.startswith(
                    "eTot(Val)")
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
                check_conv(line) for line in log_main
                if check_conv(line) is not None
            ]
            self._parse_dict["bands_e_fermi"] = np.array(
                [
                    float(line.split()[3])
                    for line in log_main
                    if line.startswith("| Fermi energy:")
                ]
            )
            line_vol = np.where(["Omega:" in line for line in log_file])[0][0]
            volume = float(log_file[line_vol].split()[2]) \
                * BOHR_TO_ANGSTROM ** 3
            self._parse_dict["bands_occ"] = [
                line.split()[3:]
                for line in log_main
                if line.startswith("| final focc:")
            ]
            self._parse_dict["bands_occ_initial"] = [
                line.split()[3:] for line in log_main
                if line.startswith("| focc:")
            ]
            self._parse_dict["bands_eigen_values"] = [
                line.split()[4:]
                for line in log_main
                if line.startswith("| final eig [eV]:")
            ]
            self._parse_dict["bands_eigen_values_initial"] = [
                line.split()[4:] for line in log_main
                if line.startswith("| eig [eV]:")
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
                    return np.array(
                        [float(ff) for f in arr for ff in f]
                    ).reshape(
                        -1, 2, len_k_points, len(arr[0])
                    )
                else:
                    return np.array(
                        [float(ff) for f in arr for ff in f]
                    ).reshape(
                        -1, len_k_points, len(arr[0])
                    )

            self._parse_dict["bands_occ"] = eig_converter(
                self._parse_dict["bands_occ"])
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
                forces = np.array(forces).reshape(
                    -1, len(self._job.structure), 3)
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
        if len(self._parse_dict["scf_energy_int"]) == 0 and \
            len(energy_int_lst) != 0:
            self._parse_dict["scf_energy_int"] = energy_int_lst
        if len(self._parse_dict["scf_energy_free"]) == 0 and \
            len(energy_free_lst) != 0:
            self._parse_dict["scf_energy_free"] = energy_free_lst
        if len(self._parse_dict["forces"]) == 0 and len(forces) != 0:
            self._parse_dict["forces"] = forces
        if len(self._parse_dict["scf_magnetic_forces"]) == 0 and \
            len(magnetic_forces) != 0:
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
                [ff[self._job.id_spx_to_pyi]
                for ff in self._parse_dict["positions"]]
            )
            force = np.array(
                [
                    json.loads(line.split("=")[1].split(";")[0])
                    for line in file_content
                    if "force" in line
                ]
            )
            self._parse_dict["forces"] = (
                force.reshape(-1, natoms, 3) *
                HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM
            )
            self._parse_dict["forces"] = np.array(
                [ff[self._job.id_spx_to_pyi]
                for ff in self._parse_dict["forces"]]
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

    def collect_charge_density(self, file_name, cwd):
        if (
            file_name in os.listdir(cwd)
            and os.stat(posixpath.join(cwd, file_name)).st_size != 0
        ):
            self.charge_density.from_file(
                filename=posixpath.join(cwd, file_name),
                normalize=True
            )

    def collect_electrostatic_potential(self, file_name, cwd):
        if (
            file_name in os.listdir(cwd) and
            os.stat(posixpath.join(cwd, file_name)).st_size != 0
        ):
            self.electrostatic_potential.from_file(
                filename=posixpath.join(cwd, file_name),
                normalize=False
            )

    def collect(self, directory=os.getcwd()):
        """
        The collect function, collects all the output from a Sphinx simulation.

        Args:
            directory (str): the directory to collect the output from.
        """
        self.collect_sphinx_log(file_name="sphinx.log", cwd=directory)
        self.collect_energy_dat(file_name="energy.dat", cwd=directory)
        self.collect_residue_dat(file_name="residue.dat", cwd=directory)
        self.collect_eps_dat(file_name="eps.dat", cwd=directory)
        self.collect_spins_dat(file_name="spins.dat", cwd=directory)
        self.collect_energy_struct(file_name="energy-structOpt.dat",
                                   cwd=directory)
        self.collect_relaxed_hist(file_name="relaxHist.sx", cwd=directory)
        self.collect_electrostatic_potential(file_name="vElStat-eV.sxb",
                                             cwd=directory)
        self.collect_charge_density(file_name="rho.sxb",
                                    cwd=directory)
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
                        and "atom_spin_constraints" not in
                        hdf5_dft.list_nodes()
                    ):
                        hdf5_dft["atom_spin_constraints"] = [
                            self._parse_dict["atom_spin_constrains"]
                        ]

        with hdf.open("output") as hdf5_output:
            if "sphinx" not in hdf5_output.list_groups():
                hdf5_output.create_group("sphinx")
            with hdf5_output.open("sphinx") as hdf5_sphinx:
                hdf5_sphinx["bands_occ_initial"] = \
                    self._parse_dict["bands_occ_initial"]
                hdf5_sphinx["bands_eigen_values_initial"] = self._parse_dict[
                    "bands_eigen_values_initial"
                ]
            if self.electrostatic_potential.total_data is not None:
                self.electrostatic_potential.to_hdf(
                    hdf5_output, group_name="electrostatic_potential"
                )
            if self.charge_density.total_data is not None:
                self.charge_density.to_hdf(
                    hdf5_output, group_name="charge_density"
                )
            with hdf5_output.open("generic") as hdf5_generic:
                if "dft" not in hdf5_generic.list_groups():
                    hdf5_generic.create_group("dft")
                with hdf5_generic.open("dft") as hdf5_dft:
                    hdf5_dft["scf_convergence"] = \
                        self._parse_dict["scf_convergence"]
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
                        "n_valence",
                    ]:
                        if len(self._parse_dict[k]) > 0:
                            hdf5_dft[k] = self._parse_dict[k]
                            if "scf_" in k:
                                hdf5_dft[k.replace("scf_", "")] = np.array(
                                    [vv[-1] for vv in self._parse_dict[k]]
                                )
                if len(self._parse_dict["scf_computation_time"]) > 0:
                    hdf5_generic["computation_time"] = np.array(
                        [tt[-1] for tt in
                         self._parse_dict["scf_computation_time"]]
                    )
                if len([en[-1] for en in
                    self._parse_dict["scf_energy_free"]]) > 0:
                    hdf5_generic["energy_tot"] = np.array(
                        [en[-1] for en in self._parse_dict["scf_energy_free"]]
                    )
                    hdf5_generic["energy_pot"] = np.array(
                        [en[-1] for en in self._parse_dict["scf_energy_free"]]
                    )
                hdf5_generic["volume"] = self._parse_dict["volume"]
                if "positions" not in hdf5_generic.list_nodes() or \
                    force_update:
                    if len(self._parse_dict["positions"]) > 0:
                        hdf5_generic["positions"] = np.array(
                            self._parse_dict["positions"]
                        )
                    elif len(self._parse_dict["scf_convergence"]) == 1:
                        hdf5_generic["positions"] = np.array(
                            [self._job.structure.positions]
                        )
                if ("forces" not in hdf5_generic.list_nodes() or force_update)\
                    and len(
                    self._parse_dict["forces"]
                ) > 0:
                    hdf5_generic["forces"] = \
                        np.array(self._parse_dict["forces"])
                if "cells" not in hdf5_generic.list_nodes() or force_update:
                    if len(self._parse_dict["cell"]) > 0:
                        hdf5_generic["cells"] = np.array(
                            self._parse_dict["cell"])
                    elif len(self._parse_dict["scf_convergence"]) == 1:
                        hdf5_generic["cells"] = np.array(
                            [self._job.structure.cell])


    def from_hdf(self, hdf):
        """
        Load output from an HDF5 file
        """
