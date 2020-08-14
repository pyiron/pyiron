# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, unicode_literals
import os
import posixpath

import sys
import h5py
import numpy as np
import pandas as pd
import warnings
from io import StringIO

from pyiron.lammps.potential import LammpsPotentialFile, PotentialAvailable
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.base.settings.generic import Settings
from pyiron.base.pyio.parser import Logstatus, extract_data_from_file
from pyiron.lammps.control import LammpsControl
from pyiron.lammps.potential import LammpsPotential
from pyiron.lammps.structure import LammpsStructure, UnfoldingPrism

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH "
    "- Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


s = Settings()


class LammpsBase(AtomisticGenericJob):
    """
    Class to setup and run and analyze LAMMPS simulations which is a derivative of
    atomistics.job.generic.GenericJob. The functions in these modules are written in such the function names and
    attributes are very generic (get_structure(), molecular_dynamics(), version) but the functions are written to handle
    LAMMPS specific input/output.

    Args:
        project (pyiron.project.Project instance):  Specifies the project path among other attributes
        job_name (str): Name of the job

    Attributes:
        input (lammps.Input instance): Instance which handles the input
    """

    def __init__(self, project, job_name):
        super(LammpsBase, self).__init__(project, job_name)
        self.input = Input()
        self._cutoff_radius = None
        self._is_continuation = None
        self._compress_by_default = True
        self._prism = None
        s.publication_add(self.publication)

    @property
    def bond_dict(self):
        """
        A dictionary which defines the nature of LAMMPS bonds that are to be drawn between atoms. To set the values, use
        the function `define_bonds`.

        Returns:
            dict: Dictionary of the bond properties for every species

        """
        return self.input.bond_dict

    def define_bonds(self, species, element_list, cutoff_list, max_bond_list, bond_type_list, angle_type_list=None):
        """
        Define the nature of bonds between different species. Make sure that the bonds between two species are defined
        only once (no double counting).

        Args:
            species (str): Species for which the bonds are to be drawn (e.g. O, H, C ..)
            element_list (list): List of species to which the bonds are to be made (e.g. O, H, C, ..)
            cutoff_list (list): Draw bonds only for atoms within this cutoff distance
            max_bond_list (list): Maximum number of bonds drawn from each molecule
            bond_type_list (list): Type of the bond as defined in the LAMMPS potential file
            angle_type_list (list): Type of the angle as defined in the LAMMPS potential file

        Example:
            The command below defined bonds between O and H atoms within a cutoff raduis of 2 $\AA$ with the bond and
            angle types 1 defined in the potential file used

            >> job_lammps.define_bonds(species="O", element_list-["H"], cutoff_list=[2.0], bond_type_list=[1],
            angle_type_list=[1])

        """
        if isinstance(species, str):
            if len(element_list) == len(cutoff_list) == bond_type_list == max_bond_list:
                self.input.bond_dict[species] = dict()
                self.input.bond_dict[species]["element_list"] = element_list
                self.input.bond_dict[species]["cutoff_list"] = cutoff_list
                self.input.bond_dict[species]["bond_type_list"] = bond_type_list
                self.input.bond_dict[species]["max_bond_list"] = max_bond_list
                if angle_type_list is not None:
                    self.input.bond_dict[species]["angle_type_list"] = angle_type_list
                else:
                    self.input.bond_dict[species]["angle_type_list"] = [None]
            else:
                raise ValueError("The element list, cutoff list, max bond list, and the bond type list"
                                 " must have the same length")

    @property
    def cutoff_radius(self):
        """

        Returns:

        """
        return self._cutoff_radius

    @cutoff_radius.setter
    def cutoff_radius(self, cutoff):
        """

        Args:
            cutoff:

        Returns:

        """
        self._cutoff_radius = cutoff

    @property
    def potential(self):
        """
        Execute view_potential() or list_potential() in order to see the pre-defined potential files

        Returns:

        """
        return self.input.potential.df

    @potential.setter
    def potential(self, potential_filename):
        """
        Execute view_potential() or list_potential() in order to see the pre-defined potential files

        Args:
            potential_filename:

        Returns:

        """
        stringtypes = str
        if isinstance(potential_filename, stringtypes):
            if ".lmp" in potential_filename:
                potential_filename = potential_filename.split(".lmp")[0]
            potential_db = LammpsPotentialFile()
            potential = potential_db.find_by_name(potential_filename)
        elif isinstance(potential_filename, pd.DataFrame):
            potential = potential_filename
        else:
            raise TypeError("Potentials have to be strings or pandas dataframes.")
        if self.structure:
            structure_elements = self.structure.get_species_symbols()
            potential_elements = list(potential["Species"])[0]
            if not set(structure_elements).issubset(potential_elements):
                raise ValueError("Potential {} does not support elements "
                                 "in structure {}.".format(
                                     potential_elements,
                                     structure_elements
                                ))
        self.input.potential.df = potential
        for val in ["units", "atom_style", "dimension"]:
            v = self.input.potential[val]
            if v is not None:
                self.input.control[val] = v
                if val == "units" and v != "metal":
                    warnings.warn(
                        "WARNING: Non-'metal' units are not fully supported. Your calculation should run OK, but "
                        "results may not be saved in pyiron units."
                    )
        self.input.potential.remove_structure_block()

    @property
    def potential_available(self):
        return PotentialAvailable(list_of_potentials=self.potential_list)

    @property
    def potential_list(self):
        """
        List of interatomic potentials suitable for the current atomic structure.

        use self.potentials_view() to get more details.

        Returns:
            list: potential names
        """

        return self.list_potentials()

    @property
    def potential_view(self):
        """
        List all interatomic potentials for the current atomistic sturcture including all potential parameters.

        To quickly get only the names of the potentials you can use: self.potentials_list()

        Returns:
            pandas.Dataframe: Dataframe including all potential parameters.
        """
        return self.view_potentials()

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(LammpsBase, self).set_input_to_read_only()
        self.input.control.read_only = True
        self.input.potential.read_only = True

    def validate_ready_to_run(self):
        """
        Validating input parameters before LAMMPS run
        """
        super(LammpsBase, self).validate_ready_to_run()
        if self.potential is None:
            raise ValueError(
                "This job does not contain a valid potential: {}".format(self.job_name)
            )
        scaled_positions = self.structure.get_scaled_positions(wrap=False)
        # Check if atoms located outside of non periodic box
        conditions = [(np.min(scaled_positions[:, i]) < 0.0 or
                       np.max(scaled_positions[:, i]) > 1.0) and not self.structure.pbc[i] for i in range(3)]
        if any(conditions):
            raise ValueError("You have atoms located outside the non-periodic boundaries "
                             "of the defined simulation box")

    def get_potentials_for_structure(self):
        """

        Returns:

        """
        return self.list_potentials()

    def get_final_structure(self):
        """

        Returns:

        """
        warnings.warn(
            "get_final_structure() is deprecated - please use get_structure() instead.",
            DeprecationWarning,
        )
        return self.get_structure(iteration_step=-1)

    def view_potentials(self):
        """
        List all interatomic potentials for the current atomistic sturcture including all potential parameters.

        To quickly get only the names of the potentials you can use: self.potentials_list()

        Returns:
            pandas.Dataframe: Dataframe including all potential parameters.
        """
        from pyiron.lammps.potential import LammpsPotentialFile

        if not self.structure:
            raise ValueError("No structure set.")
        list_of_elements = set(self.structure.get_chemical_symbols())
        list_of_potentials = LammpsPotentialFile().find(list_of_elements)
        if list_of_potentials is not None:
            return list_of_potentials
        else:
            raise TypeError(
                "No potentials found for this kind of structure: ",
                str(list_of_elements),
            )

    def list_potentials(self):
        """
        List of interatomic potentials suitable for the current atomic structure.

        use self.potentials_view() to get more details.

        Returns:
            list: potential names
        """
        return list(self.view_potentials()["Name"].values)

    def enable_h5md(self):
        """

        Returns:

        """
        del self.input.control["dump_modify___1"]
        del self.input.control["dump___1"]
        self.input.control[
            "dump___1"
        ] = "all h5md ${dumptime} dump.h5 position force create_group yes"

    def write_input(self):
        """
        Call routines that generate the code specific input files

        Returns:

        """
        if self.structure is None:
            raise ValueError("Input structure not set. Use method set_structure()")
        lmp_structure = self._get_lammps_structure(
            structure=self.structure, cutoff_radius=self.cutoff_radius
        )
        lmp_structure.write_file(file_name="structure.inp", cwd=self.working_directory)
        version_int_lst = self._get_executable_version_number()
        if (
            version_int_lst is not None
            and "dump_modify" in self.input.control._dataset["Parameter"]
            and (
                version_int_lst[0] < 2016
                or (version_int_lst[0] == 2016 and version_int_lst[1] < 11)
            )
        ):
            self.input.control["dump_modify"] = self.input.control[
                "dump_modify"
            ].replace(" line ", " ")
        if not all(self.structure.pbc):
            self.input.control["boundary"] = " ".join(
                ["p" if coord else "f" for coord in self.structure.pbc]
            )
        self._set_selective_dynamics()
        self.input.control.write_file(
            file_name="control.inp", cwd=self.working_directory
        )
        self.input.potential.write_file(
            file_name="potential.inp", cwd=self.working_directory
        )
        self.input.potential.copy_pot_files(self.working_directory)

    def _get_executable_version_number(self):
        """
        Get the version of the executable

        Returns:
            list: List of integers defining the version number
        """
        if self.executable.version:
            return [
                l
                for l in [
                    [int(i) for i in sv.split(".") if i.isdigit()]
                    for sv in self.executable.version.split("/")[-1].split("_")
                ]
                if len(l) > 0
            ][0]
        else:
            return None

    @property
    def publication(self):
        return {
            "lammps": {
                "lammps": {
                    "title": "Fast Parallel Algorithms for Short-Range Molecular Dynamics",
                    "journal": "Journal of Computational Physics",
                    "volume": "117",
                    "number": "1",
                    "pages": "1-19",
                    "year": "1995",
                    "issn": "0021-9991",
                    "doi": "10.1006/jcph.1995.1039",
                    "url": "http://www.sciencedirect.com/science/article/pii/S002199918571039X",
                    "author": ["Steve Plimpton"],
                }
            }
        }

    def collect_output(self):
        """

        Returns:

        """
        self.input.from_hdf(self._hdf5)
        if os.path.isfile(
            self.job_file_name(file_name="dump.h5", cwd=self.working_directory)
        ):
            self.collect_h5md_file(file_name="dump.h5", cwd=self.working_directory)
        else:
            self.collect_dump_file(file_name="dump.out", cwd=self.working_directory)
        self.collect_output_log(file_name="log.lammps", cwd=self.working_directory)
        final_structure = self.get_structure(iteration_step=-1)
        with self.project_hdf5.open("output") as hdf_output:
            final_structure.to_hdf(hdf_output)

    def convergence_check(self):
        if self._generic_input["calc_mode"] == "minimize":
            if (
                self._generic_input["max_iter"] + 1
                <= len(self["output/generic/energy_tot"])
                or len(
                    [l for l in self["log.lammps"] if "linesearch alpha is zero" in l]
                )
                != 0
            ):
                return False
            else:
                return True
        else:
            return True

    def collect_logfiles(self):
        """

        Returns:

        """
        return

    # TODO: make rotation of all vectors back to the original as in self.collect_dump_file
    def collect_h5md_file(self, file_name="dump.h5", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        prism = UnfoldingPrism(self.structure.cell, digits=15)
        if np.matrix.trace(prism.R) != 3:
            raise RuntimeError("The Lammps output will not be mapped back to pyiron correctly.")
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
        with h5py.File(file_name, mode="r", libver="latest", swmr=True) as h5md:
            positions = [
                pos_i.tolist() for pos_i in h5md["/particles/all/position/value"]
            ]
            steps = [steps_i.tolist() for steps_i in h5md["/particles/all/position/step"]]
            forces = [for_i.tolist() for for_i in h5md["/particles/all/force/value"]]
            # following the explanation at: http://nongnu.org/h5md/h5md.html
            cell = [
                np.eye(3) * np.array(cell_i.tolist())
                for cell_i in h5md["/particles/all/box/edges/value"]
            ]
            indices = [indices_i.tolist() for indices_i in h5md["/particles/all/indices/value"]]
        with self.project_hdf5.open("output/generic") as h5_file:
            h5_file["forces"] = np.array(forces)
            h5_file["positions"] = np.array(positions)
            h5_file["steps"] = np.array(steps)
            h5_file["cells"] = cell
            h5_file["indices"] = self.remap_indices(indices)

    def remap_indices(self, lammps_indices):
        """
        Give the Lammps-dumped indices, re-maps these back onto the structure's indices to preserve the species.

        The issue is that for an N-element potential, Lammps dumps the chemical index from 1 to N based on the order
        that these species are written in the Lammps input file. But the indices for a given structure are based on the
        order in which chemical species were added to that structure, and run from 0 up to the number of species
        currently in that structure. Therefore we need to be a little careful with mapping.

        Args:
            indices (numpy.ndarray/list): The Lammps-dumped integers.

        Returns:
            numpy.ndarray: Those integers mapped onto the structure.
        """
        lammps_symbol_order = np.array(self.input.potential.get_element_lst())

        # If new Lammps indices are present for which we have no species, extend the species list
        unique_lammps_indices = np.unique(lammps_indices)
        if len(unique_lammps_indices) > len(np.unique(self.structure.indices)):
            unique_lammps_indices -= 1  # Convert from Lammps start counting at 1 to python start counting at 0
            new_lammps_symbols = lammps_symbol_order[unique_lammps_indices]
            self.structure.set_species([self.structure.convert_element(el) for el in new_lammps_symbols])

        # Create a map between the lammps indices and structure indices to preserve species
        structure_symbol_order = np.array([el.Abbreviation for el in self.structure.species])
        map_ = np.array([int(np.argwhere(lammps_symbol_order == symbol)[0]) + 1 for symbol in structure_symbol_order])

        structure_indices = np.array(lammps_indices)
        for i_struct, i_lammps in enumerate(map_):
            np.place(structure_indices, lammps_indices == i_lammps, i_struct)
        # TODO: Vectorize this for-loop for computational efficiency

        return structure_indices

    def collect_errors(self, file_name, cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
        error = extract_data_from_file(file_name, tag="ERROR", num_args=1000)
        if len(error) > 0:
            error = " ".join(error[0])
            raise RuntimeError("Run time error occurred: " + str(error))
        else:
            return True

    def collect_output_log(self, file_name="log.lammps", cwd=None):
        """
        general purpose routine to extract static from a lammps log file

        Args:
            file_name:
            cwd:

        Returns:

        """
        self.collect_errors(file_name=file_name, cwd=cwd)
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
        with open(file_name, "r") as f:
            f = f.readlines()
            l_start = np.where([line.startswith("Step") for line in f])[0]
            l_end = np.where([line.startswith("Loop") for line in f])[0]
            if len(l_start) > len(l_end):
                l_end = np.append(l_end, [None])
            df = [
                pd.read_csv(
                    StringIO("\n".join(f[llst:llen])), delim_whitespace=True
                )
                for llst, llen in zip(l_start, l_end)
            ]
        df = df[-1]

        h5_dict = {
            "Step": "steps",
            "Temp": "temperature",
            "PotEng": "energy_pot",
            "TotEng": "energy_tot",
            "Volume": "volume",
        }

        for key in df.columns[df.columns.str.startswith('f_mean')]:
            h5_dict[key] = key.replace('f_', '')

        df = df.rename(index=str, columns=h5_dict)
        pressures = np.stack(
            (df.Pxx, df.Pxy, df.Pxz, df.Pxy, df.Pyy, df.Pyz, df.Pxz, df.Pyz, df.Pzz),
            axis=-1,
        ).reshape(-1, 3, 3).astype('float64')
        pressures *= 0.0001  # bar -> GPa

        # Rotate pressures from Lammps frame to pyiron frame if necessary
        rotation_matrix = self._prism.R.T
        if np.matrix.trace(rotation_matrix) != 3:
            pressures = rotation_matrix.T @ pressures @ rotation_matrix

        df = df.drop(
            columns=df.columns[
                ((df.columns.str.len() == 3) & df.columns.str.startswith("P"))
            ]
        )
        df["pressures"] = pressures.tolist()
        if 'mean_pressure[1]' in df.columns:
            pressures = np.stack(
                (df['mean_pressure[1]'], df['mean_pressure[4]'], df['mean_pressure[5]'],
                 df['mean_pressure[4]'], df['mean_pressure[2]'], df['mean_pressure[6]'],
                 df['mean_pressure[5]'], df['mean_pressure[6]'], df['mean_pressure[3]']),
                axis=-1,
            ).reshape(-1, 3, 3).astype('float64')
            pressures *= 0.0001  # bar -> GPa
            if np.matrix.trace(rotation_matrix) != 3:
                pressures = rotation_matrix.T @ pressures @ rotation_matrix
            df = df.drop(
                columns=df.columns[
                    (df.columns.str.startswith("mean_pressure") & df.columns.str.endswith(']'))
                ]
            )
            df["mean_pressures"] = pressures.tolist()

        with self.project_hdf5.open("output/generic") as hdf_output:
            # This is a hack for backward comparability
            for k, v in df.items():
                hdf_output[k] = np.array(v)

    def calc_minimize(
            self,
            ionic_energy_tolerance=0.0,
            ionic_force_tolerance=1e-4,
            e_tol=None,
            f_tol=None,
            max_iter=1000000,
            pressure=None,
            n_print=100,
            style='cg'
    ):
        rotation_matrix = self._get_rotation_matrix(pressure=pressure)
        # Docstring set programmatically -- Ensure that changes to signature or defaults stay consistent!
        if e_tol is not None:
            ionic_energy_tolerance = e_tol
        if f_tol is not None:
            ionic_force_tolerance = f_tol
        super(LammpsBase, self).calc_minimize(
            ionic_energy_tolerance=ionic_energy_tolerance,
            ionic_force_tolerance=ionic_force_tolerance,
            e_tol=e_tol,
            f_tol=f_tol,
            max_iter=max_iter,
            pressure=pressure,
            n_print=n_print,
        )
        self.input.control.calc_minimize(
            ionic_energy_tolerance=ionic_energy_tolerance,
            ionic_force_tolerance=ionic_force_tolerance,
            max_iter=max_iter,
            pressure=pressure,
            n_print=n_print,
            style=style,
            rotation_matrix=rotation_matrix
        )
    calc_minimize.__doc__ = LammpsControl.calc_minimize.__doc__

    def calc_static(self):
        """

        Returns:

        """
        super(LammpsBase, self).calc_static()
        self.input.control.calc_static()

    def calc_md(
        self,
        temperature=None,
        pressure=None,
        n_ionic_steps=1000,
        time_step=1.0,
        n_print=100,
        temperature_damping_timescale=100.0,
        pressure_damping_timescale=1000.0,
        seed=None,
        tloop=None,
        initial_temperature=None,
        langevin=False,
        delta_temp=None,
        delta_press=None,
    ):
        # Docstring set programmatically -- Ensure that changes to signature or defaults stay consistent!
        if self.server.run_mode.interactive_non_modal:
            warnings.warn(
                "calc_md() is not implemented for the non modal interactive mode use calc_static()!"
            )
        rotation_matrix = self._get_rotation_matrix(pressure=pressure)
        super(LammpsBase, self).calc_md(
            temperature=temperature,
            pressure=pressure,
            n_ionic_steps=n_ionic_steps,
            time_step=time_step,
            n_print=n_print,
            temperature_damping_timescale=temperature_damping_timescale,
            pressure_damping_timescale=pressure_damping_timescale,
            seed=seed,
            tloop=tloop,
            initial_temperature=initial_temperature,
            langevin=langevin,
        )
        self.input.control.calc_md(
            temperature=temperature,
            pressure=pressure,
            n_ionic_steps=n_ionic_steps,
            time_step=time_step,
            n_print=n_print,
            temperature_damping_timescale=temperature_damping_timescale,
            pressure_damping_timescale=pressure_damping_timescale,
            seed=seed,
            tloop=tloop,
            initial_temperature=initial_temperature,
            langevin=langevin,
            delta_temp=delta_temp,
            delta_press=delta_press,
            job_name=self.job_name,
            rotation_matrix=rotation_matrix
        )
    calc_md.__doc__ = LammpsControl.calc_md.__doc__

    def calc_vcsgc(
        self,
        mu=None,
        target_concentration=None,
        kappa=1000.,
        mc_step_interval=100,
        swap_fraction=0.1,
        temperature_mc=None,
        window_size=None,
        window_moves=None,
        temperature=None,
        pressure=None,
        n_ionic_steps=1000,
        time_step=1.0,
        n_print=100,
        temperature_damping_timescale=100.0,
        pressure_damping_timescale=1000.0,
        seed=None,
        initial_temperature=None,
        langevin=False
    ):
        """
        Run variance-constrained semi-grand-canonical MD/MC for a binary system. In addition to VC-SGC arguments, all
        arguments for a regular MD calculation are also accepted.

        https://vcsgc-lammps.materialsmodeling.org

        Note:
            For easy visualization later (with `get_structure`), it is highly recommended that the initial structure
            contain at least one atom of each species.

        Warning:
            - The fix does not yet support non-orthogonal simulation boxes; using one will give a runtime error.

        Args:
            mu (dict): A dictionary of chemical potentials, one for each element the potential treats, where the
                dictionary keys are just the chemical symbol. Note that only the *relative* chemical potentials are used
                here, such that the swap acceptance probability is influenced by the chemical potential difference
                between the two species (a more negative value increases the odds of swapping *to* that element.)
                (Default is None, all elements have the same chemical potential.)
            target_concentration: A dictionary of target simulation domain concentrations for each species *in the
                potential*. Dictionary keys should be the chemical symbol of the corresponding species, and the sum of
                all concentrations must be 1. (Default is None, which runs regular semi-grand-canonical MD/MC without
                any variance constraint.)
            kappa: Variance constraint for the MC. Larger value means a tighter adherence to the target concentrations.
                (Default is 1000.)
            mc_step_interval (int): How many steps of MD between each set of MC moves. (Default is 100.) Must divide the
                number of ionic steps evenly.
            swap_fraction (float): The fraction of atoms whose species is swapped at each MC phase. (Default is 0.1.)
            temperature_mc (float): The temperature for accepting MC steps. (Default is None, which uses the MD
                temperature.)
            window_size (float): The size of the sampling window for parallel calculations as a fraction of something
                unspecified in the VC-SGC docs, but it must lie between 0.5 and 1. (Default is None, window is
                determined automatically.)
            window_moves (int): The number of times the sampling window is moved during one MC cycle. (Default is None,
                number of moves is determined automatically.)
        """
        if mu is None:
            mu = {}
            for el in self.input.potential.get_element_lst():
                mu[el] = 0.

        self._generic_input["calc_mode"] = "vcsgc"
        self._generic_input["temperature"] = temperature
        self._generic_input["n_ionic_steps"] = n_ionic_steps
        self._generic_input["n_print"] = n_print
        self._generic_input.remove_keys(["max_iter"])
        self.input.control.calc_vcsgc(
            mu=mu,
            ordered_element_list=self.input.potential.get_element_lst(),
            target_concentration=target_concentration,
            kappa=kappa,
            mc_step_interval=mc_step_interval,
            swap_fraction=swap_fraction,
            temperature_mc=temperature_mc,
            window_size=window_size,
            window_moves=window_moves,
            temperature=temperature,
            pressure=pressure,
            n_ionic_steps=n_ionic_steps,
            time_step=time_step,
            n_print=n_print,
            temperature_damping_timescale=temperature_damping_timescale,
            pressure_damping_timescale=pressure_damping_timescale,
            seed=seed,
            initial_temperature=initial_temperature,
            langevin=langevin,
            job_name=self.job_name,
        )

    # define hdf5 input and output
    def to_hdf(self, hdf=None, group_name=None):
        """

        Args:
            hdf:
            group_name:

        Returns:

        """
        super(LammpsBase, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()
        self.input.to_hdf(self._hdf5)

    def from_hdf(self, hdf=None, group_name=None):  # TODO: group_name should be removed
        """

        Args:
            hdf:
            group_name:

        Returns:

        """
        super(LammpsBase, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()
        self.input.from_hdf(self._hdf5)

    def write_restart_file(self, filename="restart.out"):
        """

        Args:
            filename:

        Returns:

        """
        self.input.control.modify(write_restart=filename, append_if_not_present=True)

    def compress(self, files_to_compress=None):
        """
        Compress the output files of a job object.

        Args:
            files_to_compress (list):
        """
        if files_to_compress is None:
            files_to_compress = [
                f for f in list(self.list_files()) if f not in ["restart.out"]
            ]
        super(LammpsBase, self).compress(files_to_compress=files_to_compress)

    def read_restart_file(self, filename="restart.out"):
        """

        Args:
            filename:

        Returns:

        """
        self._is_continuation = True
        self.input.control.set(read_restart=filename)
        self.input.control["reset_timestep"] = 0
        self.input.control.remove_keys(
            ["dimension", "read_data", "boundary", "atom_style", "velocity"]
        )

    def collect_dump_file(self, file_name="dump.out", cwd=None):
        """
        general purpose routine to extract static from a lammps dump file

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = self.job_file_name(file_name=file_name, cwd=cwd)
        output = {}
        with open(file_name, "r") as ff:
            dump = ff.readlines()

        steps = np.genfromtxt(
            [
                dump[nn]
                for nn in np.where([ll.startswith("ITEM: TIMESTEP") for ll in dump])[0]
                + 1
            ],
            dtype=int,
        )
        steps = np.array([steps]).flatten()
        output["steps"] = steps

        natoms = np.genfromtxt(
            [
                dump[nn]
                for nn in np.where(
                    [ll.startswith("ITEM: NUMBER OF ATOMS") for ll in dump]
                )[0]
                + 1
            ],
            dtype=int,
        )
        natoms = np.array([natoms]).flatten()

        prism = self._prism
        rotation_lammps2orig = self._prism.R.T
        cells = np.genfromtxt(
            " ".join(
                (
                    [
                        " ".join(dump[nn:nn + 3])
                        for nn in np.where(
                            [ll.startswith("ITEM: BOX BOUNDS") for ll in dump]
                        )[0]
                        + 1
                    ]
                )
            ).split()
        ).reshape(len(natoms), -1)
        lammps_cells = np.array([to_amat(cc) for cc in cells])
        unfolded_cells = np.array([prism.unfold_cell(cell) for cell in lammps_cells])
        output["cells"] = unfolded_cells


        l_start = np.where([ll.startswith("ITEM: ATOMS") for ll in dump])[0]
        l_end = l_start + natoms + 1
        content = [
            pd.read_csv(
                StringIO("\n".join(dump[llst:llen]).replace("ITEM: ATOMS ", "")),
                delim_whitespace=True,
            )
            for llst, llen in zip(l_start, l_end)
        ]

        indices = np.array([cc["type"] for cc in content], dtype=int)
        output["indices"] = self.remap_indices(indices)

        forces = np.array(
            [np.stack((cc["fx"], cc["fy"], cc["fz"]), axis=-1) for cc in content]
        )
        output["forces"] = np.matmul(forces, rotation_lammps2orig)

        if 'f_mean_forces[1]' in content[0].keys():
            forces = np.array(
                [np.stack((cc["f_mean_forces[1]"],
                           cc["f_mean_forces[2]"],
                           cc["f_mean_forces[3]"]),
                          axis=-1) for cc in content]
            )
            output["mean_forces"] = np.matmul(forces, rotation_lammps2orig)

        if np.all([flag in content[0].columns.values for flag in ["vx", "vy", "vz"]]):
            velocities = np.array(
                [np.stack((cc["vx"], cc["vy"], cc["vz"]), axis=-1) for cc in content]
            )
            output["velocities"] = np.matmul(velocities, rotation_lammps2orig)

        if 'f_mean_velocities[1]' in content[0].keys():
            velocities = np.array(
                [np.stack((cc["f_mean_velocities[1]"],
                           cc["f_mean_velocities[2]"],
                           cc["f_mean_velocities[3]"]),
                          axis=-1) for cc in content]
            )
            output["mean_velocities"] = np.matmul(velocities, rotation_lammps2orig)
        direct_unwrapped_positions = np.array(
            [np.stack((cc["xsu"], cc["ysu"], cc["zsu"]), axis=-1) for cc in content]
        )
        unwrapped_positions = np.matmul(direct_unwrapped_positions, lammps_cells)
        output["unwrapped_positions"] = np.matmul(unwrapped_positions, rotation_lammps2orig)
        if 'f_mean_positions[1]' in content[0].keys():
            direct_unwrapped_positions = np.array(
                [np.stack((cc["f_mean_positions[1]"],
                           cc["f_mean_positions[2]"],
                           cc["f_mean_positions[3]"]),
                          axis=-1) for cc in content]
            )
            unwrapped_positions = np.matmul(direct_unwrapped_positions, lammps_cells)
            output["mean_unwrapped_positions"] = np.matmul(unwrapped_positions, rotation_lammps2orig)

        direct_positions = direct_unwrapped_positions - np.floor(direct_unwrapped_positions)
        positions = np.matmul(direct_positions, lammps_cells)
        output["positions"] = np.matmul(positions, rotation_lammps2orig)

        keys = content[0].keys()
        for kk in keys[keys.str.startswith('c_')]:
            output[kk.replace('c_', '')] = np.array([cc[kk] for cc in content], dtype=float)

        with self.project_hdf5.open("output/generic") as hdf_output:
            for k, v in output.items():
                hdf_output[k] = v

    # Outdated functions:
    def set_potential(self, file_name):
        """

        Args:
            file_name:

        Returns:

        """
        print("This function is outdated use the potential setter instead!")
        self.potential = file_name

    def next(self, job_name=None, job_type=None):
        """
        Restart a new job created from an existing Lammps calculation.
        Args:
            project (pyiron.project.Project instance): Project instance at which the new job should be created
            job_name (str): Job name
            job_type (str): Job type. If not specified a Lammps job type is assumed

        Returns:
            new_ham (lammps.lammps.Lammps instance): New job
        """
        return super(LammpsBase, self).restart(
            job_name=job_name, job_type=job_type
        )

    def restart(self, job_name=None, job_type=None):
        """
        Restart a new job created from an existing Lammps calculation.
        Args:
            project (pyiron.project.Project instance): Project instance at which the new job should be created
            job_name (str): Job name
            job_type (str): Job type. If not specified a Lammps job type is assumed

        Returns:
            new_ham (lammps.lammps.Lammps instance): New job
        """
        new_ham = super(LammpsBase, self).restart(
            job_name=job_name, job_type=job_type
        )
        if new_ham.__name__ == self.__name__:
            new_ham.potential = self.potential
            if os.path.isfile(os.path.join(self.working_directory, "restart.out")):
                new_ham.read_restart_file(filename="restart.out")
                new_ham.restart_file_list.append(
                    posixpath.join(self.working_directory, "restart.out")
                )
        return new_ham

    def _get_lammps_structure(self, structure=None, cutoff_radius=None):
        lmp_structure = LammpsStructure(bond_dict=self.input.bond_dict)
        lmp_structure.potential = self.input.potential
        lmp_structure.atom_type = self.input.control["atom_style"]
        if cutoff_radius is not None:
            lmp_structure.cutoff_radius = cutoff_radius
        else:
            lmp_structure.cutoff_radius = self.cutoff_radius
        lmp_structure.el_eam_lst = self.input.potential.get_element_lst()

        def structure_to_lammps(structure):
            """
            Converts a structure to the Lammps coordinate frame

            Args:
                structure (pyiron.atomistics.structure.atoms.Atoms): Structure to convert.

            Returns:
                pyiron.atomistics.structure.atoms.Atoms: Structure with the LAMMPS coordinate frame.
            """
            prism = UnfoldingPrism(structure.cell)
            lammps_structure = structure.copy()
            lammps_structure.set_cell(prism.A)
            lammps_structure.positions = np.matmul(structure.positions, prism.R)
            return lammps_structure

        if structure is not None:
            lmp_structure.structure = structure_to_lammps(structure)
        else:
            lmp_structure.structure = structure_to_lammps(self.structure)
        if not set(lmp_structure.structure.get_species_symbols()).issubset(
            set(lmp_structure.el_eam_lst)
        ):
            raise ValueError(
                "The selected potentials do not support the given combination of elements."
            )
        return lmp_structure

    def _set_selective_dynamics(self):
        if "selective_dynamics" in self.structure._tag_list.keys():
            if self.structure.selective_dynamics._default is None:
                self.structure.selective_dynamics._default = [True, True, True]
            sel_dyn = np.logical_not(self.structure.selective_dynamics.list())
            # Enter loop only if constraints present
            if len(np.argwhere(np.any(sel_dyn, axis=1)).flatten()) != 0:
                all_indices = np.arange(len(self.structure), dtype=int)
                constraint_xyz = np.argwhere(np.all(sel_dyn, axis=1)).flatten()
                not_constrained_xyz = np.setdiff1d(all_indices, constraint_xyz)
                # LAMMPS starts counting from 1
                constraint_xyz += 1
                ind_x = np.argwhere(sel_dyn[not_constrained_xyz, 0]).flatten()
                ind_y = np.argwhere(sel_dyn[not_constrained_xyz, 1]).flatten()
                ind_z = np.argwhere(sel_dyn[not_constrained_xyz, 2]).flatten()
                constraint_xy = not_constrained_xyz[np.intersect1d(ind_x, ind_y)] + 1
                constraint_yz = not_constrained_xyz[np.intersect1d(ind_y, ind_z)] + 1
                constraint_zx = not_constrained_xyz[np.intersect1d(ind_z, ind_x)] + 1
                constraint_x = (
                    not_constrained_xyz[np.setdiff1d(np.setdiff1d(ind_x, ind_y), ind_z)]
                    + 1
                )
                constraint_y = (
                    not_constrained_xyz[np.setdiff1d(np.setdiff1d(ind_y, ind_z), ind_x)]
                    + 1
                )
                constraint_z = (
                    not_constrained_xyz[np.setdiff1d(np.setdiff1d(ind_z, ind_x), ind_y)]
                    + 1
                )
                if len(constraint_xyz) > 0:
                    self.input.control["group___constraintxyz"] = "id " + " ".join(
                        [str(ind) for ind in constraint_xyz]
                    )
                    self.input.control[
                        "fix___constraintxyz"
                    ] = "constraintxyz setforce 0.0 0.0 0.0"
                    if self._generic_input["calc_mode"] == "md":
                        self.input.control[
                            "velocity___constraintxyz"
                        ] = "set 0.0 0.0 0.0"
                if len(constraint_xy) > 0:
                    self.input.control["group___constraintxy"] = "id " + " ".join(
                        [str(ind) for ind in constraint_xy]
                    )
                    self.input.control[
                        "fix___constraintxy"
                    ] = "constraintxy setforce 0.0 0.0 NULL"
                    if self._generic_input["calc_mode"] == "md":
                        self.input.control[
                            "velocity___constraintxy"
                        ] = "set 0.0 0.0 NULL"
                if len(constraint_yz) > 0:
                    self.input.control["group___constraintyz"] = "id " + " ".join(
                        [str(ind) for ind in constraint_yz]
                    )
                    self.input.control[
                        "fix___constraintyz"
                    ] = "constraintyz setforce NULL 0.0 0.0"
                    if self._generic_input["calc_mode"] == "md":
                        self.input.control[
                            "velocity___constraintyz"
                        ] = "set NULL 0.0 0.0"
                if len(constraint_zx) > 0:
                    self.input.control["group___constraintxz"] = "id " + " ".join(
                        [str(ind) for ind in constraint_zx]
                    )
                    self.input.control[
                        "fix___constraintxz"
                    ] = "constraintxz setforce 0.0 NULL 0.0"
                    if self._generic_input["calc_mode"] == "md":
                        self.input.control[
                            "velocity___constraintxz"
                        ] = "set 0.0 NULL 0.0"
                if len(constraint_x) > 0:
                    self.input.control["group___constraintx"] = "id " + " ".join(
                        [str(ind) for ind in constraint_x]
                    )
                    self.input.control[
                        "fix___constraintx"
                    ] = "constraintx setforce 0.0 NULL NULL"
                    if self._generic_input["calc_mode"] == "md":
                        self.input.control[
                            "velocity___constraintx"
                        ] = "set 0.0 NULL NULL"
                if len(constraint_y) > 0:
                    self.input.control["group___constrainty"] = "id " + " ".join(
                        [str(ind) for ind in constraint_y]
                    )
                    self.input.control[
                        "fix___constrainty"
                    ] = "constrainty setforce NULL 0.0 NULL"
                    if self._generic_input["calc_mode"] == "md":
                        self.input.control[
                            "velocity___constrainty"
                        ] = "set NULL 0.0 NULL"
                if len(constraint_z) > 0:
                    self.input.control["group___constraintz"] = "id " + " ".join(
                        [str(ind) for ind in constraint_z]
                    )
                    self.input.control[
                        "fix___constraintz"
                    ] = "constraintz setforce NULL NULL 0.0"
                    if self._generic_input["calc_mode"] == "md":
                        self.input.control[
                            "velocity___constraintz"
                        ] = "set NULL NULL 0.0"

    @staticmethod
    def _modify_structure_to_allow_requested_deformation(structure, pressure, prism=None):
        """
        Lammps will not allow xy/xz/yz cell deformations in minimization or MD for non-triclinic cells. In case the
        requested pressure for a calculation has these non-diagonal entries, we need to make sure it will run. One way
        to do this is by invoking the lammps `change_box` command, but it is easier to just force our box to to be
        triclinic by adding a very small cell perturbation (in the case where it isn't triclinic already).

        Args:
            pressure (float/int/list/numpy.ndarray/tuple): Between three and six pressures for the x, y, z, xy, xz, and
                yz directions, in that order, or a single value.
        """
        if hasattr(pressure, '__len__'):
            non_diagonal_pressures = np.any([p is not None for p in pressure[3:]])

            if prism is None:
                prism = UnfoldingPrism(structure.cell)

            if non_diagonal_pressures:
                try:
                    if not prism.is_skewed():
                        skew_structure = structure.copy()
                        skew_structure.cell[0, 1] += 2 * prism.acc
                        return skew_structure
                except AttributeError:
                    warnings.warn(
                        "WARNING: Setting a calculation type which uses pressure before setting the structure risks " +
                        "constraining your cell shape evolution if non-diagonal pressures are used but the structure " +
                        "is not triclinic from the start of the calculation."
                    )
        return structure

    def _get_rotation_matrix(self, pressure):
        """

        Args:
            pressure:

        Returns:

        """
        if self.structure is not None:
            if self._prism is None:
                self._prism = UnfoldingPrism(self.structure.cell)

            self.structure = self._modify_structure_to_allow_requested_deformation(
                pressure=pressure,
                structure=self.structure,
                prism=self._prism
            )
            rotation_matrix = self._prism.R
        else:
            warnings.warn(
                "No structure set, can not validate the simulation cell!"
            )
            rotation_matrix = None
        return rotation_matrix


class Input:
    def __init__(self):
        self.control = LammpsControl()
        self.potential = LammpsPotential()
        self.bond_dict = dict()
        # Set default bond parameters
        self._load_default_bond_params()

    def _load_default_bond_params(self):
        """
        Function to automatically load a few default bond params (wont automatically write them)

        """
        # Default bond properties of a water molecule
        self.bond_dict["O"] = dict()
        self.bond_dict["O"]["element_list"] = ["H"]
        self.bond_dict["O"]["cutoff_list"] = [2.0]
        self.bond_dict["O"]["max_bond_list"] = [2]
        self.bond_dict["O"]["bond_type_list"] = [1]
        self.bond_dict["O"]["angle_type_list"] = [1]

    def to_hdf(self, hdf5):
        """

        Args:
            hdf5:

        Returns:

        """
        with hdf5.open("input") as hdf5_input:
            self.control.to_hdf(hdf5_input)
            self.potential.to_hdf(hdf5_input)

    def from_hdf(self, hdf5):
        """

        Args:
            hdf5:

        Returns:

        """
        with hdf5.open("input") as hdf5_input:
            self.control.from_hdf(hdf5_input)
            self.potential.from_hdf(hdf5_input)
            if "bond_dict" in hdf5_input.list_nodes():
                self.bond_dict = hdf5_input["bond_dict"]


def to_amat(l_list):
    """

    Args:
        l_list:

    Returns:

    """
    lst = np.reshape(l_list, -1)
    if len(lst) == 9:
        xlo_bound, xhi_bound, xy, ylo_bound, yhi_bound, xz, zlo_bound, zhi_bound, yz = (
            lst
        )

    elif len(lst) == 6:
        xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound = lst
        xy, xz, yz = 0.0, 0.0, 0.0
    else:
        raise ValueError("This format for amat not yet implemented: " + str(len(lst)))

    # > xhi_bound - xlo_bound = xhi -xlo  + MAX(0.0, xy, xz, xy + xz) - MIN(0.0, xy, xz, xy + xz)
    # > xhili = xhi -xlo   = xhi_bound - xlo_bound - MAX(0.0, xy, xz, xy + xz) + MIN(0.0, xy, xz, xy + xz)
    xhilo = (
        (xhi_bound - xlo_bound)
        - max([0.0, xy, xz, xy + xz])
        + min([0.0, xy, xz, xy + xz])
    )

    # > yhilo = yhi -ylo = yhi_bound -ylo_bound - MAX(0.0, yz) + MIN(0.0, yz)
    yhilo = (yhi_bound - ylo_bound) - max([0.0, yz]) + min([0.0, yz])

    # > zhi - zlo = zhi_bound- zlo_bound
    zhilo = zhi_bound - zlo_bound

    cell = [[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]]
    return cell
