# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from ast import literal_eval
import numpy as np
import pandas as pd
import shutil
import os
from pyiron_base import Settings, GenericParameters
from pyiron.atomistics.job.potentials import PotentialAbstract, find_potential_file_base

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


class LammpsPotential(GenericParameters):

    """
    This module helps write commands which help in the control of parameters related to the potential used in LAMMPS
    simulations
    """

    def __init__(self, input_file_name=None):
        super(LammpsPotential, self).__init__(
            input_file_name=input_file_name,
            table_name="potential_inp",
            comment_char="#",
        )
        self._potential = None
        self._attributes = {}
        self._df = None

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, new_dataframe):
        self._df = new_dataframe
        # ToDo: In future lammps should also support more than one potential file - that is currently not implemented.
        try:
            self.load_string("".join(list(new_dataframe["Config"])[0]))
        except IndexError:
            raise ValueError(
                "Potential not found! "
                "Validate the potential name by self.potential in self.list_potentials()."
            )

    def remove_structure_block(self):
        self.remove_keys(["units"])
        self.remove_keys(["atom_style"])
        self.remove_keys(["dimension"])

    @property
    def files(self):
        if len(self._df["Filename"].values[0]) > 0 and self._df["Filename"].values[0] != ['']:
            absolute_file_paths = [
                files for files in list(self._df["Filename"])[0] if os.path.isabs(files)
            ]
            relative_file_paths = [
                files
                for files in list(self._df["Filename"])[0]
                if not os.path.isabs(files)
            ]
            env = os.environ
            resource_path_lst = s.resource_paths
            for conda_var in ["CONDA_PREFIX", "CONDA_DIR"]:
                if conda_var in env.keys():  # support iprpy-data package
                    resource_path_lst += [os.path.join(env[conda_var], "share", "iprpy")]
            for path in relative_file_paths:
                absolute_file_paths.append(find_potential_file_base(
                    path=path,
                    resource_path_lst=resource_path_lst,
                    rel_path=os.path.join("lammps", "potentials")
                ))
            if len(absolute_file_paths) != len(list(self._df["Filename"])[0]):
                raise ValueError("Was not able to locate the potentials.")
            else:
                return absolute_file_paths

    def copy_pot_files(self, working_directory):
        if self.files is not None:
            _ = [shutil.copy(path_pot, working_directory) for path_pot in self.files]

    def get_element_lst(self):
        return list(self._df["Species"])[0]

    def _find_line_by_prefix(self, prefix):
        """
        Find a line that starts with the given prefix.  Differences in white
        space are ignored.  Raises a ValueError if not line matches the prefix.

        Args:
            prefix (str): line prefix to search for

        Returns:
            list: words of the matching line

        Raises:
            ValueError: if not matching line was found
        """

        def isprefix(prefix, lst):
            if len(prefix) > len(lst): return False
            return all(n == l for n, l in zip(prefix, lst))

        # compare the line word by word to also match lines that differ only in
        # whitespace
        prefix = prefix.split()
        for parameter, value in zip(self._dataset["Parameter"],
                                    self._dataset["Value"]):
            words = (parameter + " " + value).strip().split()
            if isprefix(prefix, words):
                return words

        raise ValueError("No line with prefix \"{}\" found.".format(
                            " ".join(prefix)))

    def get_element_id(self, element_symbol):
        """
        Return numeric element id for element. If potential does not contain
        the element raise a :class:NameError.  Only makes sense for potentials
        with pair_style "full".

        Args:
            element_symbol (str): short symbol for element

        Returns:
            int: id matching the given symbol

        Raise:
            NameError: if potential does not contain this element
        """

        try:
            line = "group {} type".format(element_symbol)
            return int(self._find_line_by_prefix(line)[3])

        except ValueError:
            msg = "potential does not contain element {}".format(
                    element_symbol)
            raise NameError(msg) from None

    def get_charge(self, element_symbol):
        """
        Return charge for element. If potential does not specify a charge,
        raise a :class:NameError.  Only makes sense for potentials
        with pair_style "full".

        Args:
            element_symbol (str): short symbol for element

        Returns:
            float: charge speicified for the given element

        Raises:
            NameError: if potential does not specify charge for this element
        """

        try:
            line = "set group {} charge".format(element_symbol)
            return float(self._find_line_by_prefix(line)[4])

        except ValueError:
            msg = "potential does not specify charge for element {}".format(
                    element_symbol)
            raise NameError(msg) from None

    def to_hdf(self, hdf, group_name=None):
        if self._df is not None:
            with hdf.open("potential") as hdf_pot:
                hdf_pot["Config"] = self._df["Config"].values[0]
                hdf_pot["Filename"] = self._df["Filename"].values[0]
                hdf_pot["Name"] = self._df["Name"].values[0]
                hdf_pot["Model"] = self._df["Model"].values[0]
                hdf_pot["Species"] = self._df["Species"].values[0]
                if "Citations" in self._df.columns.values:
                    hdf_pot["Citations"] = self._df["Citations"].values[0]
        super(LammpsPotential, self).to_hdf(hdf, group_name=group_name)

    def from_hdf(self, hdf, group_name=None):
        with hdf.open("potential") as hdf_pot:
            try:
                entry_dict = {
                    "Config": [hdf_pot["Config"]],
                    "Filename": [hdf_pot["Filename"]],
                    "Name": [hdf_pot["Name"]],
                    "Model": [hdf_pot["Model"]],
                    "Species": [hdf_pot["Species"]],
                }
                if "Citations" in hdf_pot.list_nodes():
                    entry_dict["Citations"] = [hdf_pot["Citations"]]
                self._df = pd.DataFrame(entry_dict)
            except ValueError:
                pass
        super(LammpsPotential, self).from_hdf(hdf, group_name=group_name)


class LammpsPotentialFile(PotentialAbstract):
    """
    The Potential class is derived from the PotentialAbstract class, but instead of loading the potentials from a list,
    the potentials are loaded from a file.

    Args:
        potential_df:
        default_df:
        selected_atoms:
    """

    def __init__(self, potential_df=None, default_df=None, selected_atoms=None):
        if potential_df is None:
            potential_df = self._get_potential_df(
                plugin_name="lammps",
                file_name_lst={"potentials_lammps.csv"},
                backward_compatibility_name="lammpspotentials",
            )
        super(LammpsPotentialFile, self).__init__(
            potential_df=potential_df,
            default_df=default_df,
            selected_atoms=selected_atoms,
        )

    def default(self):
        if self._default_df is not None:
            atoms_str = "_".join(sorted(self._selected_atoms))
            return self._default_df[
                (self._default_df["Name"] == self._default_df.loc[atoms_str].values[0])
            ]
        return None

    def find_default(self, element):
        """
        Find the potentials

        Args:
            element (set, str): element or set of elements for which you want the possible LAMMPS potentials
            path (bool): choose whether to return the full path to the potential or just the potential name

        Returns:
            list: of possible potentials for the element or the combination of elements

        """
        if isinstance(element, set):
            element = element
        elif isinstance(element, list):
            element = set(element)
        elif isinstance(element, str):
            element = set([element])
        else:
            raise TypeError("Only, str, list and set supported!")
        element_lst = list(element)
        if self._default_df is not None:
            merged_lst = list(set(self._selected_atoms + element_lst))
            atoms_str = "_".join(sorted(merged_lst))
            return self._default_df[
                (self._default_df["Name"] == self._default_df.loc[atoms_str].values[0])
            ]
        return None

    def __getitem__(self, item):
        potential_df = self.find(element=item)
        selected_atoms = self._selected_atoms + [item]
        return LammpsPotentialFile(
            potential_df=potential_df,
            default_df=self._default_df,
            selected_atoms=selected_atoms,
        )


class PotentialAvailable(object):
    def __init__(self, list_of_potentials):
        self._list_of_potentials = {
            "pot_" + v.replace("-", "_").replace('.', '_'): v for v in list_of_potentials
        }

    def __getattr__(self, name):
        if name in self._list_of_potentials.keys():
            return self._list_of_potentials[name]
        else:
            raise AttributeError

    def __dir__(self):
        return list(self._list_of_potentials.keys())

    def __repr__(self):
        return str(dir(self))
