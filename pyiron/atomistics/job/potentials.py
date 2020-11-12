# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
An abstract Potential class to provide an easy access for the available potentials. Currently implemented for the
OpenKim https://openkim.org database.
"""
import pandas
import os
from pyiron_base import Settings

__author__ = "Martin Boeckmann, Jan Janssen"
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


class PotentialAbstract(object):
    """
    The PotentialAbstract class loads a list of available potentials and sorts them. Afterwards the potentials can be
    accessed through:
        PotentialAbstract.<Element>.<Element> or PotentialAbstract.find_potentials_set({<Element>, <Element>}

    Args:
        potential_df:
        default_df:
        selected_atoms:
    """

    def __init__(self, potential_df, default_df=None, selected_atoms=None):
        self._potential_df = potential_df
        self._default_df = default_df
        if selected_atoms is not None:
            self._selected_atoms = selected_atoms
        else:
            self._selected_atoms = []

    def find(self, element):
        """
        Find the potentials

        Args:
            element (set, str): element or set of elements for which you want the possible LAMMPS potentials

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
        return self._potential_df[
            [
                True if set(element).issubset(species) else False
                for species in self._potential_df["Species"].values
            ]
        ]

    def find_by_name(self, potential_name):
        return self._potential_df[(self._potential_df["Name"] == potential_name)]

    def list(self):
        """
        List the available potentials

        Returns:
            list: of possible potentials for the element or the combination of elements
        """
        return self._potential_df

    def __getattr__(self, item):
        return self[item]

    def __getitem__(self, item):
        potential_df = self.find(element=item)
        selected_atoms = self._selected_atoms + [item]
        return PotentialAbstract(
            potential_df=potential_df,
            default_df=self._default_df,
            selected_atoms=selected_atoms,
        )

    def __str__(self):
        return str(self.list())

    @staticmethod
    def _get_potential_df(plugin_name, file_name_lst, backward_compatibility_name):
        """

        Args:
            plugin_name (str):
            file_name_lst (set):
            backward_compatibility_name (str):

        Returns:
            pandas.DataFrame:
        """
        env = os.environ
        resource_path_lst = s.resource_paths
        for conda_var in ["CONDA_PREFIX", "CONDA_DIR"]:
            if conda_var in env.keys():  # support iprpy-data package
                resource_path_lst += [os.path.join(env[conda_var], "share", "iprpy")]
        for resource_path in resource_path_lst:
            if os.path.exists(os.path.join(resource_path, plugin_name, "potentials")):
                resource_path = os.path.join(resource_path, plugin_name, "potentials")
            if "potentials" in resource_path:
                for path, folder_lst, file_lst in os.walk(resource_path):
                    for periodic_table_file_name in file_name_lst:
                        if (
                            periodic_table_file_name in file_lst
                            and periodic_table_file_name.endswith(".csv")
                        ):
                            return pandas.read_csv(
                                os.path.join(path, periodic_table_file_name),
                                index_col=0,
                                converters={
                                    "Species": lambda x: x.replace("'", "")
                                    .strip("[]")
                                    .split(", "),
                                    "Config": lambda x: x.replace("'", "")
                                    .replace("\\n", "\n")
                                    .strip("[]")
                                    .split(", "),
                                    "Filename": lambda x: x.replace("'", "")
                                    .strip("[]")
                                    .split(", "),
                                },
                            )
                        elif (
                            periodic_table_file_name in file_lst
                            and periodic_table_file_name.endswith(".h5")
                        ):
                            return pandas.read_hdf(
                                os.path.join(path, periodic_table_file_name), mode="r"
                            )
        raise ValueError("Was not able to locate the potential files.")

    @staticmethod
    def _get_potential_default_df(
        plugin_name,
        file_name_lst={"potentials_vasp_pbe_default.csv"},
        backward_compatibility_name="defaultvasppbe",
    ):
        """

        Args:
            plugin_name (str):
            file_name_lst (set):
            backward_compatibility_name (str):

        Returns:
            pandas.DataFrame:
        """
        for resource_path in s.resource_paths:
            pot_path = os.path.join(resource_path, plugin_name, "potentials")
            if os.path.exists(pot_path):
                resource_path = pot_path
            if "potentials" in resource_path:
                for path, folder_lst, file_lst in os.walk(resource_path):
                    for periodic_table_file_name in file_name_lst:
                        if (
                            periodic_table_file_name in file_lst
                            and periodic_table_file_name.endswith(".csv")
                        ):
                            return pandas.read_csv(
                                os.path.join(path, periodic_table_file_name),
                                index_col=0,
                            )
                        elif (
                            periodic_table_file_name in file_lst
                            and periodic_table_file_name.endswith(".h5")
                        ):
                            return pandas.read_hdf(
                                os.path.join(path, periodic_table_file_name), mode="r"
                            )
        raise ValueError("Was not able to locate the potential files.")


def find_potential_file_base(path, resource_path_lst, rel_path):
    if path is not None:
        for resource_path in resource_path_lst:
            path_direct = os.path.join(resource_path, path)
            path_indirect = os.path.join(resource_path, rel_path, path)
            if os.path.exists(path_direct):
                return path_direct
            elif os.path.exists(path_indirect):
                return path_indirect
    raise ValueError("Either the filename or the functional has to be defined.",
                     path, resource_path_lst)
