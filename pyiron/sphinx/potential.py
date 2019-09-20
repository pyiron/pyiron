# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import pandas
from pyiron.base.settings.generic import Settings
from pyiron.vasp.potential import VaspPotentialAbstract

__author__ = "Osamu Waseda"
__copyright__ = (
    "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Osamu Waseda"
__email__ = "waseda@mpie.de"
__status__ = "development"
__date__ = "Sep 20, 2019"

s = Settings()


class SphinxPotentialFile(VaspPotentialAbstract):
    """
    The Potential class is derived from the PotentialAbstract class, but instead of loading the potentials from a list,
    the potentials are loaded from a file.

    Args:
        xc (str): Exchange correlation functional ['PBE', 'LDA']
    """
    def __init__(self, xc=None, selected_atoms=None):
        potential_df = self._get_potential_df(
            plugin_name="sphinx",
            file_name_lst={"potentials_sphinx.csv"},
            backward_compatibility_name="sphinxpotentials",
        )
        if xc == "PBE":
            default_df = self._get_potential_default_df(
                plugin_name="vasp",
                file_name_lst={"potentials_vasp_pbe_default.csv"},
                backward_compatibility_name="defaultvasppbe",
            )
            potential_df = potential_df[(potential_df["Model"] == "gga-pbe")]
        elif xc == "GGA":
            default_df = self._get_potential_default_df(
                plugin_name="vasp",
                file_name_lst={"potentials_vasp_pbe_default.csv"},
                backward_compatibility_name="defaultvasppbe",
            )
            potential_df = potential_df[(potential_df["Model"] == "gga-pbe")]
        elif xc == "LDA":
            default_df = self._get_potential_default_df(
                plugin_name="vasp",
                file_name_lst={"potentials_vasp_lda_default.csv"},
                backward_compatibility_name="defaultvasplda",
            )
            potential_df = potential_df[(potential_df["Model"] == "lda")]
        elif xc == "JTH":
            default_df = self._get_potential_default_df(
                plugin_name="sphinx",
                file_name_lst={"potentials_sphinx_jth_default.csv"},
                backward_compatibility_name="defaultsphinxjth",
            )
            potential_df = potential_df[(potential_df["Model"] == "jth-gga-pbe")]
        else:
            raise ValueError(
                'The exchange correlation functional has to be set and it can either be "LDA" or "PBE"'
            )
        super(SphinxPotentialFile, self).__init__(
            potential_df=potential_df,
            default_df=default_df,
            selected_atoms=selected_atoms,
        )

    def add_new_element(self, parent_element, new_element):
        """
        Adding a new user defined element with a different POTCAR file. It is assumed that the file exists

        Args:
            parent_element (str): Parent element
            new_element (str): Name of the new element (the name of the folder where the new POTCAR file exists

        """
        ds = self.find_default(element=parent_element)
        ds["Species"].values[0][0] = new_element
        path_list = ds["Filename"].values[0][0].split("/")
        path_list[-2] = new_element
        name_list = ds["Name"].values[0].split("-")
        name_list[0] = new_element
        ds["Name"].values[0] = "-".join(name_list)
        ds["Filename"].values[0][0] = "/".join(path_list)
        self._potential_df = self._potential_df.append(ds)
        ds = pandas.Series()
        ds.name = new_element
        ds["Name"] = "-".join(name_list)
        self._default_df = self._default_df.append(ds)


def find_potential_file(file_name=None, xc=None, path=None, pot_path_dict=None):
    if path is not None:
        for resource_path in s.resource_paths:
            if os.path.exists(os.path.join(resource_path, "sphinx", "potentials", path)):
                return os.path.join(resource_path, "sphinx", "potentials", path)
    elif xc is not None and file_name is not None:
        for resource_path in s.resource_paths:
            if os.path.exists(
                os.path.join(resource_path, "sphinx", "potentials", pot_path_dict[xc])
            ):
                resource_path = os.path.join(
                    resource_path, "sphinx", "potentials", pot_path_dict[xc]
                )
            if "potentials" in resource_path:
                for path, folder_lst, file_lst in os.walk(resource_path):
                    if file_name in file_lst:
                        return os.path.join(path, file_name)
    raise ValueError("Either the filename or the functional has to be defined.")

