# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import posixpath

import numpy as np
import pandas
import tables
from pyiron.base.generic.parameters import GenericParameters
from pyiron.base.settings.generic import Settings
from pyiron.atomistics.job.potentials import PotentialAbstract

__author__ = "Jan Janssen"
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


class VaspPotentialAbstract(PotentialAbstract):
    """

    Args:
        potential_df:
        default_df:
        selected_atoms:
    """

    def __init__(self, potential_df=None, default_df=None, selected_atoms=None):
        if potential_df is None:
            potential_df = self._get_potential_df(
                plugin_name="vasp",
                file_name_lst={"potentials_vasp.csv"},
                backward_compatibility_name="vasppotentials",
            )
        super(VaspPotentialAbstract, self).__init__(
            potential_df=potential_df,
            default_df=default_df,
            selected_atoms=selected_atoms,
        )

    def default(self):
        if self._default_df is not None:
            return pandas.concat(
                [
                    self._potential_df[
                        (
                            self._potential_df["Name"]
                            == self._default_df.loc[atom].values[0]
                        )
                    ]
                    for atom in self._selected_atoms
                ]
            )
        return None

    def find_default(self, element):
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
            return pandas.concat(
                [
                    self._potential_df[
                        (
                            self._potential_df["Name"]
                            == self._default_df.loc[atom].values[0]
                        )
                    ]
                    for atom in merged_lst
                ]
            )
        return None

    def find(self, element):
        if isinstance(element, set):
            element = element
        elif isinstance(element, list):
            element = set(element)
        elif isinstance(element, str):
            element = set([element])
        else:
            raise TypeError("Only, str, list and set supported!")
        element_lst = list(element)
        merged_lst = list(set(self._selected_atoms + element_lst))
        return pandas.concat(
            [super(VaspPotentialAbstract, self).find({atom}) for atom in merged_lst]
        )

    def list(self):
        if len(self._selected_atoms) != 0:
            return pandas.concat(
                [
                    super(VaspPotentialAbstract, self).find({atom})
                    for atom in self._selected_atoms
                ]
            )
        else:
            return pandas.DataFrame({})

    def list_potential_names(self):
        df = self.list()
        if len(df) != 0:
            return list(self.list()["Name"])
        else:
            return []

    @staticmethod
    def _return_potential_file(file_name):
        for resource_path in s.resource_paths:
            resource_path_potcar = os.path.join(
                resource_path, "vasp", "potentials", file_name
            )
            if os.path.exists(resource_path_potcar):
                return resource_path_potcar
        return None

    def __dir__(self):
        return [val.replace("-", "_") for val in self.list_potential_names()]

    def __getitem__(self, item):
        item_replace = item.replace("_gga_pbe", "-gga-pbe").replace("_lda", "-lda")
        if item_replace in self.list_potential_names():
            df = self.list()
            return self._return_potential_file(
                file_name=list(df[df["Name"] == item_replace]["Filename"])[0][0]
            )
        selected_atoms = self._selected_atoms + [item]
        return VaspPotentialAbstract(
            potential_df=self._potential_df,
            default_df=self._default_df,
            selected_atoms=selected_atoms,
        )


class VaspPotentialFile(VaspPotentialAbstract):
    """
    The Potential class is derived from the PotentialAbstract class, but instead of loading the potentials from a list,
    the potentials are loaded from a file.

    Args:
        xc (str): Exchange correlation functional ['PBE', 'LDA']
    """

    def __init__(self, xc=None, selected_atoms=None):
        potential_df = self._get_potential_df(
            plugin_name="vasp",
            file_name_lst={"potentials_vasp.csv"},
            backward_compatibility_name="vasppotentials",
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
        else:
            raise ValueError(
                'The exchange correlation functional has to be set and it can either be "LDA" or "PBE"'
            )
        super(VaspPotentialFile, self).__init__(
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


class VaspPotential(object):
    """
    The Potential class is derived from the PotentialAbstract class, but instead of loading the potentials from a list,
    the potentials are loaded from a file.

    Args:
        path (str): path to the potential list
    """

    def __init__(self, selected_atoms=None):
        self.pbe = VaspPotentialFile(xc="PBE", selected_atoms=selected_atoms)
        self.lda = VaspPotentialFile(xc="LDA", selected_atoms=selected_atoms)


class VaspPotentialSetter(object):
    def __init__(self, element_lst):
        super(VaspPotentialSetter, self).__setattr__("_element_lst", element_lst)
        super(VaspPotentialSetter, self).__setattr__(
            "_potential_dict", {el: None for el in element_lst}
        )

    def __getattr__(self, item):
        if item in self._element_lst:
            return item
        else:
            raise AttributeError

    def __setattr__(self, key, value):
        if key in self._element_lst:
            self._potential_dict[key] = value
        else:
            raise AttributeError

    def to_dict(self):
        return self._potential_dict

    def __repr__(self):
        return self._potential_dict.__repr__()


def find_potential_file(file_name=None, xc=None, path=None, pot_path_dict=None):
    if path is not None:
        for resource_path in s.resource_paths:
            if os.path.exists(os.path.join(resource_path, "vasp", "potentials", path)):
                return os.path.join(resource_path, "vasp", "potentials", path)
    elif xc is not None and file_name is not None:
        for resource_path in s.resource_paths:
            if os.path.exists(
                os.path.join(resource_path, "vasp", "potentials", pot_path_dict[xc])
            ):
                resource_path = os.path.join(
                    resource_path, "vasp", "potentials", pot_path_dict[xc]
                )
            if "potentials" in resource_path:
                for path, folder_lst, file_lst in os.walk(resource_path):
                    if file_name in file_lst:
                        return os.path.join(path, file_name)
    raise ValueError("Either the filename or the functional has to be defined.")


def get_enmax_among_species(symbol_lst, return_list=False, xc="PBE"):
    """
    Given a list of species symbols, finds the largest applicable encut.

    Args:
        symbol_lst (list): The list of species symbols.
        return_list (bool): Whether to return the list of all ENMAX values (in the same order as `species_lst` along with
            the largest value). (Default is False.)
        xc ("GGA"/"PBE"/"LDA"): The exchange correlation functional for which the POTCARs were generated. (Default is "PBE".)

    Returns:
        (float): The largest ENMAX among the POTCAR files for all the species.
        [optional](list): The ENMAX value corresponding to each species.
    """
    pot_path_dict = Potcar.pot_path_dict

    enmax_lst = []
    vpf = VaspPotentialFile(xc=xc)

    for symbol in symbol_lst:
        potcar_file = find_potential_file(
            path=vpf.find_default(symbol)['Filename'].values[0][0],
            pot_path_dict=pot_path_dict
        )
        with open(potcar_file) as pf:
            for i, line in enumerate(pf):
                if i == 14:
                    encut_str = line.split()[2][:-1]
                    enmax_lst.append(float(encut_str))
                    break

    if return_list:
        return max(enmax_lst), enmax_lst
    else:
        return max(enmax_lst)


class Potcar(GenericParameters):
    pot_path_dict = {"GGA": "paw-gga-pbe", "PBE": "paw-gga-pbe", "LDA": "paw-lda"}

    def __init__(self, input_file_name=None, table_name="potcar"):
        GenericParameters.__init__(
            self,
            input_file_name=input_file_name,
            table_name=table_name,
            val_only=False,
            comment_char="#",
        )
        self._structure = None
        self.electrons_per_atom_lst = list()
        self.max_cutoff_lst = list()
        self.el_path_lst = list()
        self.el_path_dict = dict()
        self.modified_elements = dict()

    def potcar_set_structure(self, structure, modified_elements):
        self._structure = structure
        self._set_default_path_dict()
        self._set_potential_paths()
        self.modified_elements = modified_elements

    def modify(self, **modify):
        if "xc" in modify:
            xc_type = modify["xc"]
            self._set_default_path_dict()
            if xc_type not in self.pot_path_dict:
                raise ValueError("xc type not implemented: " + xc_type)
        GenericParameters.modify(self, **modify)
        if self._structure is not None:
            self._set_potential_paths()

    def _set_default_path_dict(self):
        if self._structure is None:
            return
        vasp_potentials = VaspPotentialFile(xc=self.get("xc"))
        for i, el_obj in enumerate(self._structure.get_species_objects()):
            if isinstance(el_obj.Parent, str):
                el = el_obj.Parent
            else:
                el = el_obj.Abbreviation
            if isinstance(el_obj.tags, dict):
                if "pseudo_potcar_file" in el_obj.tags.keys():
                    new_element = el_obj.tags["pseudo_potcar_file"]
                    vasp_potentials.add_new_element(
                        parent_element=el, new_element=new_element
                    )
            key = vasp_potentials.find_default(el).Species.values[0][0]
            val = vasp_potentials.find_default(el).Name.values[0]
            self[key] = val

    def _set_potential_paths(self):
        element_list = (
            self._structure.get_species_symbols()
        )  # .ElementList.getSpecies()
        object_list = self._structure.get_species_objects()
        s.logger.debug("element list: {0}".format(element_list))
        self.el_path_lst = list()
        try:
            xc = self.get("xc")
        except tables.exceptions.NoSuchNodeError:
            xc = self.get("xc")
        s.logger.debug("XC: {0}".format(xc))
        vasp_potentials = VaspPotentialFile(xc=xc)
        for i, el_obj in enumerate(object_list):
            if isinstance(el_obj.Parent, str):
                el = el_obj.Parent
            else:
                el = el_obj.Abbreviation
            if (
                isinstance(el_obj.tags, dict)
                and "pseudo_potcar_file" in el_obj.tags.keys()
            ):
                new_element = el_obj.tags["pseudo_potcar_file"]
                vasp_potentials.add_new_element(
                    parent_element=el, new_element=new_element
                )
                el_path = find_potential_file(
                    path=vasp_potentials.find_default(new_element)["Filename"].values[
                        0
                    ][0],
                    pot_path_dict=self.pot_path_dict,
                )
                if not (os.path.isfile(el_path)):
                    raise ValueError("such a file does not exist in the pp directory")
            else:
                el_path = find_potential_file(
                    path=vasp_potentials.find_default(el)["Filename"].values[0][0],
                    pot_path_dict=self.pot_path_dict,
                )

            if not (os.path.isfile(el_path)):
                raise AssertionError()
            pot_name = "pot_" + str(i)

            if pot_name in self._dataset["Parameter"]:
                try:
                    ind = self._dataset["Parameter"].index(pot_name)
                except (ValueError, IndexError):
                    indices = np.core.defchararray.find(
                        self._dataset["Parameter"], pot_name
                    )
                    ind = np.where(indices == 0)[0][0]
                self._dataset["Value"][ind] = el_path
                self._dataset["Comment"][ind] = ""
            else:
                self._dataset["Parameter"].append("pot_" + str(i))
                self._dataset["Value"].append(el_path)
                self._dataset["Comment"].append("")
            if el_obj.Abbreviation in self.modified_elements.keys():
                self.el_path_lst.append(self.modified_elements[el_obj.Abbreviation])
            else:
                self.el_path_lst.append(el_path)

    def write_file(self, file_name, cwd=None):
        """
        Args:
            file_name:
            cwd:
        Returns:
        """
        self.electrons_per_atom_lst = list()
        self.max_cutoff_lst = list()
        self._set_potential_paths()
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)
        f = open(file_name, "w")
        for el_file in self.el_path_lst:
            with open(el_file) as pot_file:
                for i, line in enumerate(pot_file):
                    f.write(line)
                    if i == 1:
                        self.electrons_per_atom_lst.append(int(float(line)))
                    elif i == 14:
                        mystr = line.split()[2][:-1]
                        self.max_cutoff_lst.append(float(mystr))
        f.close()

    def load_default(self):
        file_content = """\
xc  GGA  # LDA, GGA
"""
        self.load_string(file_content)