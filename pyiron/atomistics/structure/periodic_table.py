# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, unicode_literals
import numpy as np
import os
from pyiron_base.settings.generic import Settings
from mendeleev import element
import sys
import pandas

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Martin Boeckmann"
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
pandas.options.mode.chained_assignment = None


class ChemicalElement(object):
    """
    An Object which contains the element specific parameters
    """

    def __init__(self, sub):
        """
        Constructor: assign PSE dictionary to object
        """
        self._dataset = None
        self.sub = sub
        self._mendeleev_element = None
        self._mendeleev_property_lst = None
        stringtypes = str
        if isinstance(self.sub, stringtypes):
            self._init_mendeleev(self.sub)
        elif "Parent" in self.sub.index and isinstance(self.sub.Parent, stringtypes):
            self._init_mendeleev(self.sub.Parent)
        elif len(self.sub) > 0:
            self._init_mendeleev(self.sub.Abbreviation)

        self._mendeleev_translation_dict = {
            'AtomicNumber': 'atomic_number',
            'AtomicRadius': 'covalent_radius_cordero',
            'AtomicMass': 'mass',
            'Color': 'cpk_color',
            'CovalentRadius': 'covalent_radius',
            'CrystalStructure': 'lattice_structure',
            'Density': 'density',
            'DiscoveryYear': 'discovery_year',
            'ElectronAffinity': 'electron_affinity',
            'Electronegativity': 'electronegativity',
            'Group': 'group_id',
            'Name': 'name',
            'Period': 'period',
            'StandardName': 'name',
            'VanDerWaalsRadius': 'vdw_radius',
            'MeltingPoint': 'melting_point'
        }
        self.el = None

    def _init_mendeleev(self, element_str):
        self._mendeleev_element = element(str(element_str))
        self._mendeleev_property_lst = [s for s in dir(self._mendeleev_element) if not s.startswith('_')]

    def __getattr__(self, item):
        return self[item]

    def __getitem__(self, item):
        if item in self._mendeleev_translation_dict.keys():
            item = self._mendeleev_translation_dict[item]
        if item in self._mendeleev_property_lst:
            return getattr(self._mendeleev_element, item)
        if item in self.sub.index:
            return self.sub[item]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            conditions = list()
            conditions.append(self.sub.to_dict() == other.sub.to_dict())
            return all(conditions)
        elif isinstance(other, (np.ndarray, list)):
            conditions = list()
            for sp in other:
                conditions.append(self.sub.to_dict() == sp.sub.to_dict())
            return any(conditions)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        if self != other:
            if self["AtomicNumber"] != other["AtomicNumber"]:
                return self["AtomicNumber"] > other["AtomicNumber"]
            else:
                return self["Abbreviation"] > other["Abbreviation"]
        else:
            return False

    def __ge__(self, other):
        if self != other:
            return self > other
        else:
            return True

    def __hash__(self):
        return hash(repr(self))

    @property
    def tags(self):
        if "tags" not in self.sub.keys() or self.sub["tags"] is None:
            return dict()
        return self.sub["tags"]

    def __dir__(self):
        return list(self.sub.index) + super(ChemicalElement, self).__dir__()

    def __str__(self):
        return str([self._dataset, self.sub])

    def add_tags(self, tag_dic):
        """
        Add tags to an existing element inside its specific panda series without overwriting the old tags

        Args:
            tag_dic (dict): dictionary containing e.g. key = "spin" value = "up",
                            more than one tag can be added at once

        """
        (self.sub["tags"]).update(tag_dic)

    def to_hdf(self, hdf):
        """
        saves the element with his parameters into his hdf5 job file
        Args:
            hdf (Hdfio): Hdfio object which will be used
        """
        with hdf.open(self.Abbreviation) as hdf_el:  # "Symbol of the chemical element"
            # TODO: save all parameters that are different from the parent (e.g. modified mass)
            if self.Parent is not None:
                self._dataset = {"Parameter": ["Parent"], "Value": [self.Parent]}
                hdf_el["elementData"] = self._dataset
            with hdf_el.open(
                "tagData"
            ) as hdf_tag:  # "Dictionary of element tag static"
                for key in self.tags.keys():
                    hdf_tag[key] = self.tags[key]

    def from_hdf(self, hdf):
        """
        loads an element with his parameters from the hdf5 job file and store it into its specific pandas series
        Args:
            hdf (Hdfio): Hdfio object which will be used to read a hdf5 file
        """
        pse = PeriodicTable()
        elname = self.sub.name
        with hdf.open(elname) as hdf_el:
            if "elementData" in hdf_el.list_nodes():
                element_data = hdf_el["elementData"]
                for key, val in zip(element_data["Parameter"], element_data["Value"]):
                    if key in "Parent":
                        self.sub = pse.dataframe.loc[val]
                        self.sub["Parent"] = val
                        self._init_mendeleev(val)
                    else:
                        self.sub["Parent"] = None
                        self._init_mendeleev(elname)
                    self.sub.name = elname
            if "tagData" in hdf_el.list_groups():
                with hdf_el.open(
                    "tagData"
                ) as hdf_tag:  # "Dictionary of element tag static"
                    tag_dic = {}
                    for key in hdf_tag.list_nodes():
                        tag_dic[key] = hdf_tag[key]
                        self.sub["tags"] = tag_dic


class PeriodicTable(object):
    """
    An Object which stores an elementary table which can be modified for the current session
    """

    def __init__(self, file_name=None):  # PSE_dat_file = None):
        """

        Args:
            file_name (str): Possibility to choose an source hdf5 file
        """
        self.dataframe = self._get_periodic_table_df(file_name)
        if "Abbreviation" not in self.dataframe.columns.values:
            self.dataframe["Abbreviation"] = None
        if not all(self.dataframe["Abbreviation"].values):
            for item in self.dataframe.index.values:
                if self.dataframe["Abbreviation"][item] is None:
                    self.dataframe["Abbreviation"][item] = item
        self._parent_element = None
        self.el = None

    def __getattr__(self, item):
        return self[item]

    def __getitem__(self, item):
        if item in self.dataframe.columns.values:
            return self.dataframe[item]
        if item in self.dataframe.index.values:
            return self.dataframe.loc[item]

    def from_hdf(self, hdf):
        """
        loads an element with his parameters from the hdf5 job file by creating an Object of the ChemicalElement type.
        The new element will be stored in the current periodic table.
        Changes in the tags will also be modified inside the periodic table.

        Args:
            hdf (Hdfio): Hdfio object which will be used to read the data from a hdf5 file

        Returns:

        """
        elements = hdf.list_groups()  # ["elements"]
        for el in elements:
            sub = pandas.Series()
            new_element = ChemicalElement(sub)
            new_element.sub.name = el
            new_element.from_hdf(hdf)
            new_element.sub["Abbreviation"] = el

            if "sub_tags" in new_element.tags:
                if not new_element.tags["sub_tags"]:
                    del new_element.tags["sub_tags"]

            if new_element.Parent is None:
                if not (el in self.dataframe.index.values):
                    raise AssertionError()
                if len(new_element.sub["tags"]) > 0:
                    raise ValueError("Element cannot get tag-assignment twice")
                if "tags" not in self.dataframe.keys():
                    self.dataframe["tags"] = None
                self.dataframe["tags"][el] = new_element.tags
            else:
                self.dataframe = self.dataframe.append(new_element.sub)
                self.dataframe["tags"] = self.dataframe["tags"].apply(
                    lambda x: None if pandas.isnull(x) else x
                )
                self.dataframe["Parent"] = self.dataframe["Parent"].apply(
                    lambda x: None if pandas.isnull(x) else x
                )

    def element(self, arg, **qwargs):
        """
        The method searches through the periodic table. If the table contains the element,
        it will return an Object of the type ChemicalElement containing all parameters from the periodic table.
        The option **qwargs allows a direct modification of the tag-values during the creation process
        Args:
            arg (str, ChemicalElement): sort of element
            **qwargs: e.g. a dictionary of tags

        Returns element (ChemicalElement): a element with all its properties (Abbreviation, AtomicMass, Weight, ...)

        """

        stringtypes = str
        if isinstance(arg, stringtypes):
            if arg in self.dataframe.index.values:
                self.el = arg
            else:
                raise KeyError(arg)
        elif isinstance(arg, int):

            if arg in list(self.dataframe["AtomicNumber"]):
                index = list(self.dataframe["AtomicNumber"]).index(arg)
                self.el = self.dataframe.iloc[index].name
        else:
            raise ValueError("type not defined: " + str(type(arg)))

        if qwargs is not None and "tags" not in self.dataframe.columns.values:
            self.dataframe["tags"] = None
            self.dataframe["tags"][self.el] = qwargs

        element = self.dataframe.loc[self.el]
        # element['CovalentRadius'] /= 100
        return ChemicalElement(element)

    def is_element(self, symbol):
        """
        Compares the Symbol with the Abbreviations of elements inside the periodic table
        Args:
            symbol (str): name of element, str

        Returns boolean: true for the same element, false otherwise

        """
        return symbol in self.dataframe["Abbreviation"]

    def atomic_number_to_abbreviation(self, atom_no):
        """

        Args:
            atom_no:

        Returns:

        """
        if not isinstance(atom_no, int):
            raise ValueError("type not defined: " + str(type(atom_no)))

        return self.Abbreviation[
            np.nonzero(self.AtomicNumber.to_numpy() == atom_no)[0][0]
        ]

    def add_element(
        self, parent_element, new_element, use_parent_potential=False, **qwargs
    ):
        """
        Add "additional" chemical elements to the Periodic Table. These can be used to distinguish between the various
        potentials which may exist for a given species or to introduce artificial elements such as pseudohydrogen. For
        this case set use_parent_potential = False and add in the directory containing the potential files a new file
        which is derived from the name new element.

        This function may be also used to provide additional information for the identical chemical element, e.g., to
        define a Fe_up and Fe_down to perform the correct symmetry search as well as initialization.

        Args:
            parent_element (str): name of parent element
            new_element (str): name of new element
            use_parent_potential: True: use the potential from the parent species
            **qwargs: define tags and their values, e.g. spin = "up", relax = [True, True, True]

        Returns: new element (ChemicalElement)

        """

        pandas.options.mode.chained_assignment = None
        parent_element_data_series = self.dataframe.loc[parent_element]
        parent_element_data_series["Abbreviation"] = new_element
        parent_element_data_series["Parent"] = parent_element
        parent_element_data_series.name = new_element
        if new_element not in self.dataframe.T.columns:
            self.dataframe = self.dataframe.append(parent_element_data_series)
        else:
            self.dataframe.loc[new_element] = parent_element_data_series
        if len(qwargs) != 0:
            if "tags" not in self.dataframe.columns.values:
                self.dataframe["tags"] = None
            self.dataframe["tags"][new_element] = qwargs
        if use_parent_potential:
            self._parent_element = parent_element
        return self.element(new_element)

    @staticmethod
    def _get_periodic_table_df(file_name):
        """

        Args:
            file_name:

        Returns:

        """
        if not file_name:
            for resource_path in s.resource_paths:
                if os.path.exists(os.path.join(resource_path, "atomistics")):
                    resource_path = os.path.join(resource_path, "atomistics")
                for path, folder_lst, file_lst in os.walk(resource_path):
                    for periodic_table_file_name in {"periodic_table.csv"}:
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
            raise ValueError("Was not able to locate a periodic table. ")
        else:
            if file_name.endswith(".h5"):
                return pandas.read_hdf(file_name, mode="r")
            elif file_name.endswith(".csv"):
                return pandas.read_csv(file_name, index_col=0)
            raise TypeError(
                "PeriodicTable file format not recognised: "
                + file_name
                + " supported file formats are csv, h5."
            )
