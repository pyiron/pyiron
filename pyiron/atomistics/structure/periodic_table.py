# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function, unicode_literals
import numpy as np
import os
from pyiron.base.settings.generic import Settings
import sys
import pandas

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Martin Boeckmann"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
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
        self.el = None

    def __getattr__(self, item):
        return self[item]

    def __getitem__(self, item):
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
        if 'tags' not in self.sub.keys() or self.sub['tags'] is None:
            return dict()
        return self.sub['tags']

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
        (self.sub['tags']).update(tag_dic)

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
            with hdf_el.open("tagData") as hdf_tag:  # "Dictionary of element tag static"
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
                    if key in 'Parent':
                        self.sub = pse.dataframe.loc[val]
                        self.sub['Parent'] = val
                    else:
                        self.sub['Parent'] = None
                    self.sub.name = elname
            if "tagData" in hdf_el.list_groups():
                with hdf_el.open("tagData") as hdf_tag:  # "Dictionary of element tag static"
                    tag_dic = {}
                    for key in hdf_tag.list_nodes():
                        tag_dic[key] = hdf_tag[key]
                        self.sub['tags'] = tag_dic


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
        if 'Abbreviation' not in self.dataframe.columns.values:
            self.dataframe['Abbreviation'] = None
        if not all(self.dataframe['Abbreviation'].values):
            for item in self.dataframe.index.values:
                if self.dataframe['Abbreviation'][item] is None:
                    self.dataframe['Abbreviation'][item] = item
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
            new_element.sub['Abbreviation'] = el

            if 'sub_tags' in new_element.tags:
                if not new_element.tags['sub_tags']:
                    del new_element.tags['sub_tags']

            if new_element.Parent is None:
                if not (el in self.dataframe.index.values):
                    raise AssertionError()
                if len(new_element.sub['tags']) > 0:
                    raise ValueError('Element cannot get tag-assignment twice')
                if 'tags' not in self.dataframe.keys():
                    self.dataframe['tags'] = None
                self.dataframe['tags'][el] = new_element.tags
            else:
                self.dataframe = self.dataframe.append(new_element.sub)
                self.dataframe['tags'] = self.dataframe['tags'].apply(lambda x: None if pandas.isnull(x) else x)
                self.dataframe['Parent'] = self.dataframe['Parent'].apply(lambda x: None if pandas.isnull(x) else x)

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

        if sys.version_info.major == 2:
            stringtypes = (str, unicode)
        else:
            stringtypes = str
        if isinstance(arg, stringtypes):
            if arg in self.dataframe.index.values:
                self.el = arg
            else:
                raise KeyError(arg)
        elif isinstance(arg, int):

            if arg in list(self.dataframe['AtomicNumber']):
                index = list(self.dataframe['AtomicNumber']).index(arg)
                self.el = self.dataframe.iloc[index].name
        else:
            raise ValueError("type not defined: " + str(type(arg)))

        if qwargs is not None and 'tags' not in self.dataframe.columns.values:
            self.dataframe['tags'] = None
            self.dataframe['tags'][self.el] = qwargs

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
        return symbol in self.dataframe['Abbreviation']

    def atomic_number_to_abbreviation(self, atom_no):
        """
        
        Args:
            atom_no: 

        Returns:

        """
        if not isinstance(atom_no, int):
            raise ValueError("type not defined: " + str(type(atom_no)))

        return self.Abbreviation[np.nonzero(self.AtomicNumber == atom_no)[0][0]]

    def add_element(self, parent_element, new_element, use_parent_potential=False, **qwargs):
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
        parent_element_data_series['Abbreviation'] = new_element
        parent_element_data_series['Parent'] = parent_element
        parent_element_data_series.name = new_element
        if new_element not in self.dataframe.T.columns:
            self.dataframe = self.dataframe.append(parent_element_data_series)
        else:
            self.dataframe.loc[new_element] = parent_element_data_series
        if len(qwargs) != 0:
            if 'tags' not in self.dataframe.columns.values:
                self.dataframe['tags'] = None
            self.dataframe['tags'][new_element] = qwargs
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
                if os.path.exists(os.path.join(resource_path, 'atomistics')):
                    resource_path = os.path.join(resource_path, 'atomistics')
                for path, folder_lst, file_lst in os.walk(resource_path):
                    for periodic_table_file_name in {'periodic_table.csv'}:
                        if periodic_table_file_name in file_lst and periodic_table_file_name.endswith('.csv'):
                            return pandas.read_csv(os.path.join(path, periodic_table_file_name), index_col=0)
                        elif periodic_table_file_name in file_lst and periodic_table_file_name.endswith('.h5'):
                            return pandas.read_hdf(os.path.join(path, periodic_table_file_name), mode='r')
            raise ValueError('Was not able to locate a periodic table. ')
        else:
            if file_name.endswith('.h5'):
                return pandas.read_hdf(file_name, mode='r')
            elif file_name.endswith('.csv'):
                return pandas.read_csv(file_name, index_col=0)
            raise TypeError("PeriodicTable file format not recognised: " + file_name +
                            " supported file formats are csv, h5.")


class ElementColorDictionary(object):
    """

    """
    elementColors = {'H': [1, 255, 255, 255, 255],
                     'He': [2, 217, 255, 255, 255],
                     'Li': [3, 204, 128, 255, 255],
                     'Be': [4, 194, 255, 0, 255],
                     'B': [5, 255, 181, 181, 255],
                     'C': [6, 144, 144, 144, 255],
                     'N': [7, 48, 80, 248, 255],
                     'O': [8, 255, 13, 13, 255],
                     'F': [9, 144, 224, 80, 255],
                     'Ne': [10, 179, 227, 245, 255],
                     'Na': [11, 171, 92, 242, 255],
                     'Mg': [12, 138, 255, 0, 255],
                     'Al': [13, 191, 166, 166, 255],
                     'Si': [14, 240, 200, 160, 255],
                     'P': [15, 255, 128, 0, 255],
                     'S': [16, 255, 255, 48, 255],
                     'Cl': [17, 31, 240, 31, 255],
                     'Ar': [18, 128, 209, 227, 255],
                     'K': [19, 143, 64, 212, 255],
                     'Ca': [20, 61, 255, 0, 255],
                     'Sc': [21, 230, 230, 230, 255],
                     'Ti': [22, 191, 194, 199, 255],
                     'V': [23, 166, 166, 171, 255],
                     'Cr': [24, 138, 153, 199, 255],
                     'Mn': [25, 156, 122, 199, 255],
                     'Fe': [26, 224, 102, 51, 255],
                     'Co': [27, 240, 144, 160, 255],
                     'Ni': [28, 80, 208, 80, 255],
                     'Cu': [29, 200, 128, 51, 255],
                     'Zn': [30, 125, 128, 176, 255],
                     'Ga': [31, 194, 143, 143, 255],
                     'Ge': [32, 102, 143, 143, 255],
                     'As': [33, 189, 128, 227, 255],
                     'Se': [34, 255, 161, 0, 255],
                     'Br': [35, 166, 41, 41, 255],
                     'Kr': [36, 92, 184, 209, 255],
                     'Rb': [37, 112, 46, 176, 255],
                     'Sr': [38, 0, 255, 0, 255],
                     'Y': [39, 148, 255, 255, 255],
                     'Zr': [40, 148, 224, 224, 255],
                     'Nb': [41, 115, 194, 201, 255],
                     'Mo': [42, 84, 181, 181, 255],
                     'Tc': [43, 59, 158, 158, 255],
                     'Ru': [44, 36, 143, 143, 255],
                     'Rh': [45, 10, 125, 140, 255],
                     'Pd': [46, 0, 105, 133, 255],
                     'Ag': [47, 192, 192, 192, 255],
                     'Cd': [48, 255, 217, 143, 255],
                     'In': [49, 166, 117, 115, 255],
                     'Sn': [50, 102, 128, 128, 255],
                     'Sb': [51, 158, 99, 181, 255],
                     'Te': [52, 212, 122, 0, 255],
                     'I': [53, 148, 0, 148, 255],
                     'Xe': [54, 66, 158, 176, 255],
                     'Cs': [55, 87, 23, 143, 255],
                     'Ba': [56, 0, 201, 0, 255],
                     'La': [57, 112, 212, 255, 255],
                     'Ce': [58, 255, 255, 199, 255],
                     'Pr': [59, 217, 255, 199, 255],
                     'Nd': [60, 199, 255, 199, 255],
                     'Pm': [61, 163, 255, 199, 255],
                     'Sm': [62, 143, 255, 199, 255],
                     'Eu': [63, 97, 255, 199, 255],
                     'Gd': [64, 69, 255, 199, 255],
                     'Tb': [65, 48, 255, 199, 255],
                     'Dy': [66, 31, 255, 199, 255],
                     'Ho': [67, 0, 255, 156, 255],
                     'Er': [68, 0, 230, 117, 255],
                     'Tm': [69, 0, 212, 82, 255],
                     'Yb': [70, 0, 191, 56, 255],
                     'Lu': [71, 0, 171, 36, 255],
                     'Hf': [72, 77, 194, 255, 255],
                     'Ta': [73, 77, 166, 255, 255],
                     'W': [74, 33, 148, 214, 255],
                     'Re': [75, 38, 125, 171, 255],
                     'Os': [76, 38, 102, 150, 255],
                     'Ir': [77, 23, 84, 135, 255],
                     'Pt': [78, 208, 208, 224, 255],
                     'Au': [79, 255, 209, 35, 255],
                     'Hg': [80, 184, 184, 208, 255],
                     'Tl': [81, 166, 84, 77, 255],
                     'Pb': [82, 87, 89, 97, 255],
                     'Bi': [83, 158, 79, 181, 255],
                     'Po': [84, 171, 92, 0, 255],
                     'At': [85, 117, 79, 69, 255],
                     'Rn': [86, 66, 130, 150, 255],
                     'Fr': [87, 66, 0, 102, 255],
                     'Ra': [88, 0, 125, 0, 255],
                     'Ac': [89, 112, 171, 250, 255],
                     'Th': [90, 0, 186, 255, 255],
                     'Pa': [91, 0, 161, 255, 255],
                     'U': [92, 0, 143, 255, 255],
                     'Np': [93, 0, 128, 255, 255],
                     'Pu': [94, 0, 107, 255, 255],
                     'Am': [95, 84, 92, 242, 255],
                     'Cm': [96, 120, 92, 227, 255],
                     'Bk': [97, 138, 79, 227, 255],
                     'Cf': [98, 161, 54, 212, 255],
                     'Es': [99, 179, 31, 212, 255],
                     'Fm': [100, 179, 31, 186, 255],
                     'Md': [101, 179, 13, 166, 255],
                     'No': [102, 189, 13, 135, 255],
                     'Lr': [103, 199, 0, 102, 255],
                     'Rf': [104, 204, 0, 89, 255],
                     'Db': [105, 209, 0, 79, 255],
                     'Sg': [106, 217, 0, 69, 255],
                     'Bh': [107, 224, 0, 56, 255],
                     'Hs': [108, 230, 0, 46, 255],
                     'Mt': [109, 235, 0, 38, 255]
                     }

    def to_lut(self):
        """

        Returns:

        """
        rv = np.zeros((256, 4), dtype=int)
        if sys.version_info.major > 2:
            for k, el in self.elementColors.items():
                rv[el[0], :] = np.array(el[1:5])
        else:
            for k, el in self.elementColors.iteritems():
                rv[el[0], :] = np.array(el[1:5])
        return rv



