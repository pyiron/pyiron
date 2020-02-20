# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import ast
import numpy as np

"""
General purpose output parser
"""

__author__ = "Joerg Neugebauer"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


def extract_data_from_str_lst(str_lst, tag, num_args=1):
    """
    General purpose routine to extract any static from a log (text) file

    Args:
        file_name (str): file name or path to the file, can either be absolute or relative
        tag (str): string at the beginning of the line
        num_args (int): number of arguments separated by ' ' or ',' to extract after the tag

    Returns:
        list: List of arguments extracted as strings
    """

    def multiple_delimiter_split(s, seps):
        res = [s]
        for sep in seps:
            s, res = res, []
            for seq in s:
                res += seq.split(sep)
        while "" in res:
            res.remove("")
        return res

    collector = []
    ind_start = len(tag.split())
    for line_in_file in str_lst:
        if line_in_file.startswith(tag):
            collector = []
            vals = multiple_delimiter_split(line_in_file, (" ", ","))
            if num_args == 1:
                collector.append(vals[ind_start])
            else:
                collector.append(vals[ind_start : num_args + ind_start])

    return collector


def extract_data_from_file(file_name, tag, num_args=1):
    """
    General purpose routine to extract any static from a log (text) file

    Args:
        file_name (str): file name or path to the file, can either be absolute or relative
        tag (str): string at the beginning of the line
        num_args (int): number of arguments separated by ' ' or ',' to extract after the tag

    Returns:
        list: List of arguments extracted as strings
    """
    with open(file_name) as infile:
        content = infile.readlines()
    return extract_data_from_str_lst(str_lst=content, tag=tag, num_args=num_args)


class Logstatus(object):
    """
    Generic Parser for parsing output files by searching for a specific pattern structure and extracting the data that
    follows the pattern into the status_dict dictionary.

    Args:
        iter_levels (int): Levels of iteration - default = 1
    """

    def __init__(self, h5=None, iter_levels=1):  # path = None, # path of h5 file
        if h5 is not None:
            h5.add_group("generic")
            h5.move_up()
            self.h5 = h5
            self.h5_group_data = h5.getGroup().logStatus

        self.status_dict = {}
        self.iter_levels = iter_levels
        self.iter = iter_levels * [0]
        self.store_as_vector = []
        self.h5_open = False

    def reset_iter(self, dim=0):
        """
        Reset iteration level

        Args:
            dim (int): reset value - default = 0
        """
        for i in range(dim, self.iter_levels):
            self.iter[i] = 0

    def raise_iter(self, dim=0):
        """
        Increase the iteration level

        Args:
            dim (int): position - default = 0
        """
        self.iter[dim] += 1

    def append(self, title, data_to_append, vec=False):
        """
        Append data to the LogStatus object status_dict dictionary

        Args:
            title (str): Title of the data to append
            data_to_append (list,dict): the data can be of various types
            vec (bool): [True/False] if the data is a single vector instead of a matrix or a tensor
        """
        if title in self.status_dict.keys():
            if vec:
                raise ValueError(
                    "For appending matrix rather than vector option needed!"
                )
            self.status_dict[title].append([list(self.iter), data_to_append])
        else:
            self.status_dict[title] = [[list(self.iter), data_to_append]]

    def to_hdf(self, hdf):
        """
        Store the LogStatus object status_dict dictionary in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 object to store the dictionary in.
        """
        for key, value in self.status_dict.items():
            if key in self.store_as_vector:
                if len(value) > 1:
                    raise ValueError(
                        "Multi-dimensional array cannot be saved as vector"
                    )
                hdf[key] = np.array(value[0][1])
            else:
                hdf[key] = np.array([val for _, val in value])

    def combine_xyz(self, x_key, y_key, z_key, combined_key, as_vector=False):
        """
        Combine three lists representing the x,y,z coordinates, by accessing them from the status_dict dictionary,
        combining them, store them under the combined_key and remove the other three keys.

        Args:
            x_key (str): key of the x coordinates
            y_key (str): key of the y coordinates
            z_key (str): key of the z coordinates
            combined_key (str): name of the combined coordinates
        """
        if (
            x_key in self.status_dict
            and y_key in self.status_dict
            and z_key in self.status_dict
        ):
            combined_lst = []
            if as_vector:
                time_x, val_x = self.status_dict[x_key][0]
                time_y, val_y = self.status_dict[y_key][0]
                time_z, val_z = self.status_dict[z_key][0]
                for val_t_x, val_t_y, val_t_z in zip(val_x, val_y, val_z):
                    combined_lst.append([time_x, [val_t_x, val_t_y, val_t_z]])
            else:
                for var_x, var_y, var_z in zip(
                    self.status_dict[x_key],
                    self.status_dict[y_key],
                    self.status_dict[z_key],
                ):
                    time_x, val_x = var_x
                    time_y, val_y = var_y
                    time_z, val_z = var_z
                    combined_lst.append(
                        [
                            time_x,
                            [
                                [val_t_x, val_t_y, val_t_z]
                                for val_t_x, val_t_y, val_t_z in zip(
                                    val_x, val_y, val_z
                                )
                            ],
                        ]
                    )
            del self.status_dict[x_key]
            del self.status_dict[y_key]
            del self.status_dict[z_key]
            self.status_dict[combined_key] = combined_lst

    def combine_mat(self, x_key, xy_key, xz_key, y_key, yz_key, z_key, combined_key):
        """
        Combine three lists representing the x,y,z coordinates, by accessing them from the status_dict dictionary,
        combining them, store them under the combined_key and remove the other three keys.

        Args:
            x_key (str): key of the x coordinates
            y_key (str): key of the y coordinates
            z_key (str): key of the z coordinates
            combined_key (str): name of the combined coordinates
        """
        if (
            x_key in self.status_dict
            and y_key in self.status_dict
            and z_key in self.status_dict
        ):
            combined_lst = []
            for var_xx, var_xy, var_xz, var_yy, var_yz, var_zz in zip(
                self.status_dict[x_key],
                self.status_dict[xy_key],
                self.status_dict[xz_key],
                self.status_dict[y_key],
                self.status_dict[yz_key],
                self.status_dict[z_key],
            ):
                time_xx, val_xx = var_xx
                time_xy, val_xy = var_xy
                time_xz, val_xz = var_xz
                time_yy, val_yy = var_yy
                time_yz, val_yz = var_yz
                time_zz, val_zz = var_zz
                combined_lst.append(
                    [
                        time_xx,
                        [
                            [
                                [var_t_xx, var_t_xy, var_t_xz],
                                [var_t_yx, var_t_yy, var_t_yz],
                                [var_t_zx, var_t_zy, var_t_zz],
                            ]
                            for var_t_xx, var_t_xy, var_t_xz, var_t_yx, var_t_yy, var_t_yz, var_t_zx, var_t_zy, var_t_zz in zip(
                                val_xx,
                                val_xy,
                                val_xz,
                                val_xy,
                                val_yy,
                                val_yz,
                                val_xz,
                                val_yz,
                                val_zz,
                            )
                        ],
                    ]
                )
            del self.status_dict[x_key]
            del self.status_dict[xy_key]
            del self.status_dict[xz_key]
            del self.status_dict[y_key]
            del self.status_dict[yz_key]
            del self.status_dict[z_key]
            self.status_dict[combined_key] = combined_lst

    def convert_unit(self, key, factor):
        if key in self.status_dict:
            return_lst = []
            for step in self.status_dict[key]:
                time, values = step
                return_lst.append([time, (np.array(values) * factor).tolist()])
            self.status_dict[key] = return_lst

    @staticmethod
    def extract_item(l_item):
        """
        Method to extract information from a single line - currently very specific for the Lammps output

        Args:
            l_item (str): line to extract information from

        Returns:
            str, list: the tag_string as string and the arguments as list
        """
        item_list = l_item.split()
        first_item = item_list[1]
        if first_item == "NUMBER":
            num_elements = 3
        elif first_item == "BOX":
            num_elements = 2
        else:
            num_elements = 1
        tag = item_list[1 : num_elements + 1]
        tag_string = " ".join(el for el in tag)
        if len(item_list) == num_elements + 1:
            args = None
        else:
            args = item_list[num_elements + 1 : :]
        return tag_string, args

    def extract_from_list(self, list_of_lines, tag_dict, h5_dict=None, key_dict=None):
        """
        Main function of the LogStatus class to extract data from an output file by searching for the tag dictionary

        Args:
            file_name (str): absolute path to the output file
            tag_dict (dict): Dictionary with tags/patterns as key and an additional dictionary to describe the data
                             structure. The data structure dictionary can contain the following keys:
                             - "arg": position of the argument - or dimension (":", ":,:")
                             - "type": Python data type
                             - "h5": HDF5 key to store the information
                             - "rows": number of rows from the line where the tag was found
                             - "splitTag": split the tag - [True/False]
                             - "splitArg": split the argument - [True/False]
                             - "lineSkip": skip a line
                             - "func": function to convert the data
            h5_dict (dict): Translation dictionary of output tags as keys to the tags used on the HDF5 file as values.
            key_dict (dict): Translation dictionary of python internal tags as keys to the output tags as values.
        """
        val_item = {}
        tag_vals = {}

        tag = LogTag(tag_dict, h5_dict, key_dict)

        iterate_over_lines = iter(list_of_lines)
        for line_read in iterate_over_lines:
            while True:
                if tag.is_item(line_read):  # items):
                    tag_name = tag.tag_name
                    if tag.rows() == 0:  # read single line_read
                        tag.set_item(tag_vals, self)
                        try:
                            line_read = next(iterate_over_lines)
                        except StopIteration:
                            break
                    else:
                        for _ in range(tag.line_skip()):
                            line_read = next(iterate_over_lines)
                        if isinstance(tag.rows(), str):
                            i_line = 0
                            while True:
                                try:
                                    line_read = next(iterate_over_lines)
                                except StopIteration:
                                    break
                                if line_read.find(tag.rows().strip()) > -1:
                                    break
                                if "WARNING:" in line_read:
                                    break
                                val_line = [
                                    [ast.literal_eval(l) for l in line_read.split()]
                                ]
                                if i_line == 0:
                                    val_array = np.array(val_line)
                                else:
                                    val_array = np.append(
                                        arr=val_array, values=val_line, axis=0
                                    )
                                i_line += 1

                        else:
                            for i_line in range(tag.rows()):
                                try:
                                    line_read = next(iterate_over_lines)
                                except StopIteration:
                                    break
                                val_line = [
                                    [ast.literal_eval(l) for l in line_read.split()]
                                ]
                                if i_line == 0:
                                    val_array = np.array(val_line)
                                else:
                                    val_array = np.append(
                                        arr=val_array, values=val_line, axis=0
                                    )

                        if tag.is_func():
                            val_array = tag.apply_func(val_array)

                        val_item[tag_name] = val_array
                        if np.shape(val_array) == (1, 1):
                            self.append(tag.h5(), val_array[0, 0])
                        elif tag.test_split():
                            tag_list = None
                            if tag.split_tag:
                                tag_list = tag_name.split()
                            elif tag.split_arg:
                                if "header" not in tag_dict[tag_name].keys():
                                    tag_list = tag.val_list
                                else:
                                    tag_list = tag_dict[tag_name]["header"]
                            for i, t in enumerate(tag_list):
                                if "header" not in tag_dict[tag_name].keys():
                                    self.append(
                                        tag.translate(t), np.copy(val_array[:, i])
                                    )
                                else:
                                    self.append(t, np.copy(val_array[:, i]))
                        else:
                            self.append(tag.h5(), np.copy(val_array))
                else:
                    try:
                        line_read = next(iterate_over_lines)
                    except StopIteration:
                        break

    def extract_file(self, file_name, tag_dict, h5_dict=None, key_dict=None):
        """
        Main function of the LogStatus class to extract data from an output file by searching for the tag dictionary

        Args:
            file_name (str): absolute path to the output file
            tag_dict (dict): Dictionary with tags/patterns as key and an additional dictionary to describe the data
                             structure. The data structure dictionary can contain the following keys:
                             - "arg": position of the argument - or dimension (":", ":,:")
                             - "type": Python data type
                             - "h5": HDF5 key to store the information
                             - "rows": number of rows from the line where the tag was found
                             - "splitTag": split the tag - [True/False]
                             - "splitArg": split the argument - [True/False]
                             - "lineSkip": skip a line
                             - "func": function to convert the data
            h5_dict (dict): Translation dictionary of output tags as keys to the tags used on the HDF5 file as values.
            key_dict (dict): Translation dictionary of python internal tags as keys to the output tags as values.
        """
        with open(file_name, "r") as f:
            content = f.readlines()
        self.extract_from_list(
            list_of_lines=content, tag_dict=tag_dict, h5_dict=h5_dict, key_dict=key_dict
        )


class LogTag(object):
    """
    LogTag object to parse for a specific pattern in the output file

    Args:
        tag_dict (dict): Dictionary with tags/patterns as key and an additional dictionary to describe the data
                         structure. The data structure dictionary can contain the following keys:
                         - "arg": position of the argument - or dimension (":", ":,:")
                         - "type": Python data type
                         - "h5": HDF5 key to store the information
                         - "rows": number of rows from the line where the tag was found
                         - "splitTag": split the tag - [True/False]
                         - "splitArg": split the argument - [True/False]
                         - "lineSkip": skip a line
                         - "func": function to convert the data
        h5_dict (dict): Translation dictionary of output tags as keys to the tags used on the HDF5 file as values.
        key_dict (dict): Translation dictionary of python internal tags as keys to the output tags as values.
    """

    def __init__(self, tag_dict, h5_dict=None, key_dict=None):
        self._tag_dict = None
        self._tag_first_word = None
        self._current = None
        self._dyn_tags = None
        self._key_dict = None
        self._h5_dict = None
        self._tag_name = None
        self.tag_dict = tag_dict
        self.key_dict = key_dict
        self.h5_dict = h5_dict

    @property
    def current(self):
        """
        Get the current tag

        Returns:
            dict: current tag
        """
        return self._current

    @current.setter
    def current(self, tag_name):
        """
        Set the current tag

        Args:
            tag_name (str): current tag
        """
        if tag_name not in self.tag_dict.keys():
            raise ValueError("Unknown tag_name: " + tag_name)
        self._tag_name = tag_name
        self._current = self.tag_dict[tag_name]

    @property
    def tag_name(self):
        """
        Get tag name

        Returns:
            str: tag name
        """
        return self._tag_name

    @property
    def tag_dict(self):
        """
        Get tag dictionary with tags/patterns as key and an additional dictionary to describe the data
        structure. The data structure dictionary can contain the following keys:
        - "arg": position of the argument - or dimension (":", ":,:")
        - "type": Python data type
        - "h5": HDF5 key to store the information
        - "rows": number of rows from the line where the tag was found
        - "splitTag": split the tag - [True/False]
        - "splitArg": split the argument - [True/False]
        - "lineSkip": skip a line
        - "func": function to convert the data

        Returns:
            dict: tag dictionary
        """
        return self._tag_dict

    @tag_dict.setter
    def tag_dict(self, tag_dict):
        """
        Set tag dictionary with tags/patterns as key and an additional dictionary to describe the data
        structure. The data structure dictionary can contain the following keys:
        - "arg": position of the argument - or dimension (":", ":,:")
        - "type": Python data type
        - "h5": HDF5 key to store the information
        - "rows": number of rows from the line where the tag was found
        - "splitTag": split the tag - [True/False]
        - "splitArg": split the argument - [True/False]
        - "lineSkip": skip a line
        - "func": function to convert the data

        Args:
            tag_dict (dict): tag dictionary
        """
        self._tag_dict = tag_dict
        self._tag_first_word = tuple(self.tag_dict.keys())
        self.dyn_tags = tag_dict

    @property
    def tag_first_word(self):
        """
        Get first word of the tag

        Returns:
            str: first word
        """
        return self._tag_first_word

    @property
    def dyn_tags(self):
        """
        Get dynamic tags

        Returns:
            dict: dynamic tags
        """
        return self._dyn_tags

    @dyn_tags.setter
    def dyn_tags(self, tag_dict):
        """
        Set dynamic tags

        Args:
            tag_dict (dict): tag dictionary
        """
        dyn_tags = {}
        for w in tag_dict.keys():
            items = w.split()
            if items[0][:1] == "$":
                dyn_tags[w[1:]] = w
        self._dyn_tags = dyn_tags

    @property
    def key_dict(self):
        """
        Get translation dictionary of python internal tags as keys to the output tags as values.

        Returns:
            dict: key dictionary
        """
        return self._key_dict

    @key_dict.setter
    def key_dict(self, key_dict):
        """
        Set translation dictionary of python internal tags as keys to the output tags as values.

        Args:
            key_dict (dict): key dictionary
        """
        self._key_dict = key_dict

    @property
    def h5_dict(self):
        """
        Get translation dictionary of output tags as keys to the tags used on the HDF5 file as values.

        Returns:
            dict: h5 dictionary
        """
        return self._h5_dict

    @h5_dict.setter
    def h5_dict(self, h5_dict):
        """
        Set translation dictionary of output tags as keys to the tags used on the HDF5 file as values.

        Args:
            h5_dict (dict): h5 dictionary
        """
        self._h5_dict = h5_dict

    def is_item(self, item_line, start=0):
        """
        Check if the current line - item_line - matches one of the provided tags, if that is the case set the tag to be
        the current tag and update the val_list with the corresponding values.

        Args:
            item_line (str): Line of the output file
            start (int): Character to start with when parsing the item_line - default=0

        Returns:
            bool: [True/False]
        """
        l = item_line.strip()
        if not l.startswith(
            self.tag_first_word, start
        ):  # start -> line must start with tag
            return False
        tag = None
        for tag in self.tag_first_word:
            if start == l.find(tag, start):
                break

        items = [ls.strip() for ls in l[len(tag) :].split()]
        self.current = tag
        self.val_list = items
        return True

    def get_item(self, item, default):
        """
        If item is part of the current dictionary keys the corresponding value is returned otherwise the default is
        returned.

        Args:
            item (str): dictionary key
            default (list, dict, int, float): Default value

        Returns:
            list, dict, int, float: The values connected to the key item in the current dictionary and if item is not a
                                    key in the current dictionary return the default value.
        """
        if self.current is None:
            raise ValueError("current tag not defined!")
        if item in self.current.keys():
            return self.current[item]
        else:
            return default

    def h5(self):
        """
        Translate current tag to HDF5 tag using the tag dictionary

        Returns:
            str: hdf5 key name
        """
        return self.get_item(item="h5", default=self.tag_name)

    def translate(self, item):
        """
        Translate current tag to HDF5 tag using the h5_dict dictionary

        Args:
            item (str): Python tag

        Returns:
            str: HDF5 tag
        """
        if self.h5_dict is None:
            raise ValueError("h5_dict is None!" + item)
        if item in self.h5_dict.keys():
            return self.h5_dict[item]
        else:
            raise ValueError("tag not in h5_dict: " + item)

    def arg(self):
        """
        Get tag argument

        Returns:
            str: tag arguments
        """
        l_arg = self.get_item(item="arg", default=0)
        if isinstance(l_arg, str):
            return l_arg
        else:
            return str(l_arg)

    def line_skip(self):
        """
        Check how many lines should be skipped.

        Returns:
            bool: [True/ False]
        """
        return bool(self.get_item(item="lineSkip", default=0))

    def rows(self):
        """
        Number of rows to parse

        Returns:
            int, str: number of rows
        """
        rows = self.get_item(item="rows", default=0)
        try:
            return int(rows)
        except ValueError:
            return rows

    def test_split(self):
        """
        Check if the argument or the tag should be split - if "splitArg" or "splitTag" is included in the tag_dict
        dictionary.

        Returns:
            bool: [True/ False]
        """
        self.split_arg = self.get_item(item="splitArg", default=False)
        self.split_tag = self.get_item(item="splitTag", default=False)
        return self.split_arg or self.split_tag

    def is_func(self):
        """
        Check if a function is defined to convert the data - if "func" is included in the tag_dict dictionary

        Returns:
            bool: [True/ False]
        """
        my_func = self.get_item(item="func", default=None)
        return my_func is not None

    def apply_func(self, val):
        """
        Apply the function on a given value

        Args:
            val (dict, list, float, int): value to apply the function on

        Returns:
            dict, list, float, int: result of applying the function
        """
        my_func = self.get_item(item="func", default=None)
        if my_func is not None:
            return my_func(val)

    def set_item(self, tag_vals, log_file):
        """
        Set LogTag item

        Args:
            tag_vals (dict): tag value dictionary
            log_file (Logstatus): Logstatus object

        Returns:
            list: tag name, tag values, rows, line skip [True/False]
        """
        tag_name = self.tag_name
        if self.rows() == 0:
            if not len(self.arg()) == 1:
                val = []
                for i_item in ast.literal_eval(self.arg()):
                    val.append(ast.literal_eval("self.val_list[" + i_item + "]"))
            else:  # input is an array
                val = eval("self.val_list[" + self.arg() + "]")
            if isinstance(val, str):
                val = ast.literal_eval(val)
            tag_vals[tag_name] = val
            if len(self.arg()) == 1:
                log_file.append(self.h5(), data_to_append=val)
            else:
                for i_num, i_val in enumerate(val):
                    log_file.append(self.h5()[i_num], data_to_append=i_val)
            if tag_name in self.dyn_tags.keys():
                self.resolve_dynamic_variable(val)
        return tag_name, tag_vals, self.rows(), self.line_skip()

    def resolve_dynamic_variable(self, val):
        """
        Resolve dynamic variable using the key_dict dictionary

        Args:
            val: values to resolve
        """
        d_name = self.dyn_tags[self.tag_name]
        if self.key_dict is not None:
            val = [self.key_dict[v] for v in val if v in self.key_dict.keys()]
        resolved_name = " ".join(val)

        v = self.tag_dict[d_name]
        self.tag_dict[resolved_name] = v
        del self.tag_dict[d_name]
        self.dyn_tags = self.tag_dict
        self._tag_first_word = tuple(self.tag_dict.keys())
