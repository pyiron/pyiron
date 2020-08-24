# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from collections import OrderedDict
import numpy as np
import os
import pandas
import posixpath
import warnings
from ast import literal_eval
from pyiron.base.settings.generic import Settings
from pyiron.base.generic.template import PyironObject

"""
GenericParameters class defines the typical input file with a key value structure plus an additional column for comments.
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

s = Settings()


class GenericParameters(PyironObject):
    """
    GenericParameters class defines the typical input file with a key value structure plus an additional column for comments.
    Convenience class to easily create, read, and modify input files

    Args:
        table_name (str): name of the input file inside the HDF5 file - optional
        input_file_name (str/Nonetype): name of the input file (if None default parameters are used)
        val_only (bool): input format consists of value (comments) only
        comment_char (str): separator that characterizes comment (e.g. "#" for python)
        separator_char (str): separator that characterizes the split between key and value - default=' '
        end_value_char (str): special character at the end of every line - default=''

    Attributes:

        .. attribute:: file_name

            file name of the input file

        .. attribute:: table_name

            name of the input table inside the HDF5 file

        .. attribute:: val_only

            boolean option to switch from a key value list to an value only input file

        .. attribute:: comment_char

            separator that characterizes comment

        .. attribute:: separator_char

            separator that characterizes the split between key and value

        .. attribute:: multi_word_separator

            multi word separator to have multi word keys

        .. attribute:: end_value_char

            special character at the end of every line

        .. attribute:: replace_char_dict

            dictionary to replace certain character combinations
    """

    def __init__(
        self,
        table_name=None,
        input_file_name=None,
        val_only=False,
        comment_char="#",
        separator_char=" ",
        end_value_char="",
    ):
        self.__name__ = "GenericParameters"
        self.__version__ = "0.1"

        self._file_name = None
        self._table_name = None
        self._comment_char = None
        self._val_only = None
        self._separator_char = None
        self._multi_word_separator = None
        self._end_value_char = None
        self._replace_char_dict = None
        self._block_dict = None
        self._bool_dict = {True: "True", False: "False"}
        self._dataset = OrderedDict()
        self._block_line_dict = {}
        self.end_value_char = end_value_char
        self.file_name = input_file_name
        self.table_name = table_name
        self.val_only = val_only
        self.comment_char = comment_char
        self.separator_char = separator_char
        self.multi_word_separator = "___"
        self.read_only = False
        if input_file_name is None:
            self.load_default()
        else:
            self.read_input(self.file_name)

    @property
    def file_name(self):
        """
        Get the file name of the input file

        Returns:
            str: file name
        """
        return self._file_name

    @file_name.setter
    def file_name(self, new_file_name):
        """
        Set the file name of the input file

        Args:
            new_file_name (str): file name
        """
        self._file_name = new_file_name

    @property
    def table_name(self):
        """
        Get the name of the input table inside the HDF5 file

        Returns:
            str: table name
        """
        return self._table_name

    @table_name.setter
    def table_name(self, new_table_name):
        """
        Set the name of the input table inside the HDF5 file

        Args:
            new_table_name (str): table name
        """
        self._table_name = new_table_name

    @property
    def val_only(self):
        """
        Get the boolean option to switch from a key value list to an value only input file

        Returns:
            bool: [True/False]
        """
        return self._val_only

    @val_only.setter
    def val_only(self, val_only):
        """
        Set the boolean option to switch from a key value list to an value only input file

        Args:
            val_only (bool): [True/False]
        """
        self._val_only = val_only

    @property
    def comment_char(self):
        """
        Get the separator that characterizes comment

        Returns:
            str: comment character
        """
        return self._comment_char

    @comment_char.setter
    def comment_char(self, new_comment_char):
        """
        Set the separator that characterizes comment

        Args:
            new_comment_char (str): comment character
        """
        self._comment_char = new_comment_char

    @property
    def separator_char(self):
        """
        Get the separator that characterizes the split between key and value

        Returns:
            str: separator character
        """
        return self._separator_char

    @separator_char.setter
    def separator_char(self, new_separator_char):
        """
        Set the separator that characterizes the split between key and value

        Args:
            new_separator_char (str): separator character
        """
        self._separator_char = new_separator_char

    @property
    def multi_word_separator(self):
        """
        Get the multi word separator to have multi word keys

        Returns:
            str: multi word separator
        """
        return self._multi_word_separator

    @multi_word_separator.setter
    def multi_word_separator(self, new_multi_word_separator):
        """
        Set the multi word separator to have multi word keys

        Args:
            new_multi_word_separator (str): multi word separator
        """
        self._multi_word_separator = new_multi_word_separator

    @property
    def end_value_char(self):
        """
        Get the special character at the end of every line

        Returns:
            str: end of line character
        """
        return self._end_value_char

    @end_value_char.setter
    def end_value_char(self, new_end_value_char):
        """
        Set the special character at the end of every line

        Args:
            new_end_value_char (str): end of line character
        """
        self._end_value_char = new_end_value_char

    @property
    def replace_char_dict(self):
        """
        Get the dictionary to replace certain character combinations

        Returns:
            dict: character replace dictionary
        """
        return self._replace_char_dict

    @replace_char_dict.setter
    def replace_char_dict(self, new_replace_char_dict):
        """
        Set the dictionary to replace certain character combinations

        Args:
            new_replace_char_dict (dict): character replace dictionary
        """
        self._replace_char_dict = new_replace_char_dict

    def _read_only_check_dict(self, new_dict):
        if self.read_only and new_dict != self._dataset:
            self._read_only_error()

    @staticmethod
    def _read_only_error():
        warnings.warn(
            "The input in GenericParameters changed, while the state of the job was already finished."
        )

    def load_string(self, input_str):
        """
        Load a multi line string to overwrite the current parameter settings

        Args:
            input_str (str): multi line string
        """
        new_dict = self._lines_to_dict(input_str.splitlines())
        self._read_only_check_dict(new_dict=new_dict)
        self._dataset = new_dict

    def load_default(self):
        """
        Load defaults resets the dataset in the background to be empty
        """
        new_dict = OrderedDict()
        new_dict["Parameter"] = []
        new_dict["Value"] = []
        new_dict["Comment"] = []
        self._read_only_check_dict(new_dict=new_dict)
        self._dataset = new_dict

    def keys(self):
        """
        Return keys of GenericParameters object
        """
        if self.val_only:
            return []
        else:
            return self._dataset["Parameter"]

    def read_input(self, file_name, ignore_trigger=None):
        """
        Read input file and store the data in GenericParameters - this overwrites the current parameter settings

        Args:
            file_name (str): absolute path to the input file
            ignore_trigger (str): trigger for lines to be ignored
        """
        Settings().logger.debug("file: %s %s", file_name, os.path.isfile(file_name))
        if not os.path.isfile(file_name):
            raise ValueError("file does not exist: " + file_name)
        with open(file_name, "r") as f:
            lines = f.readlines()
            new_lines = np.array(lines).tolist()
            if ignore_trigger is not None:
                del_ind = list()
                for i, line in enumerate(lines):
                    line = line.strip()
                    if len(line.split()) > 0:
                        if ignore_trigger == line.strip()[0]:
                            del_ind.append(i)
                        elif ignore_trigger in line:
                            lines[i] = line[: line.find("!")]
                lines = np.array(lines)
                new_lines = lines[np.setdiff1d(np.arange(len(lines)), del_ind)]
        new_dict = self._lines_to_dict(new_lines)
        self._read_only_check_dict(new_dict=new_dict)
        self._dataset = new_dict

    def get_pandas(self):
        """
        Output the GenericParameters object as Pandas Dataframe for human readability.

        Returns:
            pandas.DataFrame: Pandas Dataframe of the GenericParameters object
        """
        return pandas.DataFrame(self._dataset)

    def get(self, parameter_name, default_value=None):
        """
        Get the value of a specific parameter from GenericParameters - if the parameter is not available return
        default_value if that is set.

        Args:
            parameter_name (str): parameter key
            default_value (str): default value to return is the parameter is not set

        Returns:
            str: value of the parameter
        """
        i_line = self._find_line(parameter_name)
        if i_line > -1:
            val = self._dataset["Value"][i_line]
            try:
                val_v = literal_eval(val)
            except (ValueError, SyntaxError):
                val_v = val
            if callable(val_v):
                val_v = val
            return val_v
        elif default_value is not None:
            return default_value
        else:
            raise NameError("parameter not found: " + parameter_name)

    def get_attribute(self, attribute_name):
        """
        Get the value of a specific parameter from GenericParameters

        Args:
            attribute_name (str): parameter key

        Returns:
            str: value of the parameter
        """
        if "_attributes" not in dir(self):
            return None
        i_line = np.where(np.array(self._attributes["Parameter"]) == attribute_name)[0]
        if i_line > -1:
            return self._attributes["Value"][i_line]
        else:
            return None

    def modify(self, separator=None, append_if_not_present=False, **modify_dict):
        """
        Modify values for existing parameters. The command is called as modify(param1=val1, param2=val2, ...)

        Args:
            separator (str): needed if the parameter name contains special characters such as par:
                       use then as input: modify(separator=":", par=val) - optional
            append_if_not_present (bool): do not raise an exception but append the parameter in practice use set(par=val)
                                          - default=False
            **modify_dict (dict): dictionary of parameter names and values
        """
        # print ("modify: ", modify_dict)
        if separator is not None:
            modify_dict = {k + separator: v for k, v in modify_dict.items()}

        for key, val in modify_dict.items():
            i_key = self._find_line(key)
            if i_key == -1:
                if append_if_not_present:
                    self._append(**{key: val})
                    continue
                else:
                    raise ValueError("key for modify not found " + key)
            if isinstance(val, tuple):
                val, comment = val
                if self.read_only and self._dataset["Comment"][i_key] != comment:
                    self._read_only_error()
                self._dataset["Comment"][i_key] = comment
            if self.read_only and str(self._dataset["Value"][i_key]) != str(val):
                self._read_only_error()
            self._dataset["Value"][i_key] = str(val)

    def set(self, separator=None, **set_dict):
        """
        Set the value of multiple parameters or create new parameter key, if they do not exist already.

        Args:
            separator (float/int/str): separator string - optional
            **set_dict (dict): dictionary containing the parameter keys and their corresponding values to be set
        """
        self.modify(separator=separator, append_if_not_present=True, **set_dict)

    def set_value(self, line, val):
        """
        Set the value of a parameter in a specific line

        Args:
            line (float/int/str): line number - starting with 0
            val (str/bytes): value to be set
        """
        if line < len(self._dataset["Value"]):
            if self.read_only and self._dataset["Value"][line] != val:
                self._read_only_error()
            self._dataset["Value"][line] = val
        elif line >= len(self._dataset["Value"]):
            new_array = []
            new_comments = []
            new_params = []
            for el in self._dataset["Value"]:
                new_array.append(el)
                new_comments.append("")
                new_params.append("")
            new_array.append(val)
            new_comments.append("")
            new_params.append("")
            new_dict = OrderedDict()
            new_dict["Value"] = new_array
            new_dict["Comment"] = new_comments
            new_dict["Parameter"] = new_params
            self._read_only_check_dict(new_dict=new_dict)
            self._dataset = new_dict
        else:
            raise ValueError("Wrong indexing")

    def remove_keys(self, key_list):
        """
        Remove a list of keys from the GenericParameters

        Args:
            key_list (list): list of keys to be removed
        """
        if self.read_only and any([k in self._dataset["Parameter"] for k in key_list]):
            self._read_only_error()
        for key in key_list:
            params = np.array(self._dataset["Parameter"])
            i_keys = np.where(params == key)[0]
            if len(i_keys) == 0:
                continue
            if i_keys[0] == -1:
                continue
            for i_key in i_keys[::-1]:
                self._delete_line(i_key)

    def define_blocks(self, block_dict):
        """
        Define a block section within the GenericParameters

        Args:
            block_dict (dict): dictionary to define the block
        """
        if not isinstance(block_dict, OrderedDict):
            raise AssertionError()
        self._block_dict = block_dict

    def to_hdf(self, hdf, group_name=None):
        """
        Store the GenericParameters in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object
            group_name (str): HDF5 subgroup name - optional
        """
        if group_name:
            with hdf.open(group_name) as hdf_group:
                hdf_child = hdf_group.create_group(self.table_name)
        else:
            hdf_child = hdf.create_group(self.table_name)

        self._type_to_hdf(hdf_child)
        hdf_child["data_dict"] = self._dataset

    def from_hdf(self, hdf, group_name=None):
        """
        Restore the GenericParameters from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object
            group_name (str): HDF5 subgroup name - optional
        """
        if group_name:
            with hdf.open(group_name) as hdf_group:
                data = hdf_group[self.table_name]
        else:
            data = hdf[self.table_name]
        if isinstance(data, dict):
            self._dataset = data
        else:
            self._dataset = data._read("data_dict")

    def get_string_lst(self):
        """
        Get list of strings from GenericParameters to write to input file
        """
        tab_dict = self._dataset
        # assert(len(tab_dict['Value']) == len(tab_dict['Parameter']))
        if "Parameter" not in tab_dict:
            tab_dict["Parameter"] = ["" for _ in tab_dict["Value"]]

        string_lst = []
        if self.val_only:
            value_lst = tab_dict["Value"]
        else:
            try:
                value_lst = [self[p] for p in tab_dict["Parameter"]]
            except ValueError:
                value_lst = tab_dict["Value"]
        for par, v, c in zip(tab_dict["Parameter"], value_lst, tab_dict["Comment"]):
            # special treatment for values that are bool or str
            if isinstance(v, bool):
                v_str = self._bool_dict[v]
            elif isinstance(
                v, str
            ):  # TODO: introduce variable for string symbol (" or ')
                v_str = v
            else:
                v_str = str(v)

            par = " ".join(par.split(self.multi_word_separator))
            if par == "Comment":
                string_lst.append(str(v) + self.end_value_char + "\n")
            elif c.strip() == "":
                if self.val_only:
                    string_lst.append(v_str + self.end_value_char + "\n")
                else:
                    string_lst.append(
                        par + self.separator_char + v_str + self.end_value_char + "\n"
                    )
            else:
                if self.val_only:
                    string_lst.append(
                        v_str + self.end_value_char + " " + self.comment_char + c + "\n"
                    )
                else:
                    string_lst.append(
                        par
                        + self.separator_char
                        + v_str
                        + " "
                        + self.end_value_char
                        + self.comment_char
                        + c
                        + "\n"
                    )
        return string_lst

    def write_file(self, file_name, cwd=None):
        """
        Write GenericParameters to input file

        Args:
            file_name (str): name of the file, either absolute (then cwd must be None) or relative
            cwd (str): path name (default: None)
        """
        if cwd is not None:
            file_name = posixpath.join(cwd, file_name)

        with open(file_name, "w") as f:
            for line in self.get_string_lst():
                f.write(line)

    def __repr__(self):
        """
        Human readable string representation

        Returns:
            str: pandas Dataframe structure as string
        """
        return str(self.get_pandas())

    def __setitem__(self, key, value):
        """
        Set a value for the corresponding key

        Args:
            key (str): key to be set of modified
            value (float/int/str): value to be set
        """
        if isinstance(key, int):
            if self.read_only and self._dataset["Value"][key] != value:
                self._read_only_error()
            self._dataset["Value"][key] = value
        else:
            self.set(**{key: value})

    def set_dict(self, dictionary):
        """
        Set a dictionary of key value pairs

        Args:
            dictionary (dict): dictionary of key value pairs
        """
        self.set(**dictionary)

    def __getitem__(self, item):
        """
        Get a value for the corresponding key

        Args:
            item (int, str): key

        Returns:
            str: value
        """
        if isinstance(item, int):
            return self._dataset["Value"][item]
        elif item in self._dataset["Parameter"]:
            return self.get(item)

    def __delitem__(self, key):
        """
        Delete a key from GenericParameters

        Args:
            key (str): single key
        """
        self.remove_keys([key])

    def _get_block(self, block_name):
        """
        Internal helper function to get a block by name

        Args:
            block_name (str): block name

        Returns:
            dict: dictionary of the specific block
        """
        if block_name not in self._block_dict:
            raise ValueError("unknown block: " + block_name)
        keys = self._dataset["Parameter"]
        block_dict = OrderedDict()
        for key in self._dataset:
            block_dict[key] = []
        for i, tag in enumerate(keys):
            if tag in self._block_dict[block_name]:
                for key in self._dataset:
                    block_dict[key].append(self._dataset[key][i])
        return block_dict

    def _get_attributes(self):
        """
        Internal helper function to extract pyiron specific commands (start in comments with " @my_command")

        Returns:
            (dict): {"Parameter": list of tags, "Value": list of values}
        """
        tags = self._dataset["Parameter"]
        lst_tag, lst_val = [], []
        for i, tag in enumerate(tags):
            if tag not in ["Comment"]:
                continue
            c = self._dataset["Value"][i]
            s_index = c.find(" @")
            if s_index > -1:
                tag, val = c[s_index:].split()[:2]
                lst_tag.append(tag[1:])
                lst_val.append(val)
        self._attributes = {"Parameter": lst_tag, "Value": lst_val}
        return self._attributes

    def _remove_block(self, block_name):
        """
        Internal helper function to remove a block by name

        Args:
            block_name (str): block name
        """
        if block_name not in self._block_dict:
            raise ValueError("unknown block to be removed")
        self.remove_keys(self._block_dict[block_name])

    def _insert_block(self, block_dict, next_block=None):
        """
        Internal helper function to insert a block by name

        Args:
            block_dict (dict): block dictionary
            next_block (str): name of the following block - optional
        """
        if next_block is None:  # append
            for key in block_dict:
                self._dataset[key] += block_dict[key]
        else:
            for i, tag in enumerate(self._dataset["Parameter"]):
                if tag in self._block_dict[next_block]:
                    self._insert(line_number=i, data_dict=block_dict)  # , shift=1)
                    break

    def _update_block(self, block_dict):
        """
        Internal helper function to update a block by name

        Args:
            block_dict (dict): block dictionary
        """
        tag_lst = block_dict["Parameter"]
        val_lst = block_dict["Value"]
        par_dict = {}
        for t, v in zip(tag_lst, val_lst):
            par_dict[t] = v
        self.modify(**par_dict)

    def _delete_line(self, line_number):
        """
        Internal helper function to delete a single line

        Args:
            line_number (int): line number
        """
        if self.read_only:
            self._read_only_error()
        for key, val in self._dataset.items():
            if "numpy" in str(type(val)):
                val = val.tolist()
            del val[line_number]
            self._dataset[key] = val

    def _insert(self, line_number, data_dict, shift=0):
        """
        Internal helper function to insert a single line by line number

        Args:
            line_number (int): line number
            data_dict (dict): data dictionary
            shift (int): shift line number - default=0
        """
        if self.read_only:
            self._read_only_error()
        for key, val in data_dict.items():
            lst = self._dataset[key]
            val = np.array(val).tolist()
            lst = np.array(lst).tolist()
            self._dataset[key] = lst[: line_number - shift] + val + lst[line_number:]

    def _refresh_block_line_hash_table(self):
        """
        Internal helper function to refresh the block dictionary hash table
        """
        self._block_line_dict = {}
        for i_line, par in enumerate(self._dataset["Parameter"]):
            if par.strip() == "":
                continue
            for key, val in self._block_dict.items():
                par_single = par.split()[0].split(self.multi_word_separator)[0]
                if par_single in val:
                    if key in self._block_line_dict:
                        self._block_line_dict[key].append(i_line)
                    else:
                        self._block_line_dict[key] = [i_line]
                    break
        i_line_old = 0
        for key in self._block_dict:
            if key in self._block_line_dict:
                i_line_old = np.max(self._block_line_dict[key])
            else:
                self._block_line_dict[key] = [i_line_old]

    def _append_line_in_block(self, parameter_name, value):
        """
        Internal helper function to append a line within a block

        Args:
            parameter_name (str): name of the parameter
            value (str): value of the parameter

        Returns:
            bool: [True/False]
        """
        for key, val in self._block_dict.items():
            par_first = parameter_name.split()[0].split(self.multi_word_separator)[0]
            if par_first in val:
                i_last_block_line = max(self._block_line_dict[key])
                self._insert(
                    line_number=i_last_block_line + 1,
                    data_dict={
                        "Parameter": [parameter_name],
                        "Value": [str(value)],
                        "Comment": [""],
                    },
                )
                return True
        else:
            s.logger.warning(
                "Unknown parameter (does not exist in block_dict): {}".format(
                    parameter_name
                )
            )
        return False

    def _append(self, **qwargs):
        """
        Internal helper function to append data to the GenericParameters object

        Args:
            **qwargs (dict): dictionary with parameter keys and their corresponding values
        """
        if self.read_only:
            self._read_only_error()
        for par, val in qwargs.items():
            if par in self._dataset["Parameter"]:
                raise ValueError("Parameter exists already: " + par)

            if self._block_dict is not None:
                self._refresh_block_line_hash_table()
                if self._append_line_in_block(par, val):
                    continue

            for col in self._dataset:
                self._dataset[col] = np.array(self._dataset[col]).tolist()

            comment = ""
            if isinstance(val, tuple):
                val, comment = val
            self._dataset["Parameter"].append(par)
            self._dataset["Value"].append(val)
            self._dataset["Comment"].append(comment)

    def _is_multi_word_parameter(self, key):
        """
        Internal helper function to check if a parameter included multiple words

        Args:
            key (str): parameter

        Returns:
            bool: [True/False]
        """
        par_list = key.split(self.multi_word_separator)
        return len(par_list) > 1

    def _repr_html_(self):
        """
        Internal helper function to represent the GenericParameters object within the Jupyter Framework

        Returns:
            HTML: Jupyter HTML object
        """
        return self.get_pandas()._repr_html_()

    def _lines_to_dict(self, lines):
        """
        Internal helper function to convert multiple lines to a dictionary

        Args:
            lines (list): list of lines

        Returns:
            dict: GenericParameters dictionary
        """
        lst = OrderedDict()
        lst["Parameter"] = []
        lst["Value"] = []
        lst["Comment"] = []
        for line in lines:
            # print ("line: ", line)
            if self.replace_char_dict is not None:
                for key, val in self.replace_char_dict.items():
                    line = line.replace(key, val)

            sep = line.split(self.comment_char)
            if len(line.strip()) > 0 and (
                line.strip()[0] == self.comment_char
            ):  # comment line
                lst["Parameter"].append("Comment")
                lst["Value"].append(self.bool_str_to_bool(line[:-1]))
                lst["Comment"].append("")
            elif not sep[0].strip() == "":
                sep[0] = sep[0].strip()
                if self.val_only:  # Value only entries
                    val = sep[0]
                    name = ""
                else:
                    keypos = sep[0].find(self.separator_char)
                    if keypos == -1:  # Key only entries
                        name = sep[0]
                        val = ""
                    else:  # Entires with key and value
                        name = sep[0][0:keypos]
                        val = sep[0][keypos + len(self.separator_char) :]
                lst["Parameter"].append(name.strip())
                lst["Value"].append(self.bool_str_to_bool(val.strip()))
                if len(sep) > 1:  # Handle comments
                    lst["Comment"].append(sep[-1].strip())
                else:
                    lst["Comment"].append("")
            else:  # Handle empty lines
                lst["Parameter"].append("")
                lst["Value"].append("")
                lst["Comment"].append("")
        return lst

    def _type_to_hdf(self, hdf):
        """
        Internal helper function to save type and version in hdf root

        Args:
            hdf (ProjectHDFio): HDF5 group object
        """
        hdf["NAME"] = self.__name__
        hdf["TYPE"] = str(type(self))
        hdf["VERSION"] = self.__version__
        hdf["OBJECT"] = "GenericParameters"

    def _find_line(self, key_name):
        """
        Internal helper function to find a line by key name

        Args:
            key_name (str): key name

        Returns:
            list: [line index, line]
        """
        params = self._dataset["Parameter"]
        if len(params) > 0:
            i_line_lst = np.where(np.array(params) == key_name)[0]
        else:
            i_line_lst = []
        if len(i_line_lst) == 0:
            return -1
        elif len(i_line_lst) == 1:
            return i_line_lst[0]
        else:
            error_msg = list()
            error_msg.append("Multiple occurrences of key_name: " + key_name + ". They are as follows")
            for i in i_line_lst:
                error_msg.append("dataset: {}, {}, {}".format(i,
                                                              self._dataset["Parameter"][i],
                                                              self._dataset["Value"][i]))
            error_msg = "\n".join(error_msg)
            raise ValueError(error_msg)

    def clear_all(self):
        """
        Clears all fields in the object
        """
        self._dataset["Parameter"] = []
        self._dataset["Value"] = []
        self._dataset["Comment"] = []

    def bool_str_to_bool(self, val):
        for key, value in self._bool_dict.items():
            if val == value:
                return key
        return val
