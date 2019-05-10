# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
# import os
import sys
import copy
import numpy as np
from collections import OrderedDict, Sequence

__author__ = "Joerg Neugebauer"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class SparseListElement(object):
    """
    Handle single element of a sparse lisr
    Args:
        ind: index
        val: value
    """
    def __init__(self, ind, val):
        self.index = ind
        self.value = val

    def __str__(self):
        return '({}: {})'.format(self.index, self.value)


class SparseList(object):
    """
    Object to represent a single sparse list
    Internal representation like a dict
    External representation like a list
    Args:
        sparse_list: dict object with {index: val}
        default: default value for all elements not given by index in sparse_list
        length: length of the list
    """
    def __init__(self, sparse_list, default=None, length=None):
        if isinstance(sparse_list, dict):
            self._dict = sparse_list.copy()
            if "_" in self._dict.keys():
                default = self._dict["_"]
                del self._dict["_"]

            if length is None:
                raise ValueError('Length must be provided in dict input mode')
            self._length = length
        elif isinstance(sparse_list, (list, np.ndarray)):
            # self._dict = {el: [] for el in set(sparse_list)}
            self._dict = {}
            for i, el in enumerate(sparse_list):
                self._dict[i] = el
            self._length = len(sparse_list)
            if length is not None:
                if length != self._length:
                    raise ValueError('Incompatible length of new list')
        self._default = default

    def _val_data_type(self):
        """

        Returns:

        """
        if isinstance(self.values(), dict):
            pass
            print(self.values())
        data_0 = self.values()[0]
        if isinstance(data_0, list):
            if isinstance(data_0[0], bool):
                return "list_bool"
            else:
                raise ValueError('tags which have as elements lists or tensors are not implemented')
        else:
            return "scalar"

    def to_hdf(self, hdf, key):
        """

        Args:
            hdf:
            key:

        Returns:

        """
        if len(self.list()) > 0:
            # Convert to array and store
            hdf[key] = np.array(self.list())
        elif len(self.values()) > 0:
            print('sparse array: ', key, len(self.values()))
            data_type = self._val_data_type()
            my_dict = OrderedDict()
            my_dict["index"] = self.keys()
            if data_type is "list_bool":
                my_dict["values"] = [sum([2 ** i * int(v) for i, v in enumerate(val)]) for val in self.values()]
            else:
                my_dict["values"] = self.values()
            print("values: ", self.values())
            hdf[key] = my_dict

    def __len__(self):
        return self._length

    def __copy__(self):
        return SparseList(sparse_list=self._dict, default=self._default, length=self._length)

    def keys(self):
        """
        
        Returns:
            indices of non-sparse elements
        """
        return self._dict.keys()

    def values(self):
        """

        Returns:
            values of non-sparse elements
        """
        return self._dict.values()

    def items(self):
        """

        Returns:
            index, value pairs of non-sparse elements
        """
        return self._dict.items()

    def list(self):
        """
        convert sparse list into full list
        Returns:
            list representation
        """
        full_list = [self._default for _ in range(self._length)]
        for i, val in self._dict.items():
            full_list[i] = val

        return full_list

    def __iter__(self):
        if self._default is None:
            for i, val in self._dict.items():
                yield SparseListElement(i, val)
        else:
            for i, val in enumerate(self.list()):
                yield val

    def __getitem__(self, item):
        if isinstance(item, (int, np.integer)):
            if item in self._dict:
                return self._dict[item]
            return self._default

        if isinstance(item, slice):
            ind_list = range(len(self))[item]
        elif isinstance(item, (list, tuple, np.ndarray)):
            if len(item) == 0:
                ind_list = []
            else:
                if isinstance(item[0], (int, np.integer)):
                    ind_list = item
                elif isinstance(item[0], (bool, np.bool_)):
                    ind_list = []
                    for i, bo in enumerate(item):
                        if bo:
                            ind_list.append(i)
        else:
            raise ValueError('Unknown item type: ' + str(type(item)))
        sliced_dict = {j: self._dict[ind] for j, ind in enumerate(ind_list) if ind in self._dict}

        return self.__class__(sliced_dict, default=self._default, length=len(ind_list))

    def __setitem__(self, key, value):
        if isinstance(key, (int, np.integer)):
            if key > len(self):
                raise IndexError
            self._dict[key] = value
            return
        elif isinstance(key, slice):
            key = range(len(self))[key]

        if max(key) > self._length:
            raise IndexError
        for i in key:
            self._dict[i] = value

    def __delitem__(self, key):
        # programmed for simplicity, not for performance
        ind_list = list(range(len(self)))
        if isinstance(key, (list, np.ndarray, tuple)):
            indexes = sorted(list(key), reverse=True)
            for index in indexes:
                del ind_list[index]
        else:
            del ind_list[key]
        new_list = self[ind_list]
        self._dict = new_list._dict
        self._length = new_list._length
        self._default = new_list._default

    def __add__(self, other):
        if not (isinstance(other, SparseList)):
            raise AssertionError()
        if not (self._default == other._default):
            raise AssertionError()
        new_list = self.__copy__()
        shifted_dict = {i + self._length: val for i, val in other._dict.items()}
        new_list._dict.update(shifted_dict)
        new_list._length += len(other)
        return new_list

    def __mul__(self, other):
        if not isinstance(other, (int, np.integer)):
            raise ValueError('Multiplication defined only for SparseArray*integers')
        overall_list = other * np.arange(len(self)).tolist()
        new_dic = dict()
        for k in self.keys():
            for val in np.argwhere(np.array(overall_list) == k).flatten():
                new_dic[val] = self[k]
        return self.__class__(new_dic, default=self._default, length=other * len(self))

    def __rmul__(self, other):
        if isinstance(other, int):
            return self * other

    def __str__(self):
        if self._default is None:
            return "[" + " ".join([str(el) for el in self]) + "]"
        else:
            # return "[" + " ".join([str(el) + os.sep for el in self.list()]) + "]"
            return "[" + " ".join([str(el) for el in self.list()]) + "]"

    def __repr__(self):
        return str(self.list())


def sparse_index(index_list, length, default_val=True):
    """

    Args:
        index_list:
        length:
        default_val:

    Returns:

    """
    new_dict = {i: default_val for i in index_list}
    return SparseList(new_dict, length=length)


class SparseArrayElement(object):
    """
    Single element of a SparseArray 
    Args:
        **qwargs:
    """
    def __init__(self, **qwargs):
        self._lists = dict()
        if qwargs:
            self._lists = qwargs

    def __getattr__(self, item):
        if item in self._lists.keys():
            return self._lists[item]
        raise AttributeError('Object has no attribute {} {}'.format(self.__class__, item))

    def __str__(self):
        out_str = ""
        for key, val in self._lists.items():
            out_str += '{}: {}'.format(key, val)
        return out_str

    def __eq__(self, other):
        if not (isinstance(other, SparseArrayElement)):
            raise AssertionError()
        conditions = []
        for key in self._lists.keys():
            try:
                if isinstance(self._lists[key], np.ndarray):
                    conditions += list(np.equal(self._lists[key], other._lists[key]))
                else:
                    conditions.append(self._lists[key] == other._lists[key])
            except KeyError:
                conditions.append(False)
        return all(conditions)


class SparseArray(object):
    """
    Administrate object that consists of several sparse lists (tags) and full lists that have identical indices and
    length

    Args:
        **qwargs: dictionary containing lists and SparseLists (tags) (must have identical length)
    """
    def __init__(self, length=None, **qwargs):
        self._lists = dict()
        self._length = length
        for key in qwargs:
            value = qwargs[key]
            if self._length is None:
                self._length = len(value)
            else:
                if not len(self) == len(value):
                    raise ValueError('Inconsistent vector lengths {} {} {}'.format(key, len(self), len(value)))
            self._lists[key] = value

    def __setitem__(self, key, value):
        # exclude hidden variables (starting with _ from being added to _lists
        # if (not hasattr(self, '_lists')) or (key[0] == "_"):
        #     self.__dict__[key] = value
        #     return
        # el

        if isinstance(value, SparseList):
            self._lists[key] = value
            return
        elif isinstance(value, (Sequence, np.ndarray)):
            if len(value) == len(self):
                self._lists[key] = value
                return
            else:
                raise ValueError(
                    'Length of array object and new list are inconsistent: {} {} {}'.format(key, len(value), len(self)))
        raise ValueError('Unsupported argument: ' + str(type(value)))

    def __getattr__(self, item):
        # if not (item in ["_lists"]):
        #     print "item: ", item, hasattr(self, item)
        if sys.version_info.major > 2:
            if '_lists' in dir(self):  # Python 3
                if item in self._lists.keys():
                    return self._lists[item]
        else:
            if hasattr(self, '_lists'):
                if item in self._lists.keys():
                    return self._lists[item]

        return object.__getattribute__(self, item)
        # raise AttributeError("%r object has no attribute %r" %(self.__class__, item))

    def __delitem__(self, key):
        for k in self.keys():
            if len(self._lists[k]) == 0:
                # ensure ASE compatibility
                print('Empty key in SparseList: ', k, key)
                continue
            # print "del: ", k, key
            if isinstance(self._lists[k], np.ndarray):
                self._lists[k] = np.delete(self._lists[k], key, axis=0)
                self._length = len(self._lists[k])
            elif isinstance(self._lists[k], (list, tuple)):
                if isinstance(key, (list, np.ndarray, tuple)):
                    indexes = sorted(list(key), reverse=True)
                    for index in indexes:
                        del self._lists[k][index]
                else:
                    del self._lists[k][key]
            else:
                del self._lists[k][key]
                # self._length = len(self._lists[k])

    def check_consistency(self):
        """

        Returns:

        """
        for key, val in self._lists.items():
            # for val in self._lists.values():
            #     print ("consistency: ", key, len(val), len(self))
            if not (len(val) == self._length):
                raise AssertionError()

    def __str__(self):
        out_str = "\n"
        for key, val in self._lists.items():
            out_str += key + " := [" + " ".join([str(el) for el in val]) + "] \n"
        return out_str

    def __len__(self):
        if hasattr(self, '_length'):
            return self._length
        else:
            return 0

    def __getitem__(self, item):
        new_dict = {}
        if isinstance(item, int):
            for key, value in self._lists.items():
                if value[item] is not None:
                    new_dict[key] = value[item]
            return SparseArrayElement(**new_dict)
        elif isinstance(item, (str, np.str, np.str_)):
            return self._lists[item]

        elif isinstance(item, (list, np.ndarray)):
            # print("key(__getitem__) len, type, item[0]: ", len(item), type(item), item[0])
            if len(item) == len(self):
                if isinstance(item[0], (np.bool_, bool)):
                    item = np.arange(len(item))[item]
        for key, value in self._lists.items():
            # print ('key: ', key, type(value))
            if isinstance(item, slice):
                new_dict[key] = value[item]
            else:
                if isinstance(value, (list, tuple)):
                    new_dict[key] = [value[i] for i in item]
                else:
                    if len(value) > 0:
                        try:
                            new_dict[key] = value[item]
                        except IndexError:
                            print('Index error:: ', key, item, value)
                    # else:
                    #     new_dict[key] = []
        # print ("new_dict: ", new_dict, self.__class__)
        return self.__class__(**new_dict)

    def keys(self):
        """

        Returns:

        """
        return self._lists.keys()

    def items(self):
        """

        Returns:

        """
        return self._lists.items()

    def __copy__(self):
        """

        Returns:

        """
        cls = self.__class__
        result = cls.__new__(cls)
        result.__init__()
        for k, v in self.__dict__.items():
            if k == '_lists':
                result.__dict__[k] = {}
                for key, val in self._lists.items():
                    if isinstance(val, SparseList):
                        result.__dict__[k][key] = val.__copy__()
                    elif isinstance(val, list):
                        result.__dict__[k][key] = val[:]
                    else:
                        result.__dict__[k][key] = np.copy(val)
            else:
                result.__dict__[k] = v
        return result

    def __add__(self, other):
        # print "__add__.new_elements"
        # assert(isinstance(other, self.__class__))
        new_array = self.__copy__()
        for key, val in other.items():
            if key not in self.keys():
                if isinstance(val, SparseList):
                    new_array._lists[key] = SparseList({}, default=other._lists[key]._default, length=len(self))
                else:
                    raise ValueError('Incompatible lists (for non-sparse lists keys must be identical (1)' + str(key))

        new_length = len(self) + len(other)
        for key, val in new_array.items():
            # print "key: ", key, val.__class__, isinstance(new_array, SparseList)
            if key in other.keys():
                if isinstance(new_array._lists[key], np.ndarray):
                    new_array._lists[key] = np.append(new_array._lists[key], other._lists[key], axis=0)
                elif isinstance(new_array._lists[key], (list, SparseList)):
                    new_array._lists[key] += other._lists[key]
                else:
                    raise ValueError("Type not implemented " + str(type(new_array._lists[key])))
            elif isinstance(val, SparseList):
                new_array._lists[key]._length = new_length  # TODO: default extends to all elements (may be undesired)
            else:
                print("non-matching key: ", key)
                raise ValueError('Incompatible lists (for non-sparse lists keys must be identical (2)')
        new_array._length += len(other)
        return new_array

    def __mul__(self, other):
        if not isinstance(other, int):
            raise ValueError('Multiplication with SparseMatrix only implemented for integers')
        new_array = self.__copy__()
        for key, value in self.items():
            new_array._lists[key] *= other

        new_array._length *= other
        return new_array

    def __rmul__(self, other):
        if isinstance(other, int):
            return self * other

    def add_tag(self, *args, **qwargs):
        for key in args:
            self._lists[key] = SparseList({}, length=len(self))

        for key, default in qwargs.items():
            self._lists[key] = SparseList({}, default=default, length=len(self))

    def remove_tag(self, *args, **qwargs):
        """

        Args:
            *args:
            **qwargs:

        Returns:

        """
        for key in args:
            del self._lists[key]
        for key, default in qwargs.items():
            del self._lists[key]
