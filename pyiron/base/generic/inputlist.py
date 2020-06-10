"""
Lists structure for versatile input handling.
"""

import copy
from collections.abc import Sequence, Set, Mapping, MutableMapping
import json
import numpy as np

def _normalize(key):
    if isinstance(key, str):
        if key.isdecimal():
            return int(key)
        elif '/' in key:
            return tuple(key.split('/'))

    elif isinstance(key, tuple) and len(key) == 1:
        return _normalize(key[0])

    return key



class InputList(MutableMapping):
    '''
    Mutable sequence with optional keys.

    If no argument is given, the constructor creates a new empty InputList.  If
    specified init maybe a Sequence, Set or Mapping and all recursive
    occurrences of these are also wrapped by InputList.
    >>> pl = InputList([3, 2, 1, 0])
    >>> pm = InputList({'foo': 24, 'bar': 42})

    Access can be like a normal list with integers or optionally with strings
    as keys.
    >>> pl[0]
    3
    >>> pl[2]
    1
    >>> pm['foo']
    24

    Keys do not have to be present for all elements.
    >>> pl2 = InputList([1,2])
    >>> pl2['end'] = 3
    >>> pl2
    InputList({0: 1, 1: 2, 'end': 3})

    It is also allowed to set an item one past the length of the InputList,
    this is then equivalent to appending that element.  This allows to use the
    update method also with other InputLists
    >>> pl[len(pl)] = -1
    >>> pl
    InputList([3, 2, 1, 0, -1])
    >>> pl.pop(-1)
    -1

    Where strings are used they may also be used as attributes.  Getting keys
    which clash with methods of InputList must be done with item access, but
    setting them works without overwriting the instance methods, but is not
    recommended for readability.
    >>> pm.foo
    24
    >>> pm.append = 23
    >>> pm
    InputList({'foo': 24, 'bar': 42, 'append': 23})

    Keys and indices can be tuples to traverse nested InputLists.
    >>> pn = InputList({'foo': {'bar': [4, 2]}})
    >>> pn['foo', 'bar']
    InputList([4, 2])
    >>> pn['foo', 'bar', 0]
    4

    Using keys with '/' in them is equivalent to the above after splitting the
    key.
    >>> pn['foo/bar'] == pn['foo', 'bar']
    True
    >>> pn['foo/bar/0'] == pn['foo', 'bar', 0]
    True

    To make that work strings that are only decimal digits are automatically
    converted to integers before accessing the list and keys are restricted to
    not only contain digits on initialization.
    >>> pl['0'] == pl[0]
    True
    >>> InputList({1: 42})
    Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
        File "proplist.py", line 126, in __init__
            raise ValueError(
    ValueError: keys in initializer must not be int or str of decimal digits or in correct order, is 1

    When initializing from a dict, it may not have integers or decimal strings
    as keys unless they match their position in the insertion order.  This is
    to avoid ambiguities in the final order of the InputList.
    >>> InputList({0: 'foo', 1: 'bar', 2: 42})
    InputList(['foo', 'bar', 42])
    >>> InputList({0: 'foo', 2: 42, 1: 'bar'})
    Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
        File "proplist.py", line 132, in __init__
            raise ValueError(
    ValueError: keys in initializer must not be int or str of decimal digits or in correct order, is 2


    Using keys is completely optional, InputList can always be treated as a
    list, with the exception that `iter()` iterates of the keys and indices.
    This is to correctly implement the MutableMapping protocol, to convert to a
    normal list and discard the keys use `values()`.
    >>> pm[0]
    24
    >>> pn['0/0/1']
    2
    >>> list(pl)
    [0, 1, 2, 3]
    >>> list(pl.values())
    [3, 2, 1, 0]
    >>> list(pl.keys())
    [0, 1, 2, 3]
    '''

    __version__ = "0.1.0"

    def __new__(cls, *args, **kwargs):

        instance = super().__new__(cls)
        # setting these immediately after object creation ensures that they are
        # always defined and attribute access works even before __init__ is
        # called.  This is relevant on deepcopy & pickling.
        object.__setattr__(instance, '_store', [])
        object.__setattr__(instance, '_indices', {})
        object.__setattr__(instance, 'table_name', None)

        return instance

    def __init__(self, init = None, table_name = None):

        self.table_name = table_name

        if init == None: return

        if isinstance(init, (Sequence, Set)):
            for v in init:
                self.append(self._wrap_val(v))

        elif isinstance(init, Mapping):
            for i, (k, v) in enumerate(init.items()):
                k = _normalize(k)
                v = self._wrap_val(v)
                if isinstance(k, int):
                    if k == i:
                        self.append(v)
                    else:
                        raise ValueError(
                            'keys in initializer must not be int or str of '
                            'decimal digits or in correct order, '
                            'is {!r}'.format(k))
                else:
                    self[k] = v
        else:
            ValueError('init must be Sequence, Set or Mapping')

    def __len__(self):
        return len(self._store)

    def __iter__(self):

        reverse_indices = {i: k for k, i in self._indices.items()}

        for i in range(len(self)):
            yield reverse_indices.get(i, i)

    def __getitem__(self, key):

        key = _normalize(key)

        if isinstance(key, tuple):
            return self[key[0]][key[1:]]

        elif isinstance(key, int):

            try:
                return self._store[key]
            except IndexError:
                raise IndexError('list index out of range') from None
        elif isinstance(key, str):

            try:
                return self._store[self._indices[key]]
            except KeyError:
                raise KeyError(repr(key)) from None
        else:
            raise ValueError(
                    '{} is not a valid key, must be str or int'.format(key)
            )

    def __setitem__(self, key, val):

        key = _normalize(key)

        if isinstance(key, tuple):
            self[key[0]][key[1:]] = val
        elif isinstance(key, int):
            if key < len(self):
                self._store[key] = val
            elif key == len(self):
                self.append(val)
            else:
                raise IndexError("index out of range")
        elif isinstance(key, str):
            if key not in self._indices:
                self._indices[key] = len(self._store)
                self._store.append(val)
            else:
                self._store[self._indices[key]] = val
        else:
            raise ValueError(
                    '{} is not a valid key, must be str or int'.format(key)
            )

    def __delitem__(self, key):

        key = _normalize(key)

        if isinstance(key, tuple):
            del self[key[0]][key[1:]]
        elif isinstance(key, (str, int)):
            if isinstance(key, str):
                idx = self._indices[key]
                del self._indices[key]
            else:
                idx = key

            del self._store[idx]

            for k, i in self._indices.items():
                if i > idx:
                    self._indices[k] = i - 1
        else:
            raise ValueError(
                    '{} is not a valid key, must be str or int'.format(key)
            )

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name) from None

    def __setattr__(self, name, val):
        if name in self.__dict__:
            object.__setattr__(self, name, val)
        else:
            self[name] = val

    def __delattr__(self, name):
        if name in self.__dict__:
            object.__delattr__(self, name)
        else:
            del self[name]

    def __array__(self):
        """Return bare list of values to play nice with numpy."""
        return np.array(self._store)

    def __dir__(self):
        return set(super().__dir__() + list(self._indices.keys()))

    def __repr__(self):
        name = self.__class__.__name__
        if self.has_keys():
            return name + '({' + \
                    ', '.join('{!r}: {!r}'.format(k, v) for k, v in self.items()) + \
            '})'
        else:
            return name + '([' + \
                    ', '.join('{!r}'.format(v) for v in self._store) + \
            '])'

    @classmethod
    def _wrap_val(cls, val):
        if isinstance(val, (tuple, list, dict)):
            return cls(val)
        else:
            return val

    def to_builtin(self):
        '''
        Convert the list back to builtin dict's and list's recursively.
        '''

        if self.has_keys():
            dd = dict(self)
            for k, v in dd.items():
                if isinstance(v, InputList):
                    dd[k] = v.to_builtin()

            return dd
        else:
            return list(v.to_builtin() if isinstance(v, InputList) else v
                            for v in self.values())

    # allows 'nice' displays in jupyter notebooks
    _repr_json_ = to_builtin

    def append(self, val):
        '''
        Add new value to the list without a key.
        '''
        self._store.append(self._wrap_val(val))

    def extend(self, vals):
        '''
        Append all items of vals to this InputList.
        '''

        for v in vals:
            self.append(v)

    def insert(self, index, val, key = None):
        '''
        Add a new value to the list at the specified position, with an optional
        key.  If the key is already in the list it will be updated to point to
        the new value at the new index.
        '''
        if key != None:
            for k, i in self._indices.items():
                if i >= index:
                    self._indices[k] = i + 1
            self._indices[key] = index

        self._store.insert(index, val)

    def mark(self, index, key):
        '''
        Add a key to an existing item at index.  If a key already exists, it is
        overwritten.
        >>> pl = InputList([42])
        >>> pl.mark(0, 'head')
        >>> pl.head == 42
        True
        '''
        if index >= len(self):
            raise IndexError('list index out of range')

        reverse_indices = {i: k for k, i in self._indices.items()}
        if index in reverse_indices:
            del self._indices[reverse_indices[index]]

        self._indices[key] = index

    def clear(self):
        """
        Remove all items from InputList.
        """
        self._store.clear()
        self._indices.clear()

    def create_group(self, name):
        '''
        Add a new empty sublist under the given key.

        Returns:
            the newly created sublist

        >>> pl = InputList({})
        >>> pl.create_group('group_name')
        InputList([])
        >>> list(pl.group_name)
        []
        '''
        self[name] = self.__class__()
        return self[name]

    def has_keys(self):
        """
        Returns True if there is a key set for at least one item of the list.
        """
        return bool(self._indices)

    def __copy__(self):
        # by default copy.copy will use the same objects for _store and
        # _indices, which would cause the copied and the copiee to have the
        # same underlying data storage, so instead we have to do a shallow copy
        # of those manually
        copiee = type(self)()
        copiee._store  = copy.copy(self._store)
        copiee._indices = copy.copy(self._indices)
        copiee.table_name = self.table_name
        return copiee

    def copy(self):
        """
        Returns deep copy of it self.

        >>> pl = InputList([[1,2,3]])
        >>> pl.copy() == pl
        True
        >>> pl.copy() is pl
        False
        >>> all(a is not b for a, b in zip(pl.copy().values(), pl.values()))
        True
        """
        return copy.deepcopy(self)

    def to_hdf(self, hdf, group_name=None):
        """
        Store the InputList in an HDF5 file.  If ``group_name`` or
        *self.table_name* are not `None`, create a sub group in hdf prior to
        writing if not save directly to hdf.  group_name overrides
        self.table_name if both are not None.

        Args:
            hdf (ProjectHDFio): HDF5 group object
            group_name (str, optional): HDF5 subgroup name, overrides
                self.table_name
        """

        group_name = group_name or self.table_name
        if group_name:
            hdf = hdf.create_group(group_name)

        self._type_to_hdf(hdf)
        hdf["data"] = json.dumps(self.to_builtin())

    def _type_to_hdf(self, hdf):
        """
        Internal helper function to save type and version in hdf root

        Args:
            hdf (ProjectHDFio): HDF5 group object
        """
        hdf["NAME"] = self.__class__.__name__
        hdf["TYPE"] = str(type(self))
        hdf["VERSION"] = self.__version__
        hdf["OBJECT"] = "InputList"

    def from_hdf(self, hdf, group_name=None):
        """
        Restore the InputList from an HDF5 file.  If group_name or
        self.table_name are not None, open a sub group in hdf prior to reading
        if not read directly from hdf.  group_name overrides self.table_name if
        both are not None.

        Args:
            hdf (ProjectHDFio): HDF5 group object
            group_name (str, optional): HDF5 subgroup name, overrides
                self.table_name
        """

        group_name = group_name or self.table_name
        if group_name:
            with hdf.open(group_name) as hdf_group:
                data = hdf_group["data"]
        else:
            data = hdf["data"]

        self.clear()
        self.update(json.loads(data))
