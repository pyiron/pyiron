# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
ShellOption class as part of the queueing system interface
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class ShellOption(object):
    """
    The ShellOption allows to match the parameters of the queing system with the requirements of pyiron.

    Args:
        name (str): name of the ShellOption
        key (str): command line option for the specific queing system - (optional)
        value (str): value the command line option of the specific queing system should be set to - (optional)
        value_prefix (str): a prefix for the command line option for the specific queing system - (optional)
        active (bool): if a tag is not supported by a certain queing system set it to False - [True/False]
    """
    def __init__(self, name, key=None, value=None, value_prefix=None, active=True):
        self.__name__ = name
        self._key = None
        self._value = None
        self._value_prefix = None
        self._active = False

        self.key = key
        self.value = value
        self.value_prefix = value_prefix
        self.active = active

    @property
    def active(self):
        """
        Status of the ShellOption, active ?

        Returns:
            bool: [True/False]
        """
        return self._active

    @active.setter
    def active(self, new_active):
        """
        Set the stauts of the ShellOption to active = True or active = False.

        Args:
            new_active (bool): [True/False]
        """
        if not isinstance(new_active, bool):
            raise TypeError('The ShellOption active property has to be bool.')
        self._active = new_active

    @property
    def key(self):
        """
        Command line option for the specific queing system

        Returns:
            str:

        """
        return self._key

    @key.setter
    def key(self, new_key):
        """
        Set the command line option for the specific queing system

        Args:
            new_key (str):
        """
        if new_key:
            if not isinstance(new_key, str):
                raise TypeError('The ShellOption key has to be str.')
            self._key = new_key

    @property
    def value(self):
        """
        Value of the command line option for the specific queing system

        Returns:
            str:
        """
        return self._value

    @value.setter
    def value(self, new_value):
        """
        Set the value of the command line option for the specific queing system

        Args:
            new_value (str):
        """
        if new_value:
            if not isinstance(new_value, str):
                raise TypeError('The ShellOption key has to be str.')
            self._value = new_value

    @property
    def value_prefix(self):
        """
        Value prefix of the command line option for the specific queing system

        Returns:
            str:
        """
        return self._value_prefix

    @value_prefix.setter
    def value_prefix(self, new_value_prefix):
        """
        Set the value prefix of the command line option for the specific queing system

        Args:
            new_value_prefix (str):
        """
        if new_value_prefix:
            if not isinstance(new_value_prefix, str):
                raise TypeError('The ShellOption key has to be str.')
            self._value_prefix = new_value_prefix

    @property
    def value_only(self):
        """
        Verify if the Shelloption is active and only contains a value but no key.

        Returns:
            bool: [True/False]
        """
        if self.active and not self.key and self.value:
            return True
        else:
            return False

    @property
    def key_only(self):
        """
        Verify if the Shelloption is active and only contains a key but no value.

        Returns:
            bool: [True/False]
        """
        if self.active and self.key and not self.value:
            return True
        else:
            return False

    def output_as_list(self):
        """
        Output the shell option as list, so it can be used in combination with a python sub process.

        Returns:
            list: list of strings
        """
        if self.value_only:
            if self.value_prefix:
                return [self.value_prefix+self.value]
            else:
                return [self.value]
        elif self.key_only:
            return [self.key]
        elif self.active:
            if self.value_prefix:
                return [self.key, self.value_prefix+self.value]
            else:
                return [self.key, self.value]
        else:
            return []
