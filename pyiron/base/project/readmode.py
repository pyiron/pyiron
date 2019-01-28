# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
This class is currently not used !
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"

class ReadMode(object):
    def __init__(self, mode='object'):
        self._inspect = False
        self._object = False
        self.mode = mode

    @property
    def inspect(self):
        return self._inspect

    @inspect.setter
    def inspect(self, boolean):
        if isinstance(boolean, bool):
            if boolean:
                self.mode = 'inspect'
            else:
                self.mode = 'object'
        else:
            raise TypeError('To set the ReadMode use boolean.')

    @property
    def object(self):
        return self._object

    @object.setter
    def object(self, boolean):
        if isinstance(boolean, bool):
            if boolean:
                self.mode = 'object'
            else:
                self.mode = 'inspect'
        else:
            raise TypeError('To set the ReadMode use boolean.')

    @property
    def mode(self):
        if self._inspect:
            return 'inspect'
        elif self._object:
            return 'object'
        else:
            raise ValueError('The ReadMode is currently in an undefined state.')

    @mode.setter
    def mode(self, new_mode):
        if isinstance(new_mode, str):
            if new_mode == 'object':
                self._object = True
                self._inspect = False
            elif new_mode == 'inspect':
                self._object = False
                self._inspect = True
            else:
                raise ValueError('The ReadMode can either be "object" or "inspect".')
        else:
            raise TypeError('The ReadMode mode can only be set using strings.')
