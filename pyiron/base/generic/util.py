# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
Utility functions used in pyiron.
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


def static_isinstance(obj, obj_type):
    """
    A static implementation of isinstance() - instead of comparing an object and a class, the object is compared to a
    string, like 'pyiron.base.job.generic.GenericJob' or a list of strings.

    Args:
        obj: the object to check
        obj_type (str/list): object type as string or a list of object types as string.

    Returns:
        bool: [True/False]
    """
    if not hasattr(obj, '__mro__'):
        obj = obj.__class__
    obj_class_lst = ['.'.join([subcls.__module__, subcls.__name__]) for subcls in obj.__mro__]
    if isinstance(obj_type, list):
        return any([obj_type_element in obj_class_lst for obj_type_element in obj_type])
    elif isinstance(obj_type, str):
        return obj_type in obj_class_lst
    else:
        raise TypeError()
