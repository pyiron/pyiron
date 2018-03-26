# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
ClassTemplate for pyiron Configurations
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class GenericConfig(object):
    """
    The GenericConfig object defines the interface every configuration parser should implement.
    """

    @property
    def pyiron_envs(self):
        """
        Get the list of pyiron environments, where each environment consists of a dictionary with the entries:
                'system': str,
                'type': str [Postgres/ SQLite],
                'database': str,
                'user': str (optional),
                'password': str (optional),
                'table_name': str,
                'host': str (optional),
                'top_level_dirs': str,
                'ssh_host': ssh host object (optional)

        Returns:
            list: list of pyiron environments
        """
        raise ValueError("Pyiron Environment property(pyiron_env) has to be set in derived Configuration.")

    @property
    def login_user(self):
        """
        Get the username of the current user

        Returns:
            str: username
        """
        raise ValueError("Login user property(login_user) has to be set in derived Configuration.")

    @property
    def path_bin(self):
        """
        Get the path where the Hamilton executables are located

        Returns:
            str: path
        """
        raise ValueError("Path bin property(path_bin) has to be set in derived Configuration.")

    @property
    def path_potentials(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            str: path
        """
        raise ValueError("Path potentials property(path_potentials) has to be set in derived Configuration.")

    @property
    def path_pyiron(self):
        """
        Get the path where the PyIron sourcecode is installed

        Returns:
            str: path
        """
        raise ValueError("Pyiron Sourcecode property(pyth_pyiron) has to be set in derived Configuration.")

    @property
    def resource_paths(self):
        """
        Get the paths where the PyIron sourcecode is installed

        Returns:
            list: path
        """
        raise ValueError("Pyiron Sourcecode property(pyth_pyiron) has to be set in derived Configuration.")