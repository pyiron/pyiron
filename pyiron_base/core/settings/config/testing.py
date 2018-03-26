# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron_base.core.settings.config.template import GenericConfig

"""
Class to represent the pyiron configuration for automated testing 
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class ConfigTesting(GenericConfig):
    """
    A test configuration for the automated unit tests, the configuration entries can be provided as
    input variables, no .pyiron configuration file is required to execute the unit tests.

    Args:
        login_user (str): Username of the current pyiron user
        sql_lite_database (str): Link to the SQLite database which should be used for testing
        path_bin (str): Path to the executables
        path_potentials (str): Path to the external calculation resources like emperical potentials
        path_project (str): Root path for the current pyiron user
        path_sourcecode (str): Path to the pyiron source code

    Attributes:

        .. attribute:: pyiron_envs

            pyiron environments containing the top level directories and the corresponding database tables

        .. attribute:: login_user

            Username of the current pyiron user

        .. attribute:: path_bin

            Path to the executables

        .. attribute:: path_potentials

            Path to the external calculation resources like emperical potentials

        .. attribute:: path_pyiron

            Path to the pyiron source code
    """
    def __init__(self, login_user='pyiron', sql_lite_database='./pyiron.db', path_bin='.', path_potentials='.',
                 path_project='.', path_sourcecode='.', resource_paths=None):
        if path_project[-1] != '/':
            path_project += '/'
        self._pyiron_envs = {'system': 'default',
                             'type': 'SQLite',
                             'file': sql_lite_database,
                             'table_name': 'jobs_' + login_user,
                             'top_level_dirs': {path_project: path_project},
                             'user': login_user}
        if resource_paths is None:
            self._pyiron_envs['resource_paths'] = [path_potentials, path_bin]
        else:
            self._pyiron_envs['resource_paths'] = resource_paths
        self._path_bin = path_bin
        self._path_potentials = path_potentials
        self._path_code = path_sourcecode

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
            dict: pyiron environment
        """
        return self._pyiron_envs

    # Backward compatibility
    @property
    def path_bin(self):
        """
        Get the path where the Hamilton executables are located

        Returns:
            str: path
        """
        return self._path_bin

    @property
    def path_potentials(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            str: path
        """
        return self._path_potentials

    @property
    def path_pyiron(self):
        """
        Get the path where the PyIron sourcecode is installed

        Returns:
            str: path
        """
        return self._path_code
