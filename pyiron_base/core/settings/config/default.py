# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import getpass
import inspect
import os
from pyiron_base.core.settings.config.template import GenericConfig

"""
Class for the default pyiron configuration - if no .pyiron file exists this one is used. 
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class ConfigDefault(GenericConfig):
    """
    The default configuration - if no .pyiron file exists this one is used.

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
    def __init__(self):
        module_path = os.path.dirname(os.path.abspath(inspect.getsourcefile(ConfigDefault)))
        if os.name == 'nt':
            self.pyiron_home = (self._dir_check(os.path.expanduser('~') + '/PyIron_data')).replace('\\', '/')
            self.pyiron_code = (module_path.split("\\pyiron")[0]).replace('\\', '/')
        else:
            self.pyiron_home = self._dir_check(os.path.expanduser('~') + '/PyIron_data')
            pyiron_path_lst = module_path.split("/pyiron_base")
            if len(pyiron_path_lst) == 2:
                self.pyiron_code = module_path.split("/pyiron_base")[0]
            else:
                self.pyiron_code = module_path.split("/pyiron_base")[0]
                for pyiron_str in pyiron_path_lst[1:-1]:
                    self.pyiron_code += "/pyiron"
        top_level_name_1 = self.pyiron_code + '/examples/'
        top_level_name_2 = self.pyiron_code + '/tests/'
        self._pyiron_envs = {'system': 'default',
                             'type': 'SQLite',
                             'file': self.pyiron_home + '/pyiron.db',
                             'table_name': 'jobs_' + getpass.getuser(),
                             'user': getpass.getuser(),
                             'top_level_dirs': {top_level_name_1: top_level_name_1, top_level_name_2: top_level_name_2},
                             'resource_paths': [self._dir_check(self.pyiron_code + '/static')]}

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

    @property
    def path_bin(self):
        """
        Get the path where the Hamilton executables are located

        Returns:
            str: path
        """
        return self._dir_check(self.pyiron_home + '/bin')

    @property
    def path_potentials(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            str: path
        """
        return self._dir_check(self.pyiron_code + '/static/potentials')

    @property
    def path_pyiron(self):
        """
        Get the path where the PyIron sourcecode is installed

        Returns:
            str: path
        """
        return self.pyiron_code

    @staticmethod
    def _dir_check(directory):
        """
        Validate that the directory exists and if it does not create a new one.

        Args:
            directory (str): relative path to the new directory

        Returns:
            str: relative path to the new directory
        """
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
