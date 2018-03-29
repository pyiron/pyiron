# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import sys
from pyiron_base.core.settings.config.template import GenericConfig

"""
Class for parsing the .pyiron configuration file 
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

# Depending on the Python Version - a specific version of the config parser is required.
if sys.version_info.major == 2:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser


class ConfigFile(GenericConfig):
    """
    Configuration which is loaded from the .pyiron file. This is the most common usecase.

    Args:
        config_file (str): path

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
    def __init__(self, config_file):
        if sys.version_info.major == 2:
            config = ConfigParser()
        else:
            config = ConfigParser(inline_comment_prefixes=(';',))
        if not os.path.isfile(config_file):
            raise ValueError("Configuration file missing", os.path.abspath(os.path.curdir))
        config.read(config_file)
        self.parser = config

    @property
    def login_user(self):
        """
        Get the username of the current user

        Returns:
            str: username
        """
        if self.parser.has_option('DEFAULT', 'LOGIN_USER'):
            return self.parser.get('DEFAULT', 'LOGIN_USER')
        elif self.parser.has_option('DEFAULT', 'USER'):
            return self.parser.get('DEFAULT', 'USER')
        else: 
            return 'pyiron'

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
        for section_name in self.parser.sections():
            if "database_" in section_name:
                return self._env_config(section=section_name)
        return self._env_config(section='DEFAULT')

    @property
    def path_bin(self):
        """
        Get the path where the Hamilton executables are located

        Returns:
            str: path
        """
        return self.parser.get('DEFAULT', 'LOCAL_PATH_BIN')

    @property
    def path_potentials(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            str: path
        """
        return self.parser.get('DEFAULT', 'LOCAL_PATH_POTS')

    @property
    def resource_paths(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            str: path
        """
        if self.parser.has_option('DEFAULT', 'RESOURCE_PATHS'):
            return [os.path.expanduser(rpath) 
                    for rpath in self.parser.get('DEFAULT', 'RESOURCE_PATHS').split(",")]
        else:
            return None

    @property
    def path_pyiron(self):
        """
        Get the path where the PyIron sourcecode is installed

        Returns:
            str: path
        """
        file_path = os.path.realpath(__file__).replace('\\', '/')
        return '/'.join(file_path.split('/')[:-5])

    # private functions
    def _top_level_dirs(self, section):
        """
        Get the local top level directories which are connected to the specific section.

        Args:
            section (str): section in the config file

        Returns:
            dict: dictionary with the local files and their remote counterparts.
        """
        top_dir_dict = {}
        top_dir_lst = []
        if self.parser.has_option(section, "TOP_LEVEL_DIRS"):
            top_dir_lst += [c.strip() for c in self.parser.get(section, "TOP_LEVEL_DIRS").split(",")]
        if self.parser.has_option(section, "PROJECT_PATHS"):
            top_dir_lst += [c.strip() for c in self.parser.get(section, "PROJECT_PATHS").split(",")]
        for top_dir in top_dir_lst:
            top_dir = [d.strip() for d in top_dir.split("@@")]
            if len(top_dir) == 2:
                local_path, db_path = top_dir
            else:
                local_path, db_path = top_dir[0], top_dir[0]
            if local_path[-1] != '/':
                local_path += '/'
            if db_path[-1] != '/':
                db_path += '/'
            top_dir_dict[os.path.expanduser(db_path).replace('\\', '/')] = os.path.expanduser(local_path).replace('\\', '/')
        return top_dir_dict

    def _env_config(self, section):
        """
        Read section in configuration file and return a dictionary with the corresponding parameters.

        Args:
            section(str): section in the config file as dictionary

        Returns:
            dict: dictionary with the environment configuration
        """
        if self.parser.has_option(section, "TYPE"):
            dbtype = self.parser.get(section, "TYPE")
        else:
            dbtype = 'SQLite'
        db_dict = {'system': section[9:],
                   'top_level_dirs': self._top_level_dirs(section=section),
                   'resource_paths': [os.path.expanduser(c.strip()).replace('\\', '/') 
                                      for c in self.parser.get(section, 'RESOURCE_PATHS').split(",")]}
        if self.parser.has_option(section, "USER"):
            db_dict['user'] = self.parser.get(section, "USER")
        else:
            db_dict['user'] = 'pyiron'
        if dbtype == 'Postgres':
            db_dict['type'] = dbtype
            db_dict['database'] = self.parser.get(section, "NAME")
            db_dict['password'] = self.parser.get(section, "PASSWD")
            db_dict['table_name'] = self.parser.get(section, "JOB_TABLE")
            db_dict['host'] = self.parser.get(section, "HOST")
            if self.parser.has_option(section, "VIEWERUSER") & self.parser.has_option(section, "VIEWERPASSWD") \
                    & self.parser.has_option(section, "VIEWER_TABLE"):
                db_dict['vieweruser'] = self.parser.get(section, "VIEWERUSER")
                db_dict['viewerpassword'] = self.parser.get(section, "VIEWERPASSWD")
                db_dict['viewertable'] = self.parser.get(section, "VIEWER_TABLE")
        elif dbtype == 'SQLite':
            db_dict['type'] = dbtype
            if self.parser.has_option(section, "FILE"):
                db_dict['file'] = os.path.expanduser(self.parser.get(section, "FILE")).replace('\\', '/')
            if self.parser.has_option(section, "DATABASE_FILE"):
                db_dict['file'] = os.path.expanduser(self.parser.get(section, "DATABASE_FILE")).replace('\\', '/')   
            else:
                db_dict['file'] = os.path.join(db_dict['resource_paths'][0], 'sqlite.db').replace('\\', '/')
            if self.parser.has_option(section, "JOB_TABLE"):
                db_dict['table_name'] = self.parser.get(section, "JOB_TABLE")
            else:
                db_dict['table_name'] = 'jobs_pyiron'
        elif dbtype == 'SQLalchemy':
            db_dict['type'] = dbtype
            db_dict['connection_string'] = self.parser.get(section, "CONNECTION")
            db_dict['table_name'] = self.parser.get(section, "JOB_TABLE")
        elif dbtype == 'MySQL':
            db_dict['type'] = dbtype
            db_dict['database'] = self.parser.get(section, "NAME")
            db_dict['user'] = self.parser.get(section, "USER")
            db_dict['password'] = self.parser.get(section, "PASSWD")
            db_dict['table_name'] = self.parser.get(section, "JOB_TABLE")
            db_dict['host'] = self.parser.get(section, "HOST")
        else:
            raise ValueError('Database configuration unreadable!')
        return db_dict
