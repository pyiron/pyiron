# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.


import os
import warnings
from six import with_metaclass
from pyironbase.core.settings.logger import setup_logger
from pyironbase.core.settings.config.default import ConfigDefault
from pyironbase.core.settings.config.parser import ConfigFile
from pyironbase.core.settings.config.template import GenericConfig
from pyironbase.core.settings.database import DatabaseAccess

"""
The settings file provides the attributes of the configuration as properties.
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Singleton(type):
    """
    Implemented with suggestions from

    http://stackoverflow.com/questions/6760685/creating-a-singleton-in-python

    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Settings(with_metaclass(Singleton)):
    """
    The settings object can either search for an configuration file and use the default configuration only when no
    other configuration file is found, or it can be forced to use the default configuration file.

    Args:
        config (GenericConfig object instance): Provide a GenericConfig Instance to force pyiron to use that specific
                                                configuration (optional).
    """
    def __init__(self, config=None):
        self.db_name = None
        self.db_connection_string = None
        self.db_connection_table = None
        self.db_connection_name = None
        self.db_dict = {}
        self.ssh_dict = {}
        self.db_translate_dict = {}
        self.top_path_dict = {}

        # Load config file if it exists or otherwise load default configuration
        config_file = os.path.expanduser(os.path.join("~", ".pyiron"))
        if config:
            if isinstance(config, GenericConfig):
                self._config = config
            else:
                raise TypeError('The config parameter has to be an object instance dereived from GenericConfig.')
        elif os.path.isfile(config_file):
            self._config = ConfigFile(config_file=config_file)
        else:
            self._config = ConfigDefault()

        self._viewer_conncetion_string = None
        self._viewer_connection_table = None
        self._env_load()
        self.logger = setup_logger()

    @property
    def login_user(self):
        """
        Get the username of the current user

        Returns:
            str: username
        """
        if 'user' in self._config.pyiron_envs.keys():
            return self._config.pyiron_envs['user']
        else:
            warnings.warn("Please update your .pyiron config file!", DeprecationWarning)
            return self._config.login_user

    # local paths
    @property
    def path_pyiron(self):
        """
        Get the path where the PyIron sourcecode is installed

        Returns:
            str: path
        """
        warnings.warn("Please update your PYTHONPATH.", DeprecationWarning)
        return self._config.path_pyiron

    @property
    def path_bin(self):
        """
        Get the path where the Hamilton executables are located

        Returns:
            str: path
        """
        warnings.warn("Please update your .pyiron config file!", DeprecationWarning)
        return self._config.path_bin

    @property
    def path_potentials(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            str: path
        """
        warnings.warn("Please update your .pyiron config file!", DeprecationWarning)
        return self._config.path_potentials

    @property
    def resource_paths(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            list: path of paths
        """
        if 'resource_paths' in self._config.pyiron_envs.keys():
            return self._config.pyiron_envs['resource_paths']
        else:
            warnings.warn("Please update your .pyiron config file!", DeprecationWarning)
            return self._config.resource_paths

    def __del__(self):
        """
        Close database connection
        """
        self.close_connection()

    def open_connection(self):
        """
        Internal function to open the connection to the database. Only after this function is called the database is
        accessable.
        """
        if not self.db_dict:
            self.db_dict[self.db_connection_name] = DatabaseAccess(self.db_connection_string, self.db_connection_table)

    def switch_to_viewer_mode(self):
        """
        Switch from user mode to viewer mode - if viewer_mode is enable pyiron has read only access to the database.
        """
        if self._viewer_conncetion_string is not None and self._viewer_connection_table is not None:
            if not self.db_dict[self.db_connection_name].viewer_mode:
                self.close_connection()
                database = DatabaseAccess(self._viewer_conncetion_string, self._viewer_connection_table)
                database.viewer_mode = True
                self.db_dict[self.db_connection_name] = database
            else:
                print('Database is already in viewer mode!')
        else:
            print('Viewer Mode is not available on this pyiron installation.')

    def switch_to_user_mode(self):
        """
        Switch from viewer mode to user mode - if viewer_mode is enable pyiron has read only access to the database.
        """
        if self._viewer_conncetion_string is not None and self._viewer_connection_table is not None:
            if self.db_dict[self.db_connection_name]._viewer_mode:
                self.close_connection()
                database = DatabaseAccess(self.db_connection_string, self.db_connection_table)
                database._viewer_mode = False
                self.db_dict[self.db_connection_name] = database
            else:
                print('Database is already in user mode!')
        else:
            print('Viewer Mode is not available on this pyiron installation.')

    def close_connection(self):
        """
        Internal function to close the connection to the database.
        """
        if self.db_dict:
            database = self.db_dict[self.db_connection_name]
            database.conn.close()

    def top_path(self, full_path):
        """
        Validated that the full_path is a sub directory of one of the pyrion environments loaded.

        Args:
            full_path (str): path

        Returns:
            str: path
        """
        for path in self.top_path_dict:
            if path in full_path:
                return path
        if 'ConfigDefault' in str(type(self._config)):
            raise ValueError('no config file found - using default config '
                             '- but the current path {0} is not included.'.format(full_path))
        else:
            raise ValueError("Path '{0}' is not included in top_level_dirs: {1} ".format(full_path,
                                                                                         str(self.top_path_dict)))

    # Private functions
    def _env_load(self):
        """
        Private function to convert the PyIron environment variables to and SQLalchemy connection string.
        """

        pyiron_env = self._config.pyiron_envs
        if pyiron_env['type'] == 'SQLite':
            sql_con_str = 'sqlite:///' + pyiron_env['file']
            sql_db_table = pyiron_env['table_name']
        elif pyiron_env['type'] == 'Postgres':
            sql_con_str = 'postgresql://' + pyiron_env['user'] + ':' + pyiron_env['password'] + '@' \
                          + pyiron_env['host'] + '/' + pyiron_env['database']
            sql_db_table = pyiron_env['table_name']
            if 'vieweruser' in pyiron_env.keys():
                self._viewer_conncetion_string = 'postgresql://' + pyiron_env['vieweruser'] + ':' + \
                                                 pyiron_env['viewerpassword'] + '@' + pyiron_env['host'] + '/' + \
                                                 pyiron_env['database']
                self._viewer_connection_table = pyiron_env['viewertable']
        elif pyiron_env['type'] == 'SQLalchemy':
            sql_con_str = pyiron_env['connection_string']
            sql_db_table = pyiron_env['table_name']
        elif pyiron_env['type'] == 'MySQL':
            sql_con_str = 'mysql+pymysql://' + pyiron_env['user'] + ':' + pyiron_env['password'] + '@' \
                          + pyiron_env['host'] + '/' + pyiron_env['database']
            sql_db_table = pyiron_env['table_name']
        else:
            raise TypeError('Currently only SQLite and Postgres databases are supported.')
        self.db_connection_name = pyiron_env['system']
        self.db_connection_string = sql_con_str
        self.db_connection_table = sql_db_table
        self.db_translate_dict[pyiron_env['system']] = pyiron_env['top_level_dirs']
        for local_path in pyiron_env['top_level_dirs'].values():
            local_path = local_path.replace("\\", "/")
            self.top_path_dict[local_path] = pyiron_env['system']
        self.db_name = pyiron_env['system']
