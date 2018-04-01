# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.


import os
import warnings
from six import with_metaclass
import sys
from pathlib2 import Path
from pyiron_base.core.settings.logger import setup_logger
from pyiron_base.core.settings.database import DatabaseAccess

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

# Depending on the Python Version - a specific version of the config parser is required.
if sys.version_info.major == 2:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser


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
        config (dict): Provide a dict with the configuration.
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
            if isinstance(config, dict):
                # Setup test configuration
                if 'table_name' not in config.keys():
                    config['table_name'] = 'jobs_pyiron'
                if 'system' not in config.keys():
                    config['system'] = 'test'
                if 'type' not in config.keys():
                    config['type'] = 'SQLite'
                if 'file' not in config.keys():
                    config['file'] = 'genericsettings.db'
                if 'user' not in config.keys():
                    config['user'] = 'pyiron'
                if 'top_level_dirs' not in config.keys():
                    config['top_level_dirs'] = [os.path.abspath(os.getcwd())]
                if isinstance(config['top_level_dirs'], str):
                    config['top_level_dirs'] = {config['top_level_dirs']: config['top_level_dirs']}
                elif not isinstance(config['top_level_dirs'], dict):
                    raise TypeError('The project paths should be a dict of path pairs!')
                if 'resource_paths' not in config.keys():
                    config['resource_paths'] = [os.path.abspath(os.getcwd())]
                if isinstance(config['resource_paths'], str):
                    config['resource_paths'] = [config['resource_paths']]
                elif not isinstance(config['resource_paths'], list):
                    raise TypeError('The resource paths should be a list of strings!')
                self._config = config
            else:
                raise TypeError('The config parameter has to be an object instance dereived from GenericConfig.')
        else:
            if not os.path.isfile(config_file):
                # Write default config file
                with open(config_file, 'w') as cf:
                    cf.writelines(['[DEFAULT]\n',
                                   'PROJECT_PATHS = ~/pyiron/projects\n',
                                   'RESOURCE_PATHS = ~/pyiron/resources\n'])
                project_path = self.convert_path('~/pyiron/projects')
                if not os.path.exists(project_path):
                    os.makedirs(project_path)
                resource_path = self.convert_path('~/pyiron/resources')
                if not os.path.exists(os.path.expanduser(resource_path)):
                    os.makedirs(resource_path)
            self._config = self._env_config(config_file)

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
        return self._config['user']

    @property
    def resource_paths(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            list: path of paths
        """
        return self._config['resource_paths']

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
        raise ValueError('the current path {0} is not included in the .pyiron configuration.'.format(full_path))

    # Private functions
    def _env_load(self):
        """
        Private function to convert the PyIron environment variables to and SQLalchemy connection string.
        """

        pyiron_env = self._config
        if pyiron_env['type'] == 'SQLite':
            # SQLite is raising ugly error messages when the database directory does not exist.
            if os.path.dirname(pyiron_env['file']) != '' and not os.path.exists(os.path.dirname(pyiron_env['file'])):
                os.makedirs(os.path.dirname(pyiron_env['file']))
            sql_con_str = 'sqlite://' + pyiron_env['file']
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

    # private functions
    def _env_config(self, config_file):
        """
        Read section in configuration file and return a dictionary with the corresponding parameters.

        Args:
            config_file(str): confi file to parse

        Returns:
            dict: dictionary with the environment configuration
        """
        if sys.version_info.major == 2:
            parser = ConfigParser()
        else:
            parser = ConfigParser(inline_comment_prefixes=(';',))
        if not os.path.isfile(config_file):
            raise ValueError("Configuration file missing", os.path.abspath(os.path.curdir))
        parser.read(config_file)
        if len(parser.sections()) > 0:
            section = parser.sections()[0]
        else:
            section = 'DEFAULT'
        if parser.has_option(section, "TYPE"):
            dbtype = parser.get(section, "TYPE")
        else:
            dbtype = 'SQLite'
        top_level_dirs = {}
        project_path_lst = []
        if parser.has_option(section, "PROJECT_PATHS"):
            project_path_lst = [self.convert_path(c.strip()) for c in parser.get(section, "PROJECT_PATHS").split(",")]
        elif parser.has_option(section, "TOP_LEVEL_DIRS"):
            project_path_lst = [self.convert_path(c.strip()) for c in parser.get(section, "TOP_LEVEL_DIRS").split(",")]
        else:
            ValueError('No project path identified!')
        for top_dir in project_path_lst:
            top_dir = [d.strip() for d in top_dir.split("@@")]
            if len(top_dir) == 2:
                local_path, db_path = top_dir
            else:
                local_path, db_path = top_dir[0], top_dir[0]
            if local_path[-1] != '/':
                local_path += '/'
            if db_path[-1] != '/':
                db_path += '/'
            top_level_dirs[db_path] = local_path
        resource_paths = [self.convert_path(c.strip())
                          for c in parser.get(section, 'RESOURCE_PATHS').split(",")]
        if dbtype == 'Postgres':
            db_dict = {'system': section,
                       'type': dbtype,
                       'database': parser.get(section, "NAME"),
                       'user': parser.get(section, "USER"),
                       'password': parser.get(section, "PASSWD"),
                       'table_name': parser.get(section, "JOB_TABLE"),
                       'host': parser.get(section, "HOST")}
            if parser.has_option(section, "VIEWERUSER") & parser.has_option(section, "VIEWERPASSWD") \
                    & parser.has_option(section, "VIEWER_TABLE"):
                db_dict['vieweruser'] = parser.get(section, "VIEWERUSER")
                db_dict['viewerpassword'] = parser.get(section, "VIEWERPASSWD")
                db_dict['viewertable'] = parser.get(section, "VIEWER_TABLE")
        elif dbtype == 'SQLite':
            db_dict = {'system': section,
                       'type': dbtype}
            if parser.has_option(section, "FILE"):
                db_dict['file'] = self.convert_path(parser.get(section, "FILE"))
            if parser.has_option(section, "DATABASE_FILE"):
                db_dict['file'] = self.convert_path(parser.get(section, "DATABASE_FILE"))
            else:
                db_dict['file'] = self.convert_path(os.path.join(resource_paths[0], 'sqlite.db'))
            if parser.has_option(section, "JOB_TABLE"):
                db_dict['table_name'] = parser.get(section, "JOB_TABLE")
            else:
                db_dict['table_name'] = 'jobs_pyiron'
        elif dbtype == 'SQLalchemy':
            db_dict = {'system': section,
                       'type': dbtype,
                       'connection_string': parser.get(section, "CONNECTION"),
                       'table_name': parser.get(section, "JOB_TABLE")}
        elif dbtype == 'MySQL':
            db_dict = {'system': section,
                       'type': dbtype,
                       'database': parser.get(section, "NAME"),
                       'user': parser.get(section, "USER"),
                       'password': parser.get(section, "PASSWD"),
                       'table_name': parser.get(section, "JOB_TABLE"),
                       'host': parser.get(section, "HOST")}
        else:
            raise ValueError('Database configuration unreadable!')
        db_dict['top_level_dirs'] = top_level_dirs
        db_dict['resource_paths'] = resource_paths
        if parser.has_option(section, "USER"):
            db_dict['user'] = parser.get(section, "USER")
        elif parser.has_option('DEFAULT', "USER"):
            db_dict['user'] = parser.get('DEFAULT', 'LOGIN_USER')
        else:
            db_dict['user'] = 'pyiron'
        return db_dict

    @staticmethod
    def convert_path(path):
        return Path(path).expanduser().resolve().absolute().as_posix()