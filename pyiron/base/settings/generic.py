# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from builtins import input
import os
import importlib
from six import with_metaclass
import sys
from configparser import ConfigParser
from pathlib2 import Path
from pyiron.base.settings.logger import setup_logger
from pyiron.base.database.generic import DatabaseAccess
from pyiron.base.settings.install import install_pyiron

"""
The settings file provides the attributes of the configuration as properties.
"""

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
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
        if (
            kwargs is not None
            and "config" in kwargs.keys()
            and kwargs["config"] is not None
        ):
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
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
        # Default config dictionary
        self._configuration = {
            "user": "pyiron",
            "resource_paths": ["~/pyiron/resources"],
            "project_paths": ["~/pyiron/projects"],
            "sql_connection_string": None,
            "sql_table_name": "jobs_pyiron",
            "sql_view_connection_string": None,
            "sql_view_table_name": None,
            "sql_view_user": None,
            "sql_view_user_key": None,
            "sql_file": None,
            "sql_host": None,
            "sql_type": "SQLite",
            "sql_user_key": None,
            "sql_database": None,
            "project_check_enabled": True,
            "disable_database": False,
        }
        environment = os.environ
        if "PYIRONCONFIG" in environment.keys():
            config_file = environment["PYIRONCONFIG"]
        else:
            config_file = os.path.expanduser(os.path.join("~", ".pyiron"))
        if os.path.isfile(config_file):
            self._config_parse_file(config_file)
        elif any(["PYIRON" in e for e in environment.keys()]):
            self._configuration = self.get_config_from_environment(
                environment=environment,
                config=self._configuration
            )
        else:
            print("Fall back to default configuration: "
                  "{'resource_paths': ['~/pyiron/resources'], "
                  "'project_paths': ['~/pyiron/projects']}")

        # Take dictionary as primary source - overwrite everything
        self._read_external_config(config=config)

        self._configuration["project_paths"] = [
            convert_path(path) + "/" if path[-1] != "/" else convert_path(path)
            for path in self._configuration["project_paths"]
        ]
        self._configuration["resource_paths"] = [
            convert_path(path) for path in self._configuration["resource_paths"]
        ]

        # Build the SQLalchemy connection strings
        if not self.database_is_disabled:
            self._configuration = self.convert_database_config(
                config=self._configuration
            )
        self._database = None
        self._use_local_database = False
        self._queue_adapter = None
        self._queue_adapter = self._init_queue_adapter(
            resource_path_lst=self._configuration["resource_paths"]
        )
        self.logger = setup_logger()
        self._publication_lst = {}
        self.publication_add(self.publication)

    @property
    def database(self):
        return self._database

    @property
    def queue_adapter(self):
        return self._queue_adapter

    @property
    def project_check_enabled(self):
        return self._configuration["project_check_enabled"]

    @property
    def database_is_disabled(self):
        return self._configuration["disable_database"]

    @property
    def publication_lst(self):
        """
        List of publications currently in use.

        Returns:
            list: list of publications
        """
        all_publication = []
        for v in self._publication_lst.values():
            if isinstance(v, list):
                all_publication += v
            else:
                all_publication.append(v)
        return all_publication

    def publication_add(self, pub_dict):
        """
        Add a publication to the list of publications

        Args:
            pub_dict (dict): The key should be the name of the code used and the value a list of publications to cite.
        """
        for key, value in pub_dict.items():
            if key not in self._publication_lst.keys():
                self._publication_lst[key] = value

    @property
    def login_user(self):
        """
        Get the username of the current user

        Returns:
            str: username
        """
        return self._configuration["user"]

    @property
    def resource_paths(self):
        """
        Get the path where the potentials for the individual Hamiltons are located

        Returns:
            list: path of paths
        """
        return self._configuration["resource_paths"]

    def open_connection(self):
        """
        Internal function to open the connection to the database. Only after this function is called the database is
        accessable.
        """
        if self._database is None and not self.database_is_disabled:
            self._database = DatabaseAccess(
                self._configuration["sql_connection_string"],
                self._configuration["sql_table_name"],
            )

    def switch_to_local_database(self, file_name="pyiron.db", cwd=None):
        """
        Swtich to an local SQLite based database.

        Args:
            file_name (str): SQLite database file name
            cwd (str/None): directory where the SQLite database file is located in
        """
        if not self._use_local_database and not self.database_is_disabled:
            if cwd is None and not os.path.isabs(file_name):
                file_name = os.path.join(os.path.abspath(os.path.curdir), file_name)
            elif cwd is not None:
                file_name = os.path.join(cwd, file_name)
            self.close_connection()
            self._database = DatabaseAccess(
                "sqlite:///" + file_name, self._configuration["sql_table_name"]
            )
            self._use_local_database = True
        else:
            print("Database is already in local mode or disabled!")

    def switch_to_central_database(self):
        """
        Switch to central database
        """
        if self._use_local_database and not self.database_is_disabled:
            self.close_connection()
            self._database = DatabaseAccess(
                self._configuration["sql_connection_string"],
                self._configuration["sql_table_name"],
            )
            self._use_local_database = False
        else:
            print("Database is already in central mode or disabled!")

    def switch_to_viewer_mode(self):
        """
        Switch from user mode to viewer mode - if viewer_mode is enable pyiron has read only access to the database.
        """
        if self._configuration["sql_view_connection_string"] is not None and not self.database_is_disabled:
            if not self._database.viewer_mode:
                self.close_connection()
                self._database = DatabaseAccess(
                    self._configuration["sql_view_connection_string"],
                    self._configuration["sql_view_table_name"],
                )
                self._database.viewer_mode = True
            else:
                print("Database is already in viewer mode!")
        else:
            print("Viewer Mode is not available on this pyiron installation.")

    def switch_to_user_mode(self):
        """
        Switch from viewer mode to user mode - if viewer_mode is enable pyiron has read only access to the database.
        """
        if self._configuration["sql_view_connection_string"] is not None and not self.database_is_disabled:
            if self._database.viewer_mode:
                self.close_connection()
                self._database = DatabaseAccess(
                    self._configuration["sql_connection_string"],
                    self._configuration["sql_table_name"],
                )
                self._database.viewer_mode = True
            else:
                print("Database is already in user mode!")
        else:
            print("Viewer Mode is not available on this pyiron installation.")

    def close_connection(self):
        """
        Internal function to close the connection to the database.
        """
        if hasattr(self, "_database") and self._database is not None:
            self._database.conn.close()
            self._database = None

    def top_path(self, full_path):
        """
        Validated that the full_path is a sub directory of one of the pyrion environments loaded.

        Args:
            full_path (str): path

        Returns:
            str: path
        """
        if full_path[-1] != "/":
            full_path += "/"
        if not self.project_check_enabled:
            return None
        for path in self._configuration["project_paths"]:
            if path in full_path:
                return path
        raise ValueError(
            "the current path {0} is not included in the .pyiron configuration. {1}".format(
                full_path, self._configuration["project_paths"]
            )
        )

    # private functions
    @staticmethod
    def _init_queue_adapter(resource_path_lst):
        """
        Initialize the queue adapter if a folder queues is found in one of the resource paths which contains a
        queue configuration file (queue.yaml).

        Args:
            resource_path_lst (list): List of resource paths

        Returns:
            pysqa.QueueAdapter:
        """
        for resource_path in resource_path_lst:
            if (
                os.path.exists(resource_path)
                and "queues" in os.listdir(resource_path)
                and "queue.yaml" in os.listdir(os.path.join(resource_path, "queues"))
            ):
                queueadapter = getattr(importlib.import_module("pysqa"), "QueueAdapter")
                return queueadapter(directory=os.path.join(resource_path, "queues"))
        return None

    def _config_parse_file(self, config_file):
        """
        Read section in configuration file and return a dictionary with the corresponding parameters.

        Args:
            config_file(str): confi file to parse

        Returns:
            dict: dictionary with the environment configuration
        """
        # load config parser - depending on Python version
        parser = ConfigParser(inline_comment_prefixes=(";",))

        # read config
        parser.read(config_file)

        # load first section or default section [DEFAULT]
        if len(parser.sections()) > 0:
            section = parser.sections()[0]
        else:
            section = "DEFAULT"

        # identify SQL type
        if parser.has_option(section, "TYPE"):
            self._configuration["sql_type"] = parser.get(section, "TYPE")

        # read variables
        if parser.has_option(section, "PROJECT_PATHS"):
            self._configuration["project_paths"] = [
                convert_path(c.strip())
                for c in parser.get(section, "PROJECT_PATHS").split(",")
            ]
        elif parser.has_option(
            section, "TOP_LEVEL_DIRS"
        ):  # for backwards compatibility
            self._configuration["project_paths"] = [
                convert_path(c.strip())
                for c in parser.get(section, "TOP_LEVEL_DIRS").split(",")
            ]
        else:
            ValueError("No project path identified!")
        if parser.has_option(section, "PROJECT_CHECK_ENABLED"):
            self._configuration["project_check_enabled"] = \
                parser.getboolean(section, "PROJECT_CHECK_ENABLED")
        if parser.has_option(section, "DISABLE_DATABASE"):
            self._configuration["disable_database"] = \
                parser.getboolean(section, "DISABLE_DATABASE")
        if parser.has_option(section, "RESOURCE_PATHS"):
            self._configuration["resource_paths"] = [
                convert_path(c.strip())
                for c in parser.get(section, "RESOURCE_PATHS").split(",")
            ]
        if self._configuration["sql_type"] in ["Postgres", "MySQL"]:
            if (
                parser.has_option(section, "USER")
                & parser.has_option(section, "PASSWD")
                & parser.has_option(section, "HOST")
                & parser.has_option(section, "NAME")
            ):
                self._configuration["user"] = parser.get(section, "USER")
                self._configuration["sql_user_key"] = parser.get(section, "PASSWD")
                self._configuration["sql_host"] = parser.get(section, "HOST")
                self._configuration["sql_database"] = parser.get(section, "NAME")
                self._configuration["sql_file"] = None
            else:
                raise ValueError(
                    "If type Postgres or MySQL are selected the options USER, PASSWD, HOST and NAME are"
                    "required in the configuration file."
                )
            if (
                parser.has_option(section, "VIEWERUSER")
                & parser.has_option(section, "VIEWERPASSWD")
                & parser.has_option(section, "VIEWER_TABLE")
            ):
                self._configuration["sql_view_table_name"] = parser.get(
                    section, "VIEWER_TABLE"
                )
                self._configuration["sql_view_user"] = parser.get(section, "VIEWERUSER")
                self._configuration["sql_view_user_key"] = parser.get(
                    section, "VIEWERPASSWD"
                )
        elif self._configuration["sql_type"] == "SQLalchemy":
            self._configuration["sql_connection_string"] = parser.get(
                section, "CONNECTION"
            )
        else:  # finally we assume an SQLite connection
            if parser.has_option(section, "FILE"):
                self._configuration["sql_file"] = parser.get(section, "FILE").replace(
                    "\\", "/"
                )
            if parser.has_option(section, "DATABASE_FILE"):
                self._configuration["sql_file"] = parser.get(
                    section, "DATABASE_FILE"
                ).replace("\\", "/")
        if parser.has_option(section, "JOB_TABLE"):
            self._configuration["sql_table_name"] = parser.get(section, "JOB_TABLE")

    @staticmethod
    def convert_database_config(config):
        # Build the SQLalchemy connection strings
        if config["sql_type"] == "Postgres":
            config["sql_connection_string"] = (
                "postgresql://"
                + config["user"]
                + ":"
                + config["sql_user_key"]
                + "@"
                + config["sql_host"]
                + "/"
                + config["sql_database"]
            )
            if config["sql_view_user"] is not None:
                config["sql_view_connection_string"] = (
                    "postgresql://"
                    + config["sql_view_user"]
                    + ":"
                    + config["sql_view_user_key"]
                    + "@"
                    + config["sql_host"]
                    + "/"
                    + config["sql_database"]
                )
        elif config["sql_type"] == "MySQL":
            config["sql_connection_string"] = (
                "mysql+pymysql://"
                + config["user"]
                + ":"
                + config["sql_user_key"]
                + "@"
                + config["sql_host"]
                + "/"
                + config["sql_database"]
            )
        else:
            # SQLite is raising ugly error messages when the database directory does not exist.
            if config["sql_file"] is None:
                config["sql_file"] = "/".join(
                    [config["resource_paths"][0], "sqlite.db"]
                )
            if os.path.dirname(
                config["sql_file"]
            ) != "" and not os.path.exists(
                os.path.dirname(config["sql_file"])
            ):
                os.makedirs(os.path.dirname(config["sql_file"]))
            config[
                "sql_connection_string"
            ] = "sqlite:///" + config["sql_file"].replace("\\", "/")
        return config

    def _read_external_config(self, config):
        if isinstance(config, dict):
            for key, value in config.items():
                if key not in ["resource_paths", "project_paths"] or isinstance(
                    value, list
                ):
                    self._configuration[key] = value
                elif isinstance(value, str):
                    self._configuration[key] = [value]
                else:
                    TypeError(
                        "Config dictionary parameter type not recognized ", key, value
                    )

    @staticmethod
    def get_config_from_environment(environment, config):
        env_key_mapping = {
            "PYIRONUSER": "user",
            "PYIRONRESOURCEPATHS": "resource_paths",
            "PYIRONPROJECTPATHS": "project_paths",
            "PYIRONSQLCONNECTIONSTRING": "sql_connection_string",
            "PYIRONSQLTABLENAME": "sql_table_name",
            "PYIRONSQLVIEWCONNECTIONSTRING": "sql_view_connection_string",
            "PYIRONSQLVIEWTABLENAME": "sql_view_table_name",
            "PYIRONSQLVIEWUSER": "sql_view_user",
            "PYIRONSQLVIEWUSERKEY": "sql_view_user_key",
            "PYIRONSQLFILE": "sql_file",
            "PYIRONSQHOST": "sql_host",
            "PYIRONSQLTYPE": "sql_type",
            "PYIRONSQLUSERKEY": "sql_user_key",
            "PYIRONSQLDATABASE": "sql_database",
            "PYIRONPROJECTCHECKENABLED": "project_check_enabled",
            "PYIRONDISABLE": "disable_database",
        }
        for k, v in env_key_mapping.items():
            if k in environment.keys():
                if k in ["PYIRONPROJECTCHECKENABLED", "PYIRONDISABLE"]:
                    config[v] = environment[k].lower() in ['t', 'true', 'y', 'yes']
                elif k in ["PYIRONRESOURCEPATHS", "PYIRONPROJECTPATHS"]:
                    config[v] = environment[k].split(':')
                else:
                    config[v] = environment[k]
        return config

    @property
    def publication(self):
        return {
            "pyiron": {
                "pyiron-paper": {
                    "author": [
                        "Jan Janssen",
                        "Sudarsan Surendralal",
                        "Yury Lysogorskiy",
                        "Mira Todorova",
                        "Tilmann Hickel",
                        "Ralf Drautz",
                        "Jörg Neugebauer",
                    ],
                    "title": "pyiron: An integrated development environment for computational "
                    "materials science",
                    "journal": "Computational Materials Science",
                    "volume": "161",
                    "pages": "24 - 36",
                    "issn": "0927-0256",
                    "doi": "https://doi.org/10.1016/j.commatsci.2018.07.043",
                    "url": "http://www.sciencedirect.com/science/article/pii/S0927025618304786",
                    "year": "2019",
                }
            }
        }


def convert_path(path):
    """
    Convert path to POSIX path

    Args:
        path(str): input path

    Returns:
        str: absolute path in POSIX format
    """
    return os.path.abspath(os.path.expanduser(path)).replace("\\", "/")
