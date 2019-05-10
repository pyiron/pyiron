# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function

from copy import copy
import os
import posixpath
from pyiron.base.settings.generic import Settings
from six import string_types

"""
Classes for representing the file system path in pyiron
"""

__author__ = "Jan Janssen, Joerg Neugebauer"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class GenericPath(object):
    """
    Basic class to represent a project path in PyIron. A path consists of two parts, the root
    part which defines directory path where the project repository is located (top_level_path)
    and the project part which defines the relative path from the root_path to the project.

    This class is meant for storing and accessing a path, not for moving around which is done by
    the ProjectPath class.

    Args:
        root_path (str): absolute path name of the repository
        project_path (str): relative path to the specific project

    Attributes:

        .. attribute:: root_path

            the pyiron user directory, defined in the .pyiron configuration

        .. attribute:: project_path

            the relative path of the current project / folder starting from the root path
            of the pyiron user directory

        .. attribute:: path

            the absolute path of the current project / folder

        .. attribute:: base_name

            the name of the current project / folder

    Author: Jan Janssen
    """
    def __init__(self, root_path, project_path):
        self._root_path = None
        self._project_path = None
        self.root_path = root_path
        self.project_path = project_path

    @property
    def root_path(self):
        """
        the pyiron user directory, defined in the .pyiron configuration

        Returns:
            str: pyiron user directory of the current project
        """
        return self._root_path

    @root_path.setter
    def root_path(self, new_path):
        """
        the pyiron user directory, defined in the .pyiron configuration

        Args:
            new_path (str): new pyiron root path
        """
        self._root_path = self._windows_path_to_unix_path(new_path)

    @property
    def project_path(self):
        """
        the relative path of the current project / folder starting from the root path
        of the pyiron user directory

        Returns:
            str: relative path of the current project / folder
        """
        if self._project_path[-1] != '/':
            self._project_path += '/'
        return self._project_path


    @project_path.setter
    def project_path(self, new_path):
        """
        the relative path of the current project / folder starting from the root path
        of the pyiron user directory

        Args:
            new_path (str): new pyiron project path

        """
        self._project_path = self._windows_path_to_unix_path(posixpath.normpath(new_path))

    @property
    def path(self):
        """
        The absolute path to of the current object.

        Returns:
            str: current project path
        """
        return posixpath.join(self.root_path, self.project_path)

    @property
    def base_name(self):
        """
        The name of the current project folder

        Returns:
            str: name of the current project folder
        """
        if self.project_path[-1] in ['/', '\\']:
            return self.project_path.split("/")[-2]
        else:
            return self.project_path.split("/")[-1]

    def copy(self):
        """
        Copy the GenericPath object

        Returns:
            GenericPath: independent GenericPath object pointing to the same project folder
        """
        return copy(self)

    def __copy__(self):
        """
        Copy the GenericPath object

        Returns:
            GenericPath: independent GenericPath object pointing to the same project folder
        """
        return GenericPath(self.root_path, self.project_path)

    def __repr__(self):
        """
        String representation of the GenericPath object

        Returns:
            str: string representation of the root and the project path
        """
        project_str = "Project path: \n"
        project_str += "   root: " + str(self.root_path) + "\n"
        project_str += "   project: " + str(self.project_path) + "\n"
        return project_str

    def __str__(self):
        """
        String representation of the GenericPath object

        Returns:
            str: string representation of the absolute path
        """
        return self.path

    @staticmethod
    def _windows_path_to_unix_path(path):
        """
        Helperfunction to covert windows path into unix path

        Args:
            path (str): input path in windows or unix format

        Returns:
            str: output path in unix format
        """
        linux_path = path.replace("\\", "/")
        if linux_path[-1] != '/':
            linux_path += '/'
        return linux_path


class ProjectPath(GenericPath):
    def __init__(self, path):
        """
        Open a new or an existing project. The project is defined by providing a relative or an
        absolute path. If no path is provided the current working directory is used. This is the
        main class to walk inside the project structure, create new jobs, subprojects etc.
        Note: Changing the path has no effect on the current working directory

        Args:
            path (GenericPath, str): path of the project defined by GenericPath, absolute or relative (with respect to
                                     current working directory) path

        .. attribute:: root_path
            the pyiron user directory, defined in the .pyiron configuration
        .. attribute:: project_path
            the relative path of the current project / folder starting from the root path
            of the pyiron user directory
        .. attribute:: path
            the absolute path of the current project / folder
        .. attribute:: base_name
            the name of the current project / folder
        .. attribute:: history
            previously opened projects / folders
        """
        if path == "":
            raise ValueError('ProjectPath: path is not allowed to be empty!')
        generic_path = self._convert_str_to_generic_path(path)
        super(ProjectPath, self).__init__(generic_path.root_path, generic_path.project_path)
        self._history = []

    @property
    def history(self):
        """
        The history of the previously opened paths

        Returns:
            list: list of previously opened relative paths
        """
        return self._history

    def open(self, rel_path, history=True):
        """
        if rel_path exist set the project path to this directory
        if not create it and go there

        Args:
            rel_path (str): path relative to the current project path
            history (bool): By default pyiron stores a history of previously opened paths

        Returns:
            ProjectPath: New ProjectPath object pointing to the relative path
        """
        new_project = self.copy()
        new_project._create_path(new_project.path, rel_path)
        new_project.project_path = os.path.normpath(os.path.join(new_project.project_path, rel_path)).replace('\\', '/')
        if history:
            new_project.history.append(rel_path)
        return new_project

    def close(self):
        """
        return to the path before the last open if no history exists nothing happens
        """
        if self.history:
            path_lst = self.project_path.split("/")
            hist_lst = self.history[-1].split("/")
            self.project_path = "/".join(path_lst[:-len(hist_lst)])
            del self.history[-1]

    def copy(self):
        """
        copy the path without the history, i.e., to going back with close is not possible

        Returns:
            ProjectPath:
        """
        return ProjectPath(path=self.path)

    def removedirs(self, project_name=None):
        """
        equivalent to os.removedirs  -> remove empty dirs

        Args:
            project_name (str): relative path to the project folder to be deleted

        """
        try:
            if project_name:
                os.removedirs(os.path.join(self.path, project_name))
            else:
                os.removedirs(os.path.join(self.path))
        except OSError:
            pass

    def listdir(self):
        """
        equivalent to os.listdir
        list all files and directories in this path

        Returns:
            list: list of folders and files in the current project path
        """
        try:
            return os.listdir(self.path)
        except OSError:
            return []

    def walk(self):
        """
        equivalent to os.listdir
        list all files and directories in this path

        Returns:
            Generator: Directory tree generator.
        """
        return os.walk(self.path)

    def __repr__(self):
        """
        String representation of the ProjectPath object

        Returns:
            str: string representation of the root and the project path including the path history
        """
        project_str = super(ProjectPath, self).__repr__()
        if len(self.history) > 0:
            project_str += "   history: " + str(self.history) + "\n"
        return project_str

    def __enter__(self):
        """
        Helper function to support with calls

        Returns:
            ProjectPath: The object itself
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Helper function to support with calls

        Args:
            exc_type: Python internal
            exc_val: Python internal
            exc_tb: Python internal
        """
        self.close()

    def _convert_str_to_generic_path(self, path):
        """
        Convert path in string representation to an GenericPath object

        Args:
            path (str): absolute path

        Returns:
            GenericPath: GenericPath object pointing to the absolute path
        """
        if isinstance(path, GenericPath):
            return path
        elif isinstance(path, string_types):
            path = os.path.normpath(path)
            if not os.path.isabs(path):
                path_local = self._windows_path_to_unix_path(posixpath.abspath(os.curdir))
                self._create_path(path_local, path)
                path = posixpath.join(path_local, path)
            elif not os.path.exists(path) and os.path.exists(os.path.normpath(os.path.join(path, '..'))):
                self._create_path(path)
            # else: 
            #     raise ValueError(path, ' does not exist!')
            path = self._windows_path_to_unix_path(path)
            root_path, project_path = self._get_project_from_path(path)
            return GenericPath(root_path, project_path)
        else:
            raise TypeError('Only string and GenericPath objects are supported.')

    def _create_path(self, path, rel_path=None):
        """
        Create the directory if it does not exist already using os.makedirs()

        Args:
            path (str): absolute path
            rel_path (str): relative path starting from the absolute path (optional)

        """
        if rel_path:
            rel_path = self._windows_path_to_unix_path(rel_path)
            path = posixpath.join(path, rel_path)
        try:
            os.makedirs(path)
        except os.error:
            pass

    @staticmethod
    def _get_project_from_path(full_path):
        """
        Split the absolute path in root_path and project_path using the top_path function in Settings()

        Args:
            full_path (str): absolute path

        Returns:
            str, str: root_path, project_path
        """
        root = Settings().top_path(full_path)
        pr_path = posixpath.relpath(full_path, root)
        return root, pr_path
