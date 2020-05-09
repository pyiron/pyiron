# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os

"""
Executable class loading executables from static/bin/<code>/
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


class Executable(object):
    def __init__(
        self,
        path_binary_codes,
        codename=None,
        module=None,
        code=None,
        overwrite_nt_flag=False,
    ):
        """
        Handle the path to the executable, as well as the version selection.

        Args:
            codename (str): name of the code str
            path_binary_codes (list): path to the binary codes as an absolute path
            overwrite_nt_flag (bool):
        """
        self.__version__ = None
        if code is not None:  # Backwards compatibility
            if not isinstance(code.__name__, str):
                raise TypeError("The codename should be a string.")
            codename = code.__name__
            module = code.__module__.split(".")[1]
        if codename is not None and module is not None:
            self.__name__ = codename.lower()
            code_path_lst = [
                os.path.join(path, module, "bin") for path in path_binary_codes
            ]
            backwards_compatible_path_lst = [
                os.path.join(path, self.__name__) for path in path_binary_codes
            ]
            self._path_bin = [
                exe_path
                for exe_path in (code_path_lst + backwards_compatible_path_lst)
                if os.path.exists(exe_path)
            ]
        else:  # Backwards compatibility
            self.__name__ = codename.lower()
            self._path_bin = [
                os.path.join(path, self.__name__)
                for path in path_binary_codes
                if os.path.exists(os.path.join(path, self.__name__))
            ]
        if overwrite_nt_flag:
            self._operation_system_nt = False
        else:
            self._operation_system_nt = os.name == "nt"
        self._executable_lst = self._executable_versions_list()
        self._executable = None
        self._executable_path = None
        self._mpi = False
        if self._executable_lst:
            self.version = self.default_version

    @property
    def version(self):
        """
        Version of the Executable

        Returns:
            str: version
        """
        return self.__version__

    @property
    def default_version(self):
        """
        Default Version of the Available Executables
        i.e. specifically defined

        Returns:
            str: default_version
        """
        for executable in self._executable_lst.keys():
            if "default" in executable and "mpi" not in executable:
                return executable
        return sorted(self._executable_lst.keys())[0]

    @version.setter
    def version(self, new_version):
        """
        Version of the Executable

        Args:
            new_version (str): version
        """
        if new_version in self._executable_lst.keys():
            self.__version__ = new_version
            if "mpi" in new_version:
                self._mpi = True
            self._executable_path = None
        else:
            raise ValueError(
                "Version  [%s] is not supported, please choose one of the following versions: "
                % new_version,
                str(self.available_versions),
            )

    @property
    def mpi(self):
        """
        Check if the message processing interface is activated.

        Returns:
            bool: [True/False]
        """
        if not self._mpi and self.version and "_mpi" in self.version:
            self._mpi = True
        return self._mpi

    @mpi.setter
    def mpi(self, mpi_bool):
        """
        Activate the message processing interface.

        Args:
            mpi_bool (bool): [True/False]
        """
        if not isinstance(mpi_bool, bool):
            raise TypeError("MPI can either be enabled or disabled: [True/False]")
        if self.version and "_mpi" not in self.version:
            self.version += "_mpi"
        if self.version is None and self.executable_path is None:
            raise ValueError("No executable set!")

    @property
    def available_versions(self):
        """
        List all available exectuables in the path_binary_codes for the specified codename.

        Returns:
            list: list of the available version
        """
        return self.list_executables()

    def list_executables(self):
        """
        List all available exectuables in the path_binary_codes for the specified codename.

        Returns:
            list: list of the available version
        """
        return sorted(list(self._executable_lst.keys()))

    @property
    def executable_path(self):
        """
        Get the executable path

        Returns:
            str: absolute path
        """
        if self._executable_path is not None:
            if os.name == "nt":
                return self._executable_path.replace("\\", "/")
            else:
                return self._executable_path
        return self._executable_select()

    @executable_path.setter
    def executable_path(self, new_path):
        """
        Set the executable path

        Args:
            new_path: absolute path
        """
        self.__version__ = new_path
        self._executable_path = new_path
        if new_path and "mpi" in new_path:
            self._mpi = True
        else:
            self._mpi = False

    def __repr__(self):
        """
        Executable path
        """
        return repr(self.executable_path)

    def __str__(self):
        """
        Executable path
        """
        return self.executable_path

    def _executable_versions_list(self):
        """
        Internal function to list all available exectuables in the path_binary_codes for the specified codename.

        Returns:
            dict: list of the available version
        """
        if self._operation_system_nt:
            extension = ".bat"
        else:
            extension = ".sh"
        try:
            executable_dict = {}
            for path in self._path_bin:
                for executable in os.listdir(path):
                    if (
                        executable.startswith("run_" + self.__name__ + "_")
                        & executable.endswith(extension)
                        and executable[
                            len("run_" + self.__name__) + 1 : -len(extension)
                        ]
                        not in executable_dict.keys()
                    ):
                        executable_dict[
                            executable[
                                len("run_" + self.__name__) + 1 : -len(extension)
                            ]
                        ] = os.path.join(path, executable).replace("\\", "/")
            return executable_dict
        except OSError:  # No executable exists - This is the case for GenericJob and other abstract job classes.
            return dict()

    def _executable_select(self):
        """
        Internal function to select an executable based on the codename and the version.

        Returns:
            str: absolute executable path
        """
        try:
            return self._executable_lst[self.version]
        except KeyError:
            return ""
