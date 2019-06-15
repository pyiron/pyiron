# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import pickle
import subprocess


__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2018"


class LammpsLibrary(object):
    def __init__(self, cores=1):
        executable = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lmpmpi.py')
        # print(executable)
        self._process = subprocess.Popen(['mpiexec', '-n', str(cores), 'python', executable],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         stdin=subprocess.PIPE)

    def _send(self, command, data=None):
        """
        Send a command to the Lammps Library executable

        Args:
            command (str): command to be send to the
            data:
        """
        # print('send: ', {'c': command, 'd': data})
        pickle.dump({'c': command, 'd': data}, self._process.stdin)
        self._process.stdin.flush()

    def _receive(self):
        """
        Receive data from the Lammps library

        Returns:
            data
        """
        return pickle.load(self._process.stdout)

    def command(self, command):
        """
        Send a command to the lammps library

        Args:
            command (str):
        """
        self._send(command='command', data=command)

    def gather_atoms(self, *args):
        """
        Gather atoms from the lammps library

        Args:
            *args:

        Returns:
            np.array
        """
        self._send(command='gather_atoms', data=list(args))
        return self._receive()

    def scatter_atoms(self, *args):
        """
        Scatter atoms for the lammps library

        Args:
            *args:
        """
        self._send(command='scatter_atoms', data=list(args))

    def get_thermo(self, *args):
        """
        Get thermo from the lammps library

        Args:
            *args:

        Returns:

        """
        self._send(command='get_thermo', data=list(args))
        return self._receive()

    def extract_compute(self, *args):
        """
        Extract compute from the lammps library

        Args:
            *args:

        Returns:

        """
        self._send(command='extract_compute', data=list(args))
        return self._receive()

    def close(self):
        self._send(command='close')
