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
    def __init__(self):
        executable = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sub.py')
        # print(executable)
        self._process = subprocess.Popen(['python', executable],
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         stdin=subprocess.PIPE)

    def send(self, command, data=None):
        # print('send: ', {'c': command, 'd': data})
        pickle.dump({'c': command, 'd': data}, self._process.stdin)
        self._process.stdin.flush()

    def receive(self):
        pout = pickle.load(self._process.stdout)
        # print('receive: ', pout)
        return pout

    def command(self, command):
        self.send(command='command', data=command)

    def gather_atoms(self, *args):
        self.send(command='gather_atoms', data=list(args))
        return self.receive()

    def scatter_atoms(self, *args):
        self.send(command='scatter_atoms', data=list(args))

    def get_thermo(self, *args):
        self.send(command='get_thermo', data=list(args))
        return self.receive()

    def extract_compute(self, *args):
        self.send(command='extract_compute', data=list(args))
        return self.receive()

    def close(self):
        self.send(command='close')
