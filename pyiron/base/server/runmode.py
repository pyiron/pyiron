# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
Runmode class defines the different modes a pyiron job can be executed in
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


run_mode_lst = [
    "modal",
    "non_modal",
    "queue",
    "manual",
    "thread",
    "interactive",
    "interactive_non_modal",
]


class Runmode(object):
    """
    Run mode describes how the job is going to be executed:
    - modal: the interactive run mode
    - non_modal: sending the job to the background on the same machine
    - queue: submit the job to the queuing system
    - manual: let the user manually execute the job
    - thread: internal job mode, which is selected when the master job is send to the queue.
    - interactive: the interactive run mode
    Args:
        mode (str): ['modal', 'non_modal', 'queue', 'manual', 'thread', 'interactive']
    """

    def __init__(self, mode="modal"):
        super(Runmode, self).__init__()
        self._reset_mode()
        self.mode = mode

    @property
    def mode(self):
        """
        Get the run_mode of the job
        Returns:
            str: ['modal', 'non_modal', 'queue', 'manual', 'thread', 'interactive']
        """
        return [key for key, val in self._mode.items() if val][0]

    @mode.setter
    def mode(self, new_mode):
        """
        Set the run_mode of the job
        Args:
            new_mode (str): ['modal', 'non_modal', 'queue', 'manual', 'thread', 'interactive']
        """
        if isinstance(new_mode, str) and new_mode in self._mode.keys():
            self._reset_mode()
            self._mode[new_mode] = True

    def _reset_mode(self):
        super(Runmode, self).__setattr__(
            "_mode", {run_mode: False for run_mode in run_mode_lst}
        )

    def __repr__(self):
        return repr(self.mode)

    def __str__(self):
        return str(self.mode)

    def __getattr__(self, name):
        if name in self._mode.keys():
            return self._mode[name]
        else:
            raise AttributeError

    def __setattr__(self, name, value):
        if name in self._mode.keys():
            if not isinstance(value, bool):
                raise TypeError("A run mode can only be activated using [True].")
            if value:
                self.mode = name
            else:
                raise ValueError("A run mode can only be activated using [True].")
        else:
            super(Runmode, self).__setattr__(name, value)

    def __dir__(self):
        return list(self._mode.keys())

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return other._mode == self._mode
        elif isinstance(other, str):
            return other == self.mode
        else:
            return super(Runmode, self).__eq__(other)

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return other._mode != self._mode
        elif isinstance(other, str):
            return other != self.mode
        else:
            return super(Runmode, self).__ne__(other)
