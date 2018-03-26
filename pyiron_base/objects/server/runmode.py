# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

"""
Runmode class defines the different modes a pyiron job can be executed in
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class Runmode(str):
    """
    Run mode describes how the job is going to be executed:
    - modal: the interactive run mode
    - non_modal: sending the job to the background on the same machine
    - queue: submit the job to the queuing system
    - manual: let the user manually execute the job
    - thread: internal job mode, which is selected when the master job is send to the queue.

    Args:
        mode (str): ['modal', 'non_modal', 'queue', 'manual', 'thread']
    """
    def __init__(self, mode='modal'):
        super(Runmode, self).__init__()
        self._mode = None
        self.mode = mode

    @property
    def modal(self):
        """
        Check if job is set for 'modal' mode - modal: the interactive run mode

        Returns:
            bool: [True/False]
        """
        return self._mode['modal']

    @modal.setter
    def modal(self, var_bool):
        """
        Set the job to 'modal' mode - modal: the interactive run mode

        Args:
            var_bool (bool): [True/False]
        """
        self._mode_setter(mode='modal', var_bool=var_bool)

    @property
    def thread(self):
        """
        internal job mode: Check if job is set for 'thread' mode - thread: is selected when the master job is send to
        the queue.

        Returns:
            bool: [True/False]
        """
        return self._mode['thread']

    @thread.setter
    def thread(self, var_bool):
        """
        internal job mode: Set the job to 'thread' mode - thread: is selected when the master job is send to the queue.

        Args:
            var_bool (bool): [True/False]
        """
        self._mode_setter(mode='thread', var_bool=var_bool)

    @property
    def non_modal(self):
        """
        Check if job is set for 'non_modal' mode - non_modal: sending the job to the background on the same machine

        Returns:
            bool: [True/False]
        """
        return self._mode['non_modal']

    @non_modal.setter
    def non_modal(self, var_bool):
        """
        Set the job to 'non_modal' mode - non_modal: sending the job to the background on the same machine

        Args:
            var_bool (bool): [True/False]
        """
        self._mode_setter(mode='non_modal', var_bool=var_bool)

    @property
    def queue(self):
        """
        Check if job is set for 'queue' mode - queue: submit the job to the queuing system

        Returns:
            bool: [True/False]
        """
        return self._mode['queue']

    @queue.setter
    def queue(self, var_bool):
        """
        Set the job to 'queue' mode - queue: submit the job to the queuing system

        Args:
            var_bool (bool): [True/False]
        """
        self._mode_setter(mode='queue', var_bool=var_bool)

    @property
    def manual(self):
        """
        Check if job is set for 'manual' mode - manual: let the user manually execute the job

        Returns:
            bool: [True/False]
        """
        return self._mode['manual']

    @manual.setter
    def manual(self, var_bool):
        """
        Set the job to 'manual' mode - manual: let the user manually execute the job

        Args:
            var_bool (bool): [True/False]
        """
        self._mode_setter(mode='manual', var_bool=var_bool)

    @property
    def mode(self):
        """
        Get the run_mode of the job

        Returns:
            str: ['modal', 'non_modal', 'queue', 'manual', 'thread']
        """
        for key, val in self._mode.items():
            if val:
                return key

    @mode.setter
    def mode(self, new_mode):
        """
        Set the run_mode of the job

        Args:
            new_mode (str): ['modal', 'non_modal', 'queue', 'manual', 'thread']
        """
        self._reset()
        if isinstance(new_mode, str) and new_mode in self._mode.keys():
            self._mode[new_mode] = True

    def _mode_setter(self, mode, var_bool=True):
        """
        internal function to set the run_mode

        Args:
            var_bool (bool): True
            mode (str): ['modal', 'non_modal', 'queue', 'manual', 'thread']
        """
        if not isinstance(var_bool, bool):
            raise TypeError('A run mode can only be activated using [True].')
        if var_bool:
            self._reset()
            self._mode[mode] = True
        else:
            raise ValueError('A run mode can only be activated using [True].')

    def _reset(self):
        """
        internal function to reset the run mode - sets all run modes to false.
        """
        self._mode = {'modal': False,
                      'non_modal': False,
                      'queue': False,
                      'manual': False,
                      'thread': False}

    def __repr__(self):
        return repr(self.mode)

    def __str__(self):
        return str(self.mode)
