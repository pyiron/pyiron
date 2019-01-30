# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron.atomistics.master.serial import SerialMaster

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


def convergence_goal(self, eps=0.005):
    import numpy as np
    if len(self) > 1:
        prev_job_eng = self[-2]['output/generic/energy_tot']
        last_job_eng = self[-1]['output/generic/energy_tot']
        if np.abs(prev_job_eng-last_job_eng) < eps:
            return True
    ham_prev = self[-1]
    encut_new = ham_prev.get_encut() + 100
    job_name = self.get_initial_child_name() + "_" + str(encut_new)
    ham_next = ham_prev.restart(job_name=job_name)
    ham_next.set_encut(encut_new)
    return ham_next


class ConvEncutSerial(SerialMaster):
    """

    Args:
        project:
        job_name:
    """

    def __init__(self, project, job_name):
        super(ConvEncutSerial, self).__init__(project, job_name=job_name)
        self.__name__ = "ConvEncutSerial"
        self.__version__ = '0.0.2'
        if not self["input/convergence_goal"] or self["input/convergence_goal"] == 'None':
            self.set_goal(convergence_goal, eps=0.005)

    def create_next(self, job_name=None):
        """

        Args:
            job_name:

        Returns:

        """
        ham_prev = self[-1]
        encut_new = ham_prev.get_encut() + 100
        job_name = self.get_initial_child_name() + "_" + str(encut_new)
        ham_next = ham_prev.restart(job_name=job_name)
        ham_next.set_encut(encut_new)
        return ham_next
