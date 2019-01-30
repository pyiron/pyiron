# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron.atomistics.master.serial import SerialMaster

__author__ = "Yury Lysogorskiy"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


def convergence_goal(self, **qwargs):
    import numpy as np
    eps = 0.01
    if "eps" in qwargs:
        eps = qwargs["eps"]
    latt_in = self[-1].structure.get_volume()**(1/3)
    latt_out = self[-1].get_structure().get_volume()**(1/3)
    if np.max(np.abs(latt_out-latt_in)) < eps:
        return True
    else:
        return self.create_next()


class ConvergenceVolume(SerialMaster):
    """

    Args:
        project: 
        job_name: 
    """
    def __init__(self, project, job_name):
        super(ConvergenceVolume, self).__init__(project, job_name=job_name)
        self.__name__ = "ConvergenceVolume"
        self.__version__ = '0.0.2'
        if not self["input/convergence_goal"]:
            self.set_goal(convergence_goal, eps=0.005)

    def create_next(self, job_name=None):
        """
        
        Args:
            job_name: 

        Returns:

        """
        ham_first_job_name = self.project.db.get_item_by_id(self.child_ids[0])['job']
        if job_name is None:
            job_name = ham_first_job_name + "_" + str(len(self))
        return super(ConvergenceVolume, self).create_next(job_name=job_name)
