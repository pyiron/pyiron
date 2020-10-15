# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.job.interactivewrapper import InteractiveWrapper, ReferenceJobOutput
from pyiron_base import InputList, Settings

__author__ = "Osamu Waseda"
__copyright__ = "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH " \
                "- Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Osamu Waseda"
__email__ = "waseda@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"

s = Settings()


class ART(object):
    def __init__(self, art_id, direction, gamma=0.1, fix_layer=False, non_art_id=None):
        if int(art_id) != art_id or art_id < 0:
            raise ValueError('art_id must be a posive integer')
        if len(direction) != 3 or np.isclose(np.linalg.norm(direction), 0):
            raise ValueError('direction must be a finite 3d vector')
        if gamma < 0:
            raise ValueError('gamma must be a positive float')
        if fix_layer and non_art_id is not None:
            raise ValueError('fix_layer and non_art_id cannot be set at the same time')
        self.art_id = art_id
        self.direction = direction
        self.gamma = gamma
        self.non_art_id = non_art_id
        if non_art_id is not None:
            self.non_art_id = np.array([non_art_id]).flatten()
        self.fix_layer = fix_layer

    @property
    def _R(self):
        value = np.array(self.direction)
        value = value / np.linalg.norm(value)
        return np.outer(value, value)

    def get_forces(self, f_in):
        f = np.array(f_in)
        if len(f.shape) == 2:
            f = np.array([f])
        if self.non_art_id is None:
            self.non_art_id = np.arange(len(f[0])) != self.art_id
        f_art = (1.0+self.gamma)*np.einsum('ij,nj->ni', self._R, f[:, self.art_id])
        if self.fix_layer:
            f[:, self.non_art_id] = np.einsum('nmj,ij->nmi', f[:, self.non_art_id], np.identity(3)-self._R)
        else:
            f[:, self.non_art_id] += f_art[:, np.newaxis, :] / np.sum(self.non_art_id != False)
        f[:, self.art_id] -= f_art
        return f.reshape(np.array(f_in).shape)


class ARTInteractive(InteractiveWrapper):
    """
    Apply an artificial force according to the Activation Relaxation Technique (ART)
    DOI:https://doi.org/10.1103/PhysRevE.57.2419

    The applied force f_art is calculated from the original force f by:

    f_art = f-(1+gamma)*np.dot(n,f)*n

    where gamma is a parameter to be determined in the input (default: 0.1) and n is
    the direction along which the force is reversed (3d-vector). In order to homogenize
    the total force in the system, f_art is distributed among atoms specified by
    non_art_id (default: all atoms), or if fix_layer is defined, the forces acting on
    all the other atoms along the direction n are cancelled. Note: Since the energy
    is not compatible with the forces anymore, structure optimization methods which
    rely on the energy variation (such as conjugate gradient) cannot/should not be used.

    Input:
        - art_id (int): atom id on which ART force is applied
        - direction (list/numpy.ndarray): direction along which force is reversed
        - gamma (float): prefactor for force inversion. v.s.
        - non_art_id (list/numpy.ndarray): list of atoms to be used for the force cancellation
        - fix_layer (bool): whether or not to fix all other atoms on the layer perpendicular
            to the direction along which force is reversed.

    Example:

        # Structure creation
        >>> vacancy_position = structure.positions[0]
        >>> del structure[0]
        >>> neighbors = structure.get_neighborhood(vacancy_position)
        >>> direction = neighbors.vecs[0]
        >>> art_id = neighbors.indices[0]
        >>> structure.positions[art_id] -= 0.5*direction

        # Job creation
        >>> some_atomistic_job.structure = structure
        >>> art = ARTInteractive(job_name='art')
        >>> art.ref_job = some_atomistic_job
        >>> art.input.art_id = art_id
        >>> art.input.direction = direction
        >>> some_minimizer.ref_job = art
        >>> some_minimizer.run()


    """
    def __init__(self, project, job_name):
        super(ARTInteractive, self).__init__(project, job_name)
        self.__name__ = "ARTInteractive"
        self.input = InputList(table_name='custom_dict')
        self.input.gamma = 0.1
        self.input.fix_layer = False
        self.input.non_art_id = None
        self.input.art_id = None
        self.input.direction = None
        self.output = ARTIntOutput(job=self)
        self.server.run_mode.interactive = True
        self._interactive_interface = None
        self._art = None

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        self.input.read_only = True

    @property
    def art(self):
        if self._art is None:
            self._art = ART(
                art_id=self.input.art_id,
                direction=self.input.direction,
                gamma=self.input.gamma,
                fix_layer=self.input.fix_layer,
                non_art_id=self.input.non_art_id
            )
        return self._art

    def run_if_interactive(self):
        self._logger.debug('art status: ' + str(self.status))
        if not self.status.running:
            self.ref_job_initialize()
        self.status.running = True
        if self.ref_job.server.run_mode.interactive:
            self.ref_job.run()
        else:
            self.ref_job.run(run_again=True)
        self._logger.debug('art status: ' + str(self.status))

    def interactive_forces_getter(self):
        return self.art.get_forces(self.ref_job.output.forces[-1])

    def validate_ready_to_run(self):
        """
            check whether art_id and direction are set
        """
        if self.input.art_id is None or self.input.direction is None:
            raise AssertionError('art_id and/or direction not set')

    def interactive_close(self):
        self.status.collect = True
        if self.ref_job.server.run_mode.interactive:
            self.ref_job.interactive_close()
        self.status.finished = True


class ARTIntOutput(ReferenceJobOutput):
    def __init__(self, job):
        super(ARTIntOutput, self).__init__(job=job)

    @property
    def forces(self):
        return self._job.art.get_forces(self._job.ref_job.output.forces)
