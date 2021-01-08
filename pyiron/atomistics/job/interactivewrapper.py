# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from datetime import datetime
import warnings
from pyiron_base import InteractiveWrapper as InteractiveWrapperBase
from pyiron_base.generic.util import deprecate
from pyiron.atomistics.structure.atoms import ase_to_pyiron
from pyiron.atomistics.structure.atoms import Atoms as PAtoms

__author__ = "Osamu Waseda, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class InteractiveWrapper(InteractiveWrapperBase):
    def __init__(self, project, job_name):
        super(InteractiveWrapper, self).__init__(project, job_name)

    @property
    def structure(self):
        if self.ref_job:
            return self._ref_job.structure
        else:
            return None

    @structure.setter
    def structure(self, basis):
        if self.ref_job:
            self._ref_job.structure = basis
        else:
            raise ValueError(
                "A structure can only be set after a start job has been assinged."
            )

    @deprecate("use get_structure() instead")
    def get_final_structure(self):
        """

        Returns:

        """
        if self.ref_job:
            return self._ref_job.get_structure(iteration_step=-1)
        else:
            return None

    def db_entry(self):
        """
        Generate the initial database entry

        Returns:
            (dict): db_dict
        """
        db_dict = super(InteractiveWrapper, self).db_entry()
        if self.structure:
            if isinstance(self.structure, PAtoms):
                parent_structure = self.structure.get_parent_basis()
            else:
                parent_structure = ase_to_pyiron(self.structure).get_parent_basis()
            db_dict["ChemicalFormula"] = parent_structure.get_chemical_formula()
        return db_dict


class ReferenceJobOutput(object):
    def __init__(self, job):
        self._job = job

    @property
    def indices(self):
        return self._job.ref_job.output.indices

    @property
    def cells(self):
        return self._job.ref_job.output.cells

    @property
    def energy_pot(self):
        return self._job.ref_job.output.energy_pot

    @property
    def energy_tot(self):
        return self._job.ref_job.output.energy_tot

    @property
    def forces(self):
        return self._job.ref_job.output.forces

    @property
    def positions(self):
        return self._job.ref_job.output.positions

    @property
    def pressures(self):
        return self._job.ref_job.output.pressures

    @property
    def steps(self):
        return self._job.ref_job.output.steps

    @property
    def temperatures(self):
        return self._job.ref_job.output.temperatures

    @property
    def time(self):
        return self._job.ref_job.output.time

    @property
    def unwrapped_positions(self):
        return self._job.ref_job.output.unwrapped_positions

    @property
    def volume(self):
        return self._job.ref_job.output.volume

    def __dir__(self):
        return list(set(list(self._job.ref_job.interactive_cache.keys())))
