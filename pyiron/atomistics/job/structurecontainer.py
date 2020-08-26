# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

__author__ = "Yury Lysogorskiy, Jan Janssen, Marvin Poul"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "0.1"
__maintainer__ = "Marvin Poul"
__email__ = "poul@mpie.de"
__status__ = "development"
__date__ = "Aug 12, 2020"

from pyiron.base.generic.inputlist import InputList
from pyiron.base.job.generic import GenericJob
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron.atomistics.structure.atoms import Atoms

class StructureContainer(AtomisticGenericJob):
    """
    Container to save a list of structures in HDF5 together with tags.

    Add new structures with :meth:`.append`, they are added to
    :attr:`.structure_lst`.  The HDF5 is written when :meth:`.run` is called.
    """

    __version__ = "0.1.0"

    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        self._structure_lst = InputList(table_name = "structures")
        self.server.run_mode.interactive = True

    @property
    def structure_lst(self):
        """
        :class:`.InputList`: list of :class:`~.Atoms`
        """
        return self._structure_lst

    @staticmethod
    def _to_structure(structure_or_job):
        """
        Return structure from structure or atomic job.

        Args:
            structure_or_job (:class:`~.AtomisticGenericJob`, :class:`~.Atoms`):
                if :class:`~.AtomisticGenericJob` try to get most recent structure,
        copy it and set the job_id in :attr:`~.Atoms.info`
        
        Returns:
            :class:`~.Atoms`: structure from the job or given structure
            
        Raises:
            ValueError: if given :class:`~.AtomisticGenericJob` has no structure set
            TypeError: if structure_or_job is of invalid type
        """
        if isinstance(structure_or_job, AtomisticGenericJob):
            if structure_or_job.structure:
                s = structure_or_job.get_structure(-1).copy()
                s.info["jobid"] = structure_or_job.job_id
                return s
            else:
                raise ValueError(
                        "The job does not contain any structure to import."
                )
        elif isinstance(structure_or_job, Atoms):
            return structure_or_job
        else:
            raise TypeError(
                "You can only use a structure object or an "
                "AtomisticGenericJob object."
            )

    @property
    def structure(self):
        """
        :class:`~.Atoms`: first (or only) structure set in the container

        :setter:  :class:`~.Atoms`, :class:`~.AtomisticGenericJob`
            if a job is given take the last structure and add the job id to its
            :attr:`pyiron.atomistics.structure.Atoms.info`
        """
        return self.structure_lst.get(0, None)

    @structure.setter
    def structure(self, structure_or_job):
        item = self._to_structure(structure_or_job)
        if len(self.structure_lst) >= 1:
            self.structure_lst[0] = item
        else:
            self.structure_lst.append(item)

    def append(self, structure_or_job):
        """
        Add new structure to structure list.

        The added structure will available in :attr:`~.structure_lst`.  If the
        structure is added via a job, retrieve the latest structure and add its
        id to :attr:`pyiron.atomistics.generic.Atoms.info`.

        Args:
            structure_or_job (:class:`~.AtomisticGenericJob`/:class:`~.Atoms`):
                if :class:`~.AtomisticGenericJob` add from
                :meth:`~.AtomisticGenericJob.get_structure`,
                otherwise add just the given :class:`~.Atoms`

        Returns:
            dict: item added to :attr:`~.structure_lst`
        """
        self.structure_lst.append(self._to_structure(structure_or_job))
        return self.structure_lst[0]

    def run_static(self):
        self.status.finished = True

    def run_if_interactive(self):
        self.to_hdf()
        self.status.finished = True

    def write_input(self):
        pass

    def collect_output(self):
        pass

    def to_hdf(self, hdf = None, group_name = None):
        # skip any of the AtomisticGenericJob specific serialization, since we
        # handle the structures on our own and that method might just confuse
        # self.structure and self.structure_lst
        GenericJob.to_hdf(self, hdf = hdf, group_name = group_name)

        hdf = self.project_hdf5.create_group("structures")

        for i, structure in enumerate(self.structure_lst.values()):
            structure.to_hdf(hdf, group_name = "structure_{}".format(i))

    def from_hdf(self, hdf = None, group_name = None):
        GenericJob.from_hdf(self, hdf = hdf, group_name = group_name)

        self.structure_lst.clear()

        hdf = self.project_hdf5["structures"]
        for group in sorted(hdf.list_groups()):
            structure = Atoms()
            structure.from_hdf(hdf, group_name = group)
            self.structure_lst.append(structure)
