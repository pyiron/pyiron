# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

__author__ = "Marvin Poul"
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

class StructureList(GenericJob):
    """
    Container to save a list of structures in HDF5 together with tags.

    Add new structures with :meth:`.StructureList.append`, they are
    added to :attr:`.StructureList.structures`.  The HDF5 is written when
    :meth:`.run` is called.
    """

    __version__ = "0.1.0"

    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        self.__structures = InputList(table_name = "structures")

    @property
    def structures(self):
        """
        :class:`.InputList`: list of structures
        with meta data; each item in this list has an 'atoms' key for the
        atomic structures and then as many additional keys as passed to
        :meth:`~pyiron.
        """
        return self.__structures

    def append(self, structure_or_job, **kwargs):
        """
        Add new structure to structure list.

        The added structure will available in
        :attr:`~.StructureList.structures`.  If **kwargs contains a key "atoms"
        is ignored, since the structure is saved under this key.  If the
        structure is added via a job, its id is also saved under the "jobid"
        key and **kwargs of the same name are ignored.

        Args:
            structure_or_job (AtomisticGenericJob/Atoms):
                if :class:`~.AtomisticGenericJob` add
                :attr:`~.AtomisticGenericJob.structure`,
                otherwise add just the given :class:`~.Atoms`
            **kwargs: meta data tags to add to the structure in the list
        """

        kwargs.pop("atoms", None)
        if isinstance(structure_or_job, AtomisticGenericJob):
            kwargs.pop("jobid", None)
            self.__structures.append({
                "atoms": structure_or_job.structure,
                "jobid": structure_or_job.id,
                **kwargs
            })
        elif isinstance(structure_or_job, Atoms):
            self.__structures.append({
                "atoms": structure_or_job,
                **kwargs
            })
        else:
            raise ValueError("structure_or_job must be atomistic job or Atoms")

    def run_static(self):
        self.to_hdf()
        self.status.finished = True

    def write_input(self):
        pass

    def collect_output(self):
        pass

    def to_hdf(self, hdf = None, group_name = None):
        super().to_hdf(hdf = hdf, group_name = group_name)

        hdf = self.project_hdf5.create_group("structures")

        for i, s in enumerate(self.__structures.values()):
            sub = hdf.create_group("index_{}".format(i))
            for k, v in s.items():
                if k == "atoms":
                    # use to_hdf from Atoms class
                    v.to_hdf(hdf = sub, group_name = "atoms")
                else:
                    # otherwise use the normal hdfio implementation
                    sub[k] = v

    def from_hdf(self, hdf = None, group_name = None):
        super().from_hdf(hdf = hdf, group_name = group_name)

        self.__structures.clear()

        hdf = self.project_hdf5["structures"]
        for i in sorted(hdf.list_groups()):
            sub = hdf[i]

            a = Atoms()
            a.from_hdf(hdf = sub, group_name = "atoms")
            meta = sub.list_groups()
            meta.remove("atoms")

            s = {"atoms": a}
            self.__structures.append(s)
            for k in meta:
                s[k] = sub[k]
