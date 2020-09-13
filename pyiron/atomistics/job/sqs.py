# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from multiprocessing import cpu_count
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_base import GenericParameters
from pyiron.atomistics.structure.atoms import Atoms, ase_to_pyiron, pyiron_to_ase

try:
    from pymatgen.io.ase import AseAtomsAdaptor
    from sqsgenerator.core.sqs import ParallelSqsIterator
except ImportError:
    pass

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "0.1"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Aug 14, 2020"


def pyiron_to_pymatgen(structure):
    return AseAtomsAdaptor.get_structure(pyiron_to_ase(structure))


def pymatgen_to_pyiron(structure):
    return ase_to_pyiron(AseAtomsAdaptor.get_atoms(structure))


def get_sqs_structures(structure, mole_fractions, weights=None, objective=0.0, iterations=1e6, output_structures=10, num_threads=None):
    structure = pyiron_to_pymatgen(structure)
    if not weights:
        weights = {i: 1.0/i for i in range(1,7)}
    if not num_threads:
        # Thats default behaviour
        num_threads = cpu_count()
    iterator = ParallelSqsIterator(structure, mole_fractions, weights, num_threads=num_threads)
    structures, decmp, iter_, cycle_time = iterator.iteration(iterations=iterations, output_structures=output_structures, objective=objective)
    return [pymatgen_to_pyiron(s) for s in structures], decmp, iter_, cycle_time


class SQSJob(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(SQSJob, self).__init__(project, job_name)
        self.input = GenericParameters(table_name="input")
        self.input['mole_fractions'] = dict()
        self.input['weights'] = None
        self.input['objective'] = 0.0 
        self.input['iterations'] = 1e6 
        self.input['output_structures'] = 1
        self._python_only_job = True
        self._lst_of_struct = []
    
    @property
    def list_of_structures(self):
        return self._lst_of_struct

    def list_structures(self):
        if self.status.finished:
            return self._lst_of_struct
        else:
            return []
    
    # This function is executed 
    def run_static(self):
        self._lst_of_struct, decmp, iterations, cycle_time = get_sqs_structures(
            structure=self.structure, 
            mole_fractions=self.input['mole_fractions'], 
            weights=self.input['weights'], 
            objective=self.input['objective'], 
            iterations=self.input['iterations'], 
            output_structures=self.input['output_structures'], 
            num_threads=self.server.cores
        )
        for i, structure in enumerate(self._lst_of_struct):
            with self.project_hdf5.open("output/structures/structure_" + str(i)) as h5:
                structure.to_hdf(h5)
        with self.project_hdf5.open("output") as h5:
            h5["decmp"] = decmp
            h5["cycle_time"] = cycle_time
            h5["iterations"] = iterations
        self.status.finished = True
        self.project.db.item_update(self._runtime(), self.job_id)
        
    def to_hdf(self, hdf=None, group_name=None):
        super().to_hdf(
            hdf=hdf,
            group_name=group_name
        )
        with self.project_hdf5.open("input") as h5in:
            self.input.to_hdf(h5in)

    def from_hdf(self, hdf=None, group_name=None):
        super().from_hdf(
            hdf=hdf,
            group_name=group_name
        )
        with self.project_hdf5.open("input") as h5in:
            self.input.from_hdf(h5in)
        with self.project_hdf5.open("output/structures") as hdf5_output:
            structure_names = hdf5_output.list_groups()
        for group in structure_names:
            with self.project_hdf5.open("output/structures/" + group) as hdf5_output:
                self._lst_of_struct.append(Atoms().from_hdf(hdf5_output))
