# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from multiprocessing import cpu_count
from pyiron.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_base import InputList, GenericParameters, Settings
from pyiron.atomistics.structure.atoms import Atoms, ase_to_pyiron, pyiron_to_ase
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np

try:
    from sqsgenerator.core.sqs import ParallelSqsIterator
except ImportError:
    # sqsgenerator is not available for all systems
    # This is no reason to have pyiron fail catestrophically, it just means this job class is unavailable.
    # Pass silently for now, but throw an exception if the user tries to instantiate this job.
    ParallelSqsIterator = None

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


s = Settings()


def pyiron_to_pymatgen(structure):
    return AseAtomsAdaptor.get_structure(pyiron_to_ase(structure))


def pymatgen_to_pyiron(structure):
    return ase_to_pyiron(AseAtomsAdaptor.get_atoms(structure))


def get_sqs_structures(structure, mole_fractions, weights=None, objective=0.0, iterations=1e6, output_structures=10,
                       num_threads=None):
    structure = pyiron_to_pymatgen(structure)
    if not weights:
        weights = {i: 1.0/i for i in range(1, 7)}
    if not num_threads:
        # Thats default behaviour
        num_threads = cpu_count()
    iterator = ParallelSqsIterator(structure, mole_fractions, weights, num_threads=num_threads)
    structures, decmp, iter_, cycle_time = iterator.iteration(
        iterations=iterations, output_structures=output_structures, objective=objective
    )
    return [pymatgen_to_pyiron(s) for s in structures], decmp, iter_, cycle_time


def _fail_if_imports_missing():
    """
    Just a temporary measure as long as any imports are wrapped in a try/pass instead of being on the dependencies
    list.
    """
    if ParallelSqsIterator is None:
        raise ImportError("SQSJob relies on sqsgenerator.core.sqs.ParallelSqsIterator, but this is unavailable.")


class SQSJob(AtomisticGenericJob):
    """
    Produces special quasirandom structures designed to duplicate truly random chemical distributions as well as
    possible while using finite periodic cells.

    A pyiron wrapper for the [SQS code of Dominik Gehringer](https://github.com/dgehringer/sqsgenerator).

    'Structural models used in calculations of properties of substitutionally random $A_{1-x}B_x$ alloys are usually
    constructed by randomly occupying each of the N sites of a periodic cell by A or B. We show that it is possible to
    design ‘‘special quasirandom structures’’ (SQS’s) that mimic for small N (even N=8) the first few, physically most
    relevant radial correlation functions of a perfectly random structure far better than the standard technique does.
    We demonstrate the usefulness of these SQS’s by calculating optical and thermodynamic properties of a number of
    semiconductor alloys in the local-density formalism.'
    From the abstract of Zunger, Wei, Ferreira, and Bernard, Phys. Rev. Lett. 65 (1990) 353,
    DOI:https://doi.org/10.1103/PhysRevLett.65.353

    Input:
        - mole_fractions (dict): Approximate chemical composition for the output structure(s), using chemical symbols
            as the keys and floats as the values. Values should sum to 1, but within reason will be automatically
            adjusted to accommodate the number of atoms in the provided structure (adjustments printed to standard out).
            Vacancies can also be included using the key '0'.
        - weights (list/numpy.ndarray): Specifies the weights of the individual shells. (Default is None, which uses the
            inverse of the shell number for the weight, i.e.
            [1, 0.5, 0.3333333333, 0.25, 0.2, 0.166666667, 0.1428571429].)
        - objective (float): Specifies the value the objective functions. The program tries to reach the specified
            objective function. (Default is 0.)
        - iterations (int): How many iterations to make searching for the most special quasirandom structure using a
            random shuffling procedure. (Default is 1e6)
        - n_output_structures (int): How many different SQS structures to return (in decreasing special
            quasirandomness). (Default is 1.)

    Example:
        Case I: Get SQS for a given mole fraction:

        >>> job = SQSJob(job_name='sqs')
        >>> job.structure = structure
        >>> job.input.mole_fractions = {'Al': 0.8, 'Ni':0.2}
        >>> job.run()
        
        Case II: Get SQS for a given structure already containing at least 2 elements:

        >>> job = SQSJob(job_name='sqs')
        >>> job.structure = structure
        >>> job.run()

        In Case II, if the mole fractions will be overwritten if you specify the values (like in Case I)
    """

    publication = {
        "sqs": {
            "method": {
                "title": "Special quasirandom structures",
                "author": ["Zunger, A.", "Wei, S.-H.", "Ferreira, L.G.", "Bernard, J.E."],
                "journal": "Phys. Rev. Lett.",
                "volume": "65",
                "issue": "3",
                "pages": "353",
                "numpages": "0",
                "month": "July",
                "publisher": "American Physical Society",
                "doi": "10.1103/PhysRevLett.65.353",
                "url": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.65.353",
            }
        }
    }

    def __init__(self, project, job_name):
        super(SQSJob, self).__init__(project, job_name)
        # self.input = InputList(table_name='input')
        self.input = InputList(table_name='custom_dict')
        self.input.mole_fractions = dict()
        self.input.weights = None
        self.input.objective = 0.0
        self.input.iterations = 1e6
        self.input.n_output_structures = 1
        self._python_only_job = True
        self._lst_of_struct = []
        _fail_if_imports_missing()
        self.__hdf_version__ = "0.2.0"
        self.__name__ = "SQSJob"
        s.publication_add(self.publication)

    @property
    def list_of_structures(self):
        return self._lst_of_struct

    def validate_ready_to_run(self):
        super(SQSJob, self).validate_ready_to_run()
        if len(self.input.mole_fractions) == 0:
            chem = np.unique(self.structure.get_chemical_symbols(), return_counts=True)
            self.input.mole_fractions = dict(zip(chem[0], chem[1]/np.sum(chem[1])))
        if len(self.input.mole_fractions) == 1:
            raise ValueError('There must be at least two chemical elements')

    def db_entry(self):
        """
        Generate the initial database entry

        Returns:
            (dict): db_dict
        """
        db_dict = super(AtomisticGenericJob, self).db_entry()
        if self.structure:
            struct_len = len(self.structure)
            mole_fractions = {
                k: v for k, v in self.input.mole_fractions.items()
            }
            el_lst = sorted(mole_fractions.keys())
            atom_number_lst = [
                int(np.round(mole_fractions[el]*struct_len))
                for el in el_lst
            ]
            db_dict["ChemicalFormula"] = "".join([
                e + str(n)
                for e, n in zip(el_lst, atom_number_lst)
            ])
        return db_dict

    def list_structures(self):
        if self.status.finished:
            return self._lst_of_struct
        else:
            return []
    
    def get_structure(self, iteration_step=-1, wrap_atoms=True):
        """
        Gets the structure from a given iteration step of the simulation (MD/ionic relaxation). For static calculations
        there is only one ionic iteration step
        Args:
            iteration_step (int): Step for which the structure is requested
            wrap_atoms (bool): True if the atoms are to be wrapped back into the unit cell
        Returns:
            pyiron.atomistics.structure.atoms.Atoms: The required structure
        """
        return self.list_structures()[iteration_step]
    
    # This function is executed 
    def run_static(self):
        self._lst_of_struct, decmp, iterations, cycle_time = get_sqs_structures(
            structure=self.structure, 
            mole_fractions={
                k:v for k,v in self.input.mole_fractions.items()
            },
            weights=self.input.weights,
            objective=self.input.objective,
            iterations=self.input.iterations,
            output_structures=self.input.n_output_structures,
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
        self._structure_to_hdf()
        with self.project_hdf5.open("input") as h5in:
            self.input.to_hdf(h5in)

    def from_hdf(self, hdf=None, group_name=None):
        super().from_hdf(
            hdf=hdf,
            group_name=group_name
        )
        self._structure_from_hdf()
        self._backwards_compatible_input_from_hdf()
        with self.project_hdf5.open("output/structures") as hdf5_output:
            structure_names = hdf5_output.list_groups()
        for group in structure_names:
            with self.project_hdf5.open("output/structures/" + group) as hdf5_output:
                self._lst_of_struct.append(Atoms().from_hdf(hdf5_output))

    def _backwards_compatible_input_from_hdf(self):
        if "HDF_VERSION" in self.project_hdf5.list_nodes():
            version = self.project_hdf5["HDF_VERSION"]
        else:
            version = "0.1.0"

        if version == "0.1.0":
            with self.project_hdf5.open("input") as hdf5_input:
                gp = GenericParameters(table_name="custom_dict")
                gp.from_hdf(hdf5_input)
                for k in gp.keys():
                    self.input[k] = gp[k]
        elif version == "0.2.0":
            with self.project_hdf5.open("input") as hdf5_input:
                self.input.from_hdf(hdf5_input)
        else:
            raise ValueError("Cannot handle hdf version: {}".format(version))
