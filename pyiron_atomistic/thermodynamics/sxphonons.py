# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import os
import numpy as np
import scipy.constants
from pyiron_atomistic.atomistics.structure.atoms import Atoms
from pyiron_atomistic.sphinx.base import InputWriter, Output
from pyiron_base import GenericJob, GenericParameters, JobGenerator
from pyiron_atomistic.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistic.atomistics.master.parallel import AtomisticParallelMaster

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


BOHR_TO_ANGSTROM = (
    scipy.constants.physical_constants["Bohr radius"][0] / scipy.constants.angstrom
)
HARTREE_TO_EV = scipy.constants.physical_constants["Hartree energy in eV"][0]
HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM = HARTREE_TO_EV / BOHR_TO_ANGSTROM


class SxUniqDispl(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(SxUniqDispl, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "SxUniqDispl"
        self.input = GenericParameters(table_name="displacement")
        self.input["displacement"] = 0.01
        self.structure_lst = []
        self._id_pyi_to_spx = []
        self._id_spx_to_pyi = []

    @property
    def id_spx_to_pyi(self):
        if self.structure is None:
            return None
        if len(self._id_spx_to_pyi) == 0:
            self._initialize_order()
        return self._id_spx_to_pyi

    @property
    def id_pyi_to_spx(self):
        if self.structure is None:
            return None
        if len(self._id_pyi_to_spx) == 0:
            self._initialize_order()
        return self._id_pyi_to_spx

    def _initialize_order(self):
        for elm_species in self.structure.get_species_objects():
            self._id_pyi_to_spx.append(
                np.arange(len(self.structure))[
                    self.structure.get_chemical_symbols() == elm_species.Abbreviation
                ]
            )
        self._id_pyi_to_spx = np.array(
            [ooo for oo in self._id_pyi_to_spx for ooo in oo]
        )
        self._id_spx_to_pyi = np.array([0] * len(self._id_pyi_to_spx))
        for i, p in enumerate(self._id_pyi_to_spx):
            self._id_spx_to_pyi[p] = i

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(SxUniqDispl, self).set_input_to_read_only()
        self.input.read_only = True

    def list_structures(self):
        if self.status.finished:
            return self.structure_lst
        else:
            return []

    def write_structure(self, cwd, file_name="structure_wrapper.sx"):
        structure_file_name = "structure.sx"
        iw = InputWriter()
        iw.structure = self.structure
        iw.write_structure(file_name=structure_file_name, cwd=cwd)
        with open(os.path.join(cwd, file_name), "w") as f:
            f.writelines(["structure { include <" + structure_file_name + ">; }"])

    def extract_structure(self, working_directory):
        structure_lst = [self.structure]
        parser = Output(self)
        for f in os.listdir(working_directory):
            if "input-disp" in f:
                structure_template = self.structure.copy()
                parser.collect_relaxed_hist(file_name=f, cwd=working_directory)
                structure_template.cell = parser._parse_dict["cell"][0]
                structure_template.positions = parser._parse_dict["positions"][0]
                structure_lst.append(structure_template)
        return structure_lst

    def write_input(self):
        self.write_structure(
            cwd=self.working_directory, file_name="structure_wrapper.sx"
        )
        lines = [
            "#!/bin/bash\n",
            "sxuniqdispl --log -d "
            + str(float(self.input["displacement"]) / BOHR_TO_ANGSTROM)
            + " -i structure_wrapper.sx\n",
        ]
        with open(os.path.join(self.working_directory, "sxuniqdispl.sh"), "w") as f:
            f.writelines(lines)

    def collect_output(self):
        self.structure_lst = self.extract_structure(
            working_directory=self.working_directory
        )
        with self.project_hdf5.open("output") as hdf_out:
            for ind, struct in enumerate(self.structure_lst):
                struct.to_hdf(hdf=hdf_out, group_name="structure_" + str(ind))

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(SxUniqDispl, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(SxUniqDispl, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.from_hdf(hdf5_input)
        if "output" in self.project_hdf5.list_groups():
            with self.project_hdf5.open("output") as hdf5_output:
                self.structure_lst = [
                    Atoms().from_hdf(hdf5_output, group_name)
                    for group_name in hdf5_output.list_groups()
                ]


class SxDynMat(GenericJob):
    def __init__(self, project, job_name):
        super(SxDynMat, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "SxDynMat"
        self._child_lst = []
        self._child_id_lst = []

    @property
    def child_id_lst(self):
        return self._child_id_lst

    @child_id_lst.setter
    def child_id_lst(self, child_id_lst):
        self._child_id_lst = child_id_lst

    @property
    def child_lst(self):
        if len(self._child_lst) != len(self._child_id_lst):
            self._child_lst = [
                self.project.load(job_id) for job_id in self._child_id_lst
            ]
            forces_of_first_child = self._child_lst[0].output.forces[-1]
            self._forces_lst = [
                job.output.forces[-1] - forces_of_first_child for job in self._child_lst
            ]
            self._structure_lst = [job.get_structure() for job in self._child_lst]
        return self._child_lst

    @staticmethod
    def matrix_to_str(matrix):
        """
        Function to convert an numpy matrix to an Sphinx input compatible matrix.

        Args:
            matrix (numpy.d2type): the matrix to be converted

        Returns:
            str: the matrix representation in the Sphinx input.
        """
        output_str = "["
        for i in matrix:
            output_str += "["
            for j in i:
                output_str += str(j) + ", "
            output_str = output_str[:-2] + "], "
        output_str = output_str[:-2] + "]"
        return output_str

    @staticmethod
    def vector_to_str(vector):
        """
        Function to convert an numpy vector to an Sphinx input compatible vector.

        Args:
            vector (numpy.d2type): the vector to be converted

        Returns:
            str: the vector representation in the Sphinx input.
        """
        output_str = "["
        for i in vector:
            output_str += str(i) + ", "
        output_str = output_str[:-2] + "]"
        return output_str

    def write_sxdynmat(self, file_name="sxdynmat.sx", cwd=None):
        forces_lst = [job.output.forces[-1] for job in self.child_lst]
        structure_lst = [job.get_structure() for job in self.child_lst]
        phono_dat_str = "format phononDat;\n\n"

        initial_structure = structure_lst[0]
        phono_dat_str += "pseudoPot  {\n"
        for species_obj in initial_structure.get_species_objects():
            phono_dat_str += (
                "  species { reciprocalMass=1/" + str(species_obj.AtomicMass) + "; }\n"
            )
        phono_dat_str += "}\n\n"

        first_structure = True
        for structure, force_mat in zip(structure_lst, forces_lst):
            phono_dat_str += "structure  {\n"
            if first_structure:
                phono_dat_str += (
                    "   cell = "
                    + self.matrix_to_str(structure.cell * 1 / BOHR_TO_ANGSTROM)
                    + ";\n"
                )
                first_structure = False

            for species_obj in structure.get_species_objects():
                if species_obj.Parent:
                    species = species_obj.Parent
                else:
                    species = species_obj.Abbreviation
                phono_dat_str += "species {\n"
                for elm_pos, elm_species, elm_forces in zip(
                    structure.positions, structure.get_chemical_elements(), force_mat
                ):
                    if elm_species.Parent:
                        element = elm_species.Parent
                    else:
                        element = elm_species.Abbreviation
                    if element == species:
                        phono_dat_str += (
                            "    atom { coords = "
                            + self.vector_to_str(elm_pos * 1 / BOHR_TO_ANGSTROM)
                            + "; force = "
                            + self.vector_to_str(
                                elm_forces * 1 / HARTREE_OVER_BOHR_TO_EV_OVER_ANGSTROM
                            )
                            + "; }\n"
                        )
                phono_dat_str += "}\n"
            phono_dat_str += "}\n"
        if cwd is not None:
            file_name = os.path.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(phono_dat_str)

    def write_input(self):
        self.write_sxdynmat(cwd=self.working_directory, file_name="sxdynmat.sx")
        lines = ["#!/bin/bash\n", "sxdynmat --log -i sxdynmat.sx -H\n"]
        with open(os.path.join(self.working_directory, "sxuniqdispl.sh"), "w") as f:
            f.writelines(lines)

    def get_hesse_matrix(self):
        if "output" in self.project_hdf5.list_groups():
            return self.project_hdf5["output/hesse"]

    def collect_output(self):
        with self.project_hdf5.open("output") as hdf_out:
            hdf_out["hesse"] = np.loadtxt(
                os.path.join(self.working_directory, "HesseMatrix_sphinx")
            )

    def collect_logfiles(self):
        pass


class SxPhononsJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        return [
            ["supercell_phonon_%d" % ind, sc]
            for ind, sc in enumerate(self._master._displacement_job.structure_lst)
        ]

    @staticmethod
    def job_name(parameter):
        return parameter[0]

    def modify_job(self, job, parameter):
        job.structure = parameter[1]
        return job


class SxPhonons(AtomisticParallelMaster):
    def __init__(self, project, job_name):
        super(SxPhonons, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "SxPhonons"
        self.input["displacement"] = (0.01, "atoms displacement, Ang")
        self._displacement_job = None
        self._dynmat_job = None
        self._job_generator = SxPhononsJobGenerator(self)

    def run_static(self):
        if self._displacement_job is None:
            self._displacement_job = self.project.create_job(
                self.project.job_type.SxUniqDispl, self.job_name + "_dis"
            )
            self._displacement_job.input["displacement"] = self.input["displacement"]
            self._displacement_job.structure = self.structure
            self._displacement_job.run()
        super(SxPhonons, self).run_static()

    def collect_output(self):
        """

        Returns:

        """
        self._dynmat_job = self.project.create_job(
            self.project.job_type.SxDynMat, self.job_name + "_dyn"
        )
        self._dynmat_job.child_id_lst = self.child_ids
        self._dynmat_job.run()
        with self.project_hdf5.open("output") as hdf_out:
            hdf_out["hesse"] = self._dynmat_job.get_hesse_matrix()

    def get_hesse_matrix(self):
        if "output" in self.project_hdf5.list_groups():
            return self.project_hdf5["output/hesse"]


class SxHarmPotTst(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(SxHarmPotTst, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "SxHarmPotTst"
        self.input = GenericParameters(table_name="interaction")
        self.input["interaction_radius"] = 4.0
        self.input["maximum_noise"] = 0.26
        self._positions_lst = []
        self._forces_lst = []
        self._md_job_id = None
        self._md_job = None

    @property
    def md_job(self):
        if self._md_job is None and self._md_job_id is not None:
            self._md_job = self.project.load(self._md_job_id)
        return self._md_job

    @md_job.setter
    def md_job(self, job):
        if job.status == "finished":
            self._md_job_id = job.job_id
            self._md_job = job
            self._positions_lst = job["output/generic/positions"]
            self._forces_lst = job["output/generic/forces"]
        else:
            raise ValueError("Job not finished!")

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(SxHarmPotTst, self).set_input_to_read_only()
        self.input.read_only = True

    def write_harmpot(self, cwd, file_name="harmpot.sx"):
        harm_pot_str = (
            "format harmpot;\n\n"
            + "valenceCharge=0;\n"
            + "harmonicPotential {\n"
            + '   //include "refSym.sx";\n'
            + '   //include "equivalence.sx";\n'
            + "   maxDist="
            + str(float(self.input["maximum_noise"]) / BOHR_TO_ANGSTROM)
            + ";\n"
            + '   include "shells.sx";\n'
            + "}\n"
            + 'include "structure_wrapper.sx";'
        )
        if cwd is not None:
            file_name = os.path.join(cwd, file_name)
        with open(file_name, "w") as f:
            f.write(harm_pot_str)

    def write_structure(self, cwd, file_name="structure_wrapper.sx"):
        structure_file_name = "structure.sx"
        iw = InputWriter()
        iw.structure = self._md_job.structure
        iw.write_structure(file_name=structure_file_name, cwd=cwd)
        with open(os.path.join(cwd, file_name), "w") as f:
            f.writelines(["structure { include <" + structure_file_name + ">; }"])

    def validate_ready_to_run(self):
        if len(self._positions_lst) == 0 or len(self._forces_lst) == 0:
            raise ValueError()

    def write_input(self):
        self.write_structure(
            cwd=self.working_directory, file_name="structure_wrapper.sx"
        )
        self.write_harmpot(cwd=self.working_directory, file_name="harmpot.sx")
        pos_force_mat = np.concatenate((self._positions_lst, self._forces_lst), axis=2)
        cont_pos_force_mat = pos_force_mat.reshape(-1, pos_force_mat.shape[-1])
        np.savetxt(
            os.path.join(self.working_directory, "POSITIONs"), cont_pos_force_mat
        )
        lines = [
            "#!/bin/bash\n",
            "sxstructparam -i structure_wrapper.sx -c "
            + str(float(self.input["interaction_radius"]) / BOHR_TO_ANGSTROM)
            + " --printReduced=shells.sx --log\n",
            "sxharmpottst --param=POSITIONs --vasp --printHesse HesseMatrix_sphinx -i harmpot.sx --log --svd\n",
        ]
        with open(os.path.join(self.working_directory, "sxharmpottst.sh"), "w") as f:
            f.writelines(lines)

    def get_hesse_matrix(self):
        if "output" in self.project_hdf5.list_groups():
            return self.project_hdf5["output/hesse"]

    def collect_output(self):
        with self.project_hdf5.open("output") as hdf_out:
            hdf_out["hesse"] = np.loadtxt(
                os.path.join(self.working_directory, "HesseMatrix_sphinx")
            )

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(SxHarmPotTst, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.to_hdf(hdf5_input)
            if len(self._positions_lst) != 0:
                hdf5_input["positions"] = self._positions_lst
            if len(self._forces_lst) != 0:
                hdf5_input["forces"] = self._forces_lst
            if self._md_job_id is not None:
                hdf5_input["md_job_id"] = self._md_job_id

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(SxHarmPotTst, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.from_hdf(hdf5_input)
            if "positions" in hdf5_input.list_nodes():
                self._positions_lst = hdf5_input["positions"]
            if "forces" in hdf5_input.list_nodes():
                self._forces_lst = hdf5_input["forces"]
            if "md_job_id" in hdf5_input.list_nodes():
                self._md_job_id = hdf5_input["md_job_id"]
