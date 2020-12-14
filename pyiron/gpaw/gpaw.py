# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.gpaw.pyiron_ase import AseJob
from pyiron.dft.job.generic import GenericDFTJob
from pyiron_base import GenericParameters, Settings
from pyiron_base.generic.util import ImportAlarm

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

s = Settings()

try:
    from gpaw import GPAW as GPAWcode, PW, MethfesselPaxton
    import_alarm = ImportAlarm()
except ImportError:
    import_alarm = ImportAlarm(
        "Gpaw relies on the gpaw module but th is unavailable. Please ensure your python environment contains gpaw, "
        "e.g. by running `conda install -c conda-forge gpaw`."
    )


class Gpaw(AseJob, GenericDFTJob):
    @import_alarm
    def __init__(self, project, job_name):
        super(Gpaw, self).__init__(project, job_name)
        self.__name__ = "GpawJob"
        self.input = GpawInput()

    @property
    def plane_wave_cutoff(self):
        return self.input["encut"]

    @plane_wave_cutoff.setter
    def plane_wave_cutoff(self, val):
        self.input["encut"] = val

    def get_kpoints(self):
        return self.input["kpoints"]

    def _set_kpoints(
        self,
        mesh=None,
        scheme="MP",
        center_shift=None,
        symmetry_reduction=True,
        manual_kpoints=None,
        weights=None,
        reciprocal=True,
        n_path=None,
        path_name=None,
    ):
        if scheme != "MP":
            raise ValueError(
                "Currently only MP is supported in the pyiron wrapper."
            )
        if center_shift is not None:
            raise ValueError(
                "centershift is not implemented in the pyiron wrapper."
            )
        if not symmetry_reduction:
            raise ValueError(
                "symmetry_reduction is not implemented in the pyiron wrapper."
            )
        if manual_kpoints is not None:
            raise ValueError(
                "manual_kpoints are not implemented in the pyiron wrapper."
            )
        if weights is not None:
            raise ValueError(
                "weights are not implemented in the pyiron wrapper."
            )
        if not reciprocal:
            raise ValueError(
                "reciprocal is not implemented in the pyiron wrapper."
            )
        if n_path is not None:
            raise ValueError(
                "n_path is not implemented in the pyiron wrapper."
            )
        if path_name is not None:
            raise ValueError(
                "path_name is not implemented in the pyiron wrapper."
            )
        self.input["kpoints"] = mesh

    def set_calculator(self):
        kpoints = self.input["kpoints"]
        if isinstance(kpoints, str):
            kpoints = (
                self.input["kpoints"].replace("[", "").replace("]", "").split()
            )
        self._create_working_directory()
        calc = GPAWcode(
            mode=PW(float(self.input["encut"])),
            xc=self.input["potential"],
            occupations=MethfesselPaxton(width=float(self.input["sigma"])),
            kpts=kpoints,
            txt=self.working_directory + "/" + self.job_name + ".txt",
        )
        self.structure.set_calculator(calc)

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(Gpaw, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(Gpaw, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.from_hdf(hdf5_input)


class GpawInput(GenericParameters):
    """
    Input class for the ExampleJob based on the GenericParameters class.

    Args:
        input_file_name (str): Name of the input file - optional
    """

    def __init__(self, input_file_name=None):
        super(GpawInput, self).__init__(
            input_file_name=input_file_name, table_name="input", comment_char="#"
        )

    def load_default(self):
        """
        Loading the default settings for the input file.
        """
        input_str = """\
kpoints [6,6,6]
sigma 0.1
encut 350
potential PBE
"""
        self.load_string(input_str)
