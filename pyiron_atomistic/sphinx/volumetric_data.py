# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import scipy
from scipy.io.netcdf import netcdf_file
import os
from pyiron_base import Settings
from pyiron_atomistic.atomistics.volumetric.generic import VolumetricData

__author__ = "Sudarsan Surendralal"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


BOHR_TO_ANGSTROM = (
    scipy.constants.physical_constants["Bohr radius"][0] /
    scipy.constants.angstrom
)


class SphinxVolumetricData(VolumetricData):
    """
    General class for parsing and manipulating volumetric static within Sphinx.

    The basic idea of the Base class is adapted from the pymatgen vasp
    VolumetricData class:

    http://pymatgen.org/_modules/pymatgen/io/vasp/outputs.html#VolumetricData

    """


    def __init__(self):
        super(SphinxVolumetricData, self).__init__()
        self.atoms = None
        self._diff_data = None
        self._total_data = None


    def from_file(self, filename, normalize=True):
        """
        Parses volumetric data from a sphinx binary (.sxb) file.

        Args:
            filename (str): Path of file to parse
            normalize (boolean): Flag to normalize by the volume of the cell
        """
        try:
            vol_data_list = self._read_vol_data(
                filename=filename,
                normalize=normalize
            )
        except (ValueError, IndexError, TypeError):
            raise ValueError("Unable to parse file: {}".format(filename))
        self._total_data = vol_data_list[0]
        if len(vol_data_list) == 2:
            self._total_data = vol_data_list[0] + vol_data_list[1]
            self._diff_data = vol_data_list[0] - vol_data_list[1]


    @staticmethod
    def _read_vol_data(filename, normalize=True):
        """
        Parses the Sphinx volumetric data files (rho.sxb and vElStat-eV.sxb).

        Args:
            filename (str): File to be parsed
            normalize (bool): Normalize the data with respect to the volume
                (probably sensible for rho)

        Returns:
            list: A list of the volumetric data (length >1 for density
                files with spin)

        """
        if not os.path.getsize(filename) > 0:
            s = Settings()
            s.logger.warning("File:" + filename + "seems to be empty! ")
            return None, None

        with netcdf_file(filename, mmap=False) as f:
            dim = [int(d) for d in f.variables["dim"]]
            volume = 1.0
            if normalize:
                cell = f.variables["cell"].data * BOHR_TO_ANGSTROM
                volume = np.abs(np.linalg.det(cell))
            if "mesh" in f.variables:
                # non-spin polarized
                total_data_list = [
                    np.array(f.variables["mesh"][:]).reshape(dim) / volume
                ]
            elif "mesh-0" in f.variables and "mesh-1" in f.variables:
                # spin-polarized
                total_data_list = [
                    np.array(f.variables["mesh-0"][:]).reshape(dim) / volume,
                    np.array(f.variables["mesh-1"][:]).reshape(dim) / volume
                ]
            else:
                raise ValueError(
                    "Unexpected keys in the netcdf file's variables: neither "
                    f"'mesh' nor 'mesh-0' and 'mesh-1' found in {f.variables}."
                )

        if len(total_data_list) == 0:
            s = Settings()
            s.logger.warning(
                "File:"
                + filename
                + "seems to be corrupted/empty even after parsing!"
            )
            return None

        return total_data_list


    @property
    def total_data(self):
        """
        numpy.ndarray: Total volumtric data (3D)
        """
        return self._total_data


    @total_data.setter
    def total_data(self, val):
        self._total_data = val


    @property
    def diff_data(self):
        """
        numpy.ndarray: Volumtric difference data (3D)
        """
        return self._diff_data


    @diff_data.setter
    def diff_data(self, val):
        self._diff_data = val

    def to_hdf(self, hdf, group_name="volumetric_data"):
        """
        Writes the data as a group to a HDF5 file

        Args:
            hdf (pyiron_base.generic.hdfio.ProjectHDFio): The
                HDF file/path to write the data
            group_name (str): The name of the group under which
                the data must be stored

        """
        with hdf.open(group_name) as hdf_vd:
            hdf_vd["TYPE"] = str(type(self))
            hdf_vd["total"] = self.total_data
            if self.diff_data is not None:
                hdf_vd["diff"] = self.diff_data

    def from_hdf(self, hdf, group_name="volumetric_data"):
        """
        Extract a VolumetricData instance from an HDF5 file.

        Args:
            hdf (pyiron_base.generic.hdfio.ProjectHDFio): The HDF
                file/path to read the data
            group_name (str): The name of the group under which
                the data have been stored

        Returns:
            pyiron.atomistics.volumetric.generic.VolumetricData: The
                VolumetricData instance

        """
        with hdf.open(group_name) as hdf_vd:
            self._total_data = hdf_vd["total"]
            if len(self._total_data) == 2:
                self.is_spin_polarized = True
            if "diff" in hdf_vd.list_nodes():
                self._diff_data = hdf_vd["diff"]
