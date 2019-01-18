# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import math

import numpy as np

from pyiron.base.settings.generic import Settings
from pyiron.vasp.structure import atoms_from_string, get_species_list_from_potcar
from pyiron.atomistics.volumetric.generic import VolumetricData

__author__ = "Sudarsan Surendralal"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class VaspVolumetricData(VolumetricData):
    """
    General class for parsing and manipulating volumetric static within VASP. The basic idea of the Base class is
    adapted from the pymatgen vasp VolumtricData class

    http://pymatgen.org/_modules/pymatgen/io/vasp/outputs.html#VolumetricData

    """

    def __init__(self):
        super(VaspVolumetricData, self).__init__()
        self.atoms = None
        self._diff_data = None

    def from_file(self, filename, normalize=True):
        """
        Convenience method to parse a generic volumetric static file in the vasp like format.
        Used by subclasses for parsing the file. This routine is adapted from the pymatgen vasp VolumetricData
        class with very minor modifications

        http://pymatgen.org/_modules/pymatgen/io/vasp/outputs.html#VolumetricData.

        Args:
            filename (str): Path of file to parse
            normalize (boolean): Flag to normalize by the volume of the cell

        """
        poscar_read = False
        poscar_string = list()
        dataset = list()
        all_dataset = list()
        dim = None
        dimline = None
        read_dataset = False
        ngrid_pts = 0
        data_count = 0
        atoms = None
        volume = None
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if read_dataset:
                    toks = line.split()
                    for tok in toks:
                        if data_count < ngrid_pts:
                            # This complicated procedure is necessary because
                            # vasp outputs x as the fastest index, followed by y
                            # then z.
                            x = data_count % dim[0]
                            y = int(math.floor(data_count / dim[0])) % dim[1]
                            z = int(math.floor(data_count / dim[0] / dim[1]))
                            dataset[x, y, z] = float(tok)
                            data_count += 1
                    if data_count >= ngrid_pts:
                        read_dataset = False
                        data_count = 0
                        all_dataset.append(dataset)
                elif not poscar_read:
                    if line != "" or len(poscar_string) == 0:
                        poscar_string.append(line)
                    elif line == "":
                        try:
                            atoms = atoms_from_string(poscar_string)
                        except ValueError:
                            pot_str = filename.split("/")
                            pot_str[-1] = "POTCAR"
                            potcar_file = "/".join(pot_str)
                            species = get_species_list_from_potcar(potcar_file)
                            atoms = atoms_from_string(poscar_string, species_list=species)
                        volume = atoms.get_volume()
                        poscar_read = True
                elif not dim:
                    dim = [int(i) for i in line.split()]
                    ngrid_pts = dim[0] * dim[1] * dim[2]
                    dimline = line
                    read_dataset = True
                    dataset = np.zeros(dim)
                elif line == dimline:
                    read_dataset = True
                    dataset = np.zeros(dim)
            if not normalize:
                volume = 1.0
            if len(all_dataset) == 0:
                s = Settings()
                s.logger.warning("File:" + filename + "seems to be corrupted/empty")
                return
            if len(all_dataset) == 2:
                data = {"total": all_dataset[0] / volume, "diff": all_dataset[1] / volume}
                self.diff_data = data["diff"]
            else:
                data = {"total": all_dataset[0] / volume}
            self.atoms = atoms
            self.total_data = data["total"]

    @property
    def diff_data(self):
        """
        numpy.ndarray: Volumtric difference data (3D)
        """
        return self._diff_data

    @diff_data.setter
    def diff_data(self, val):
        self._diff_data = val

    def to_hdf(self, hdf5, group_name="volumetric_data"):
        """
        Writes the data as a group to a HDF5 file

        Args:
            hdf5 (pyiron.base.generic.hdfio.ProjectHDFio): The HDF file/path to write the data to
            group_name (str): The name of the group under which the data must be stored as

        """
        with hdf5.open(group_name) as hdf_vd:
            hdf_vd["TYPE"] = str(type(self))
            hdf_vd["total"] = self.total_data
            if self.diff_data is not None:
                hdf_vd["diff"] = self.diff_data

    def from_hdf(self, hdf5, group_name="volumetric_data"):
        """
        Recreating the VolumetricData instance by reading data from the HDF5 files

        Args:
            hdf5 (pyiron.base.generic.hdfio.ProjectHDFio): The HDF file/path to write the data to
            group_name (str): The name of the group under which the data must be stored as

        Returns:
            pyiron.atomistics.volumetric.generic.VolumetricData: The VolumetricData instance

        """
        with hdf5.open(group_name) as hdf_vd:
            self.total_data = hdf_vd["total"]
            if "diff" in hdf_vd.list_nodes():
                self.diff_data = hdf_vd["diff"]
