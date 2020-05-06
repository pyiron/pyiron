# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import ast
import logging
import numpy as np

"""
Simple python executable that resembles the behavior of a real job
"""

__author__ = "Joerg Neugebauer"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


# Set the logging behaviour of the executable
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    datefmt="%m-%d %H:%M",
    filename="info.log",
    filemode="w",
)
module_logger = logging.getLogger("exampleExecutable")


class ExampleExecutable(object):
    """
    Simple python executable that resembles the behavior of a real job, i.e., it processes input files, generates
    output files and can perform restarts. Will be mainly used to test the genericJob class.
    """

    def __init__(self):
        print("Execution started")
        self.logger = logging.getLogger("exampleExecutable.module")
        self.logger.info("creating an instance")

        self.input_file = "input.inp"
        self.output_file = "output.log"
        self._count = 0
        self._alat_0 = None
        self._potential = None
        self._alpha = None
        self._alat = None
        self._count = None

    def write_restart(self):
        """
        Write a restart file for a continous simulation divided into multiple jobs.
        """
        with open("restart.out", "w") as f:
            f.write("count {0} \n".format(str(self._count)))

    def read_restart(self):
        """
        Read a restart file for a continous simulation divided into multiple jobs.
        """
        with open("restart.inp", "r") as f:
            line = f.readline()
            _, count = line.split()
            self._count = ast.literal_eval(count)

    def get_energy(self, alat):
        """
        Based on the lattice constant a random energy is calculated.

        Args:
            alat (float): lattice constant

        Returns:
            (list): list of n random energy values, where n equals self._count
        """
        return (
            self._potential(alat - self._alat_0)
            + np.random.random(self._count) * self._alpha
        )

    def run_lib(self, input_dict):
        """
        Run lib executes the job directly in Python instead of calling a separate subprocess, this is faster for all
        Python jobs but not available for non Python jobs. No input or output files are generated when running in
        library mode, instead all input is provided as an input dictionary and the output is returned as a list.

        Args:
            input_dict (dict): input consisting of ["alpha", "alat", "count"]

        Returns:
            list: alat(float), count(int), energy(list)
        """
        n_max = 4  # max. order of polynomial describing the potential
        pot_lst = [float(input_dict["a_" + str(i)]) for i in range(n_max, -1, -1)]
        self._alat_0 = pot_lst[n_max]
        pot_lst[4] = 0
        self._potential = np.poly1d(pot_lst)
        self._alpha = float(input_dict["alpha"])
        self._alat = float(input_dict["alat"])
        self._count = int(input_dict["count"])

        # make the executable a bit more real and throw warnings and error messages
        # depending on the input parameters
        self.logger.debug("type: alpha %s", type(input_dict["alpha"]))
        if self._alpha < 0:
            raise ValueError("noise amplitude alpha < 0")
        if self._count < 1:
            raise ValueError("number of energy steps must be larger than 0")

        if self._alat < 1:
            self.logger.warning("lattice constant alat < 1")

        self.logger.info("Execute program")
        energy = self.get_energy(alat=self._alat)
        return self._alat, self._count, energy

    def run(self):
        """
        Run executes the job as part of a subprocess. The input is written to an input file, then the executable is
        executed and finally the output file is collected and converted to HDF5 format for future processing.
        """
        with open(self.input_file, mode="r") as f:
            input_dict = {}
            for line in f.readlines():
                line = line.strip().split()
                key, value = line[0], line[1]
                # key, value = line.split()
                input_dict[key] = value
                self.logger.info("-> %s %s", key, str(value))

        # parse the input into the correct format
        n_max = 4  # max. order of polynomial describing the potential
        pot_lst = [float(input_dict["a_" + str(i)]) for i in range(n_max, -1, -1)]
        self._alat_0 = pot_lst[n_max]
        pot_lst[4] = 0
        self._potential = np.poly1d(pot_lst)
        self._alpha = float(input_dict["alpha"])
        self._alat = float(input_dict["alat"])
        self._count = int(input_dict["count"])

        # make the executable a bit more real and throw warnings and error messages
        # depending on the input parameters
        self.logger.debug("type: alpha %s", type(input_dict["alpha"]))
        if self._alpha < 0:
            raise ValueError("noise amplitude alpha < 0")
        if self._count < 1:
            raise ValueError("number of energy steps must be larger than 0")

        if self._alat < 1:
            self.logger.warning("lattice constant alat < 1")

        # all values in input_dict are str
        # only the actual code knows which type is needed, i.e., the conversion is done here
        if ast.literal_eval(input_dict["read_restart"]):
            self.logger.info("read restart file")
            self.read_restart()
            self.logger.info("restart file has been successfully read")

        self.logger.info("Execute program")
        energy = self.get_energy(alat=self._alat)

        self.logger.info("Program has been successfully terminated")
        if ast.literal_eval(input_dict["write_restart"]):
            self.logger.info("Write restart file")
            self.write_restart()

        # TODO: use hdf5 output file
        with open(self.output_file, mode="w") as f:
            self.logger.info("Write output file")
            f.write("exampleExecutable logFile \n")
            f.write("alat {0} \n".format(str(self._alat)))
            f.write("count {0} \n".format(str(self._count)))
            for e in energy:
                f.write("energy  {0} \n".format(str(e)))


if __name__ == "__main__":
    ExampleExecutable().run()
