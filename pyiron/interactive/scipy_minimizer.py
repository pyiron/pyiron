# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.job.interactivewrapper import InteractiveWrapper
from pyiron_base import GenericParameters
from pyiron.atomistics.job.interactive import GenericInteractiveOutput
from scipy import optimize

__author__ = "Osamu Waseda"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Osamu Waseda"
__email__ = "waseda@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


class ScipyMinimizer(InteractiveWrapper):
    def __init__(self, project, job_name):
        super(ScipyMinimizer, self).__init__(project, job_name)
        self.__name__ = "ScipyMinimizer"
        self._ref_job = None
        self.input = Input()
        self.output = ScipyMinimizerOutput(job=self)
        self.interactive_cache = {}

    @property
    def minimizer(self):
        return self.input["minimizer"]

    @minimizer.setter
    def minimizer(self, minim):
        list_of_minimizers = ["CG", "BFGS", "simple"]
        if minim not in list_of_minimizers:
            raise ValueError("Minimizer has to be chosen from the following list:", " ".join(list_of_minimizers))
        else:
            self.input["minimizer"] = minim

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(ScipyMinimizer, self).set_input_to_read_only()
        self.input.read_only = True

    def write_input(self):
        pass

    def interactive_close(self):
        if self.interactive_is_activated():
            self._interactive_library.close()
            if len(self.interactive_cache[list(self.interactive_cache.keys())[0]]) != 0:
                self.interactive_flush(path="generic")
            super(ScipyMinimizer, self).interactive_close()

    def run_static(self):
        self.ref_job_initialize()
        self._logger.debug("cg status: " + str(self.status))
        self._delete_existing_job = True
        if self.ref_job.server.run_mode.interactive:
            self._delete_existing_job = False
        self.ref_job.run(delete_existing_job=self._delete_existing_job)
        self.status.running = True
        if self.input["minimizer"] == "CG":
            output = optimize.fmin_cg(
                f=self._update_energy,
                x0=self.ref_job.structure.positions.flatten(),
                fprime=self._update_forces,
                maxiter=self.input["ionic_steps"],
                gtol=self.input["ionic_force_tolerance"],
                disp=False,
                full_output=True,
            )
            self.output._convergence = output[4]
        elif self.input["minimizer"] == "BFGS":
            output = optimize.fmin_bfgs(
                f=self._update_energy,
                x0=self.ref_job.structure.positions.flatten(),
                fprime=self._update_forces,
                maxiter=self.input["ionic_steps"],
                gtol=self.input["ionic_force_tolerance"],
                disp=False,
                full_output=True,
            )
            self.output._hessian = output[3]
            self.output._convergence = output[6]
        elif self.input["minimizer"] == "simple":
            output = optimize.fmin(
                f=self._update_energy,
                x0=self.ref_job.structure.positions.flatten(),
                maxiter=self.input["ionic_steps"],
                gtol=self.input["ionic_force_tolerance"],
                disp=False,
                full_output=True,
            )
            self.output._hessian = output[4]
        self.status.collect = True
        self.collect_output()
        if self.ref_job.server.run_mode.interactive:
            self.ref_job.interactive_close()

    def _update_forces(self, x):
        x = np.array(x).reshape(-1, 3)
        self._logger.debug("cg ref_job status: " + str(self.ref_job.status))
        if not np.equal(x, self.ref_job.structure.positions).all():
            self.ref_job.structure.positions = x
            self.ref_job.run(delete_existing_job=self._delete_existing_job)
        f = self.ref_job.output.forces[-1].flatten()
        self._logger.debug("cg ref_job status after: " + str(self.ref_job.status))
        return -f

    def _update_energy(self, x):
        x = np.array(x).reshape(-1, 3)
        self._logger.debug("cg ref_job status: " + str(self.ref_job.status))
        if not np.equal(x, self.ref_job.structure.positions).all():
            self.ref_job.structure.positions = x
            self.ref_job.run(delete_existing_job=self._delete_existing_job)
        E = self.ref_job.output.energy_pot[-1]
        self._logger.debug("cg ref_job status after: " + str(self.ref_job.status))
        return E

    def collect_output(self):
        self.output.to_hdf(self._hdf5)

    def to_hdf(self, hdf=None, group_name=None):
        super(ScipyMinimizer, self).to_hdf(hdf=hdf, group_name=group_name)
        self.output.to_hdf(self._hdf5)


class Input(GenericParameters):
    """
    class to control the generic input for a Sphinx calculation.

    Args:
        input_file_name (str): name of the input file
        table_name (str): name of the GenericParameters table
    """

    def __init__(self, input_file_name=None, table_name="input"):
        super(Input, self).__init__(
            input_file_name=input_file_name,
            table_name=table_name,
            comment_char="//",
            separator_char="=",
            end_value_char=";",
        )

    def load_default(self):
        """
        Loads the default file content
        """
        file_content = (
            "minimizer = CG\n" "ionic_steps = 1000\n" "ionic_force_tolerance = 1.0e-8\n"
        )
        self.load_string(file_content)


class ScipyMinimizerOutput(GenericInteractiveOutput):
    def __init__(self, job):
        super(ScipyMinimizerOutput, self).__init__(job=job)
        self._convergence = None
        self._hessian = None

    def __dir__(self):
        return list(set(list(self._job.interactive_cache.keys())))

    def to_hdf(self, hdf, group_name="output"):
        with hdf.open(group_name) as hdf_output:
            hdf_output["convergence"] = self._convergence
            hdf_output["hessian"] = self._hessian

    def from_hdf(self, hdf, group_name="output"):
        if "convergence" in hdf[group_name].list_nodes():
            self._convergence = hdf[group_name]["convergence"]
        if "hessian" in hdf[group_name].list_nodes():
            self._hessian = hdf[group_name]["hessian"]
