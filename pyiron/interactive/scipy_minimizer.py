# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.job.interactivewrapper import InteractiveWrapper
from pyiron_base import InputList
from pyiron.atomistics.job.interactive import GenericInteractiveOutput
from scipy.optimize import minimize

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
        self.output._result = minimize(
            method=self.input.minimizer,
            fun=self._update_energy,
            x0=self.ref_job.structure.positions.flatten(),
            jac=self._update_forces,
            tol=self.input.ionic_force_tolerance,
            options={'maxiter': self.input.ionic_steps,
                     'return_all': True }
        )
        self.status.collect = True
        self.collect_output()
        if self.ref_job.server.run_mode.interactive:
            self.ref_job.interactive_close()

    def _update(self, x):
        if not np.allclose(x, self.ref_job.structure.positions.flatten()):
            self.ref_job.structure.positions = x.reshape(-1, 3)
            self.ref_job.run(delete_existing_job=self._delete_existing_job)

    def _update_forces(self, x):
        self._update(x)
        return -self.ref_job.output.forces[-1].flatten()

    def _update_energy(self, x):
        self._update(x)
        return self.ref_job.output.energy_pot[-1]

    def collect_output(self):
        self.output.to_hdf(self._hdf5)

    def to_hdf(self, hdf=None, group_name=None):
        super(ScipyMinimizer, self).to_hdf(hdf=hdf, group_name=group_name)
        self.output.to_hdf(self._hdf5)


class Input(InputList):
    """
    Args:
        minimizer (str): minimizer to use (currently only 'CG' and 'BFGS' run reliably)
        ionic_steps (int): max number of steps
        ionic_force_tolerance (float): maximum force tolerance
    """

    def __init__(self, input_file_name=None, table_name="input"):
        self.minimizer = 'CG'
        self.ionic_steps = 100
        self.ionic_force_tolerance = 1.0e-2


class ScipyMinimizerOutput(GenericInteractiveOutput):
    def __init__(self, job):
        super(ScipyMinimizerOutput, self).__init__(job=job)
        self._result = None

    def to_hdf(self, hdf, group_name="output"):
        if self._result is None:
            return
        with hdf.open(group_name) as hdf_output:
            hdf_output["convergence"] = self._result['success']
            if 'hess_inv' in self._result.keys():
                hdf_output["hessian"] = self._result['hess_inv']

    def from_hdf(self, hdf, group_name="output"):
        if "convergence" in hdf[group_name].list_nodes():
            self._convergence = hdf[group_name]["convergence"]
        if "hessian" in hdf[group_name].list_nodes():
            self._hessian = hdf[group_name]["hessian"]
