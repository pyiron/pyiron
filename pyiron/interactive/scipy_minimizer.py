# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.atomistics.job.interactivewrapper import InteractiveWrapper
from pyiron_base import InputList
from pyiron.atomistics.job.interactive import GenericInteractiveOutput
from scipy.optimize import minimize
import scipy
import warnings

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

GPa_to_eV_by_A3 = (
    1e21 / scipy.constants.physical_constants["joule-electron volt relationship"][0]
)


class ScipyMinimizer(InteractiveWrapper):
    """
    Structure optimization class based on Scipy minimizers.

    Example I:

    # Position optimization
    >>> pr = Project('position')
    >>> job = pr.create_job('SomeAtomisticJob', 'atomistic')
    >>> job.structure = pr.create_structure('Al', 'fcc', 4.)
    >>> # it works also in the static mode, but interactive is recommended
    >>> job.server.run_mode.interactive = True
    >>> minim = pr.create_job('ScipyMinimizer', 'scipy')
    >>> minim.ref_job = job
    >>> minim.run()

    Example II:

    # Volume optimization:
    >>> pr = Project('volume')
    >>> job = pr.create_job('SomeAtomisticJob', 'atomistic')
    >>> job.structure = pr.create_structure('Al', 'fcc', 4.)
    >>> # it works also in the static mode, but interactive is recommended
    >>> job.server.run_mode.interactive = True
    >>> minim = pr.create_job('ScipyMinimizer', 'scipy')
    >>> minim.ref_job = job
    >>> minim.calc_minimize(pressure=0, volume_only=True)
    >>> minim.run()

    By setting `volume_only`, positions are not updated, so that only the
    pressures are minimized.

    It is possible to optimize both the volume and the positions, but since
    changing the cell also changes the definition of coordinates, it is a
    mathematically ill-defined problem and therefore it might end up in a
    physically undefined state. For this reason, it is strongly recommended to
    launch several jobs using the Murnaghan class (cf. `Murnaghan`).

    In order to perform volume optimization, the child job must have
    3x3-pressure output.
    """
    def __init__(self, project, job_name):
        super(ScipyMinimizer, self).__init__(project, job_name)
        self.__name__ = "ScipyMinimizer"
        self._ref_job = None
        self.input = Input()
        self.output = ScipyMinimizerOutput(job=self)
        self.interactive_cache = {}
        self._delete_existing_job = True

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it
        has to be implement in the individual classes.
        """
        super(ScipyMinimizer, self).set_input_to_read_only()
        self.input.read_only = True

    def write_input(self):
        pass

    def run_static(self):
        self.ref_job_initialize()
        self._logger.debug("cg status: " + str(self.status))
        if self.ref_job.server.run_mode.interactive:
            self._delete_existing_job = False
        self.ref_job.run(delete_existing_job=self._delete_existing_job)
        self.status.running = True
        if self.input.pressure is not None:
            x0 = self.ref_job.structure.cell.flatten()
            if not self.input.volume_only:
                x0 = np.append(x0, self.ref_job.structure.get_scaled_positions().flatten())
        else:
            x0 = self.ref_job.structure.positions.flatten()
        self.output._result = minimize(
            method=self.input.minimizer,
            fun=self._get_value,
            x0=x0,
            jac=self._get_gradient,
            tol=1.0e-10,
            options={'maxiter': self.input.ionic_steps,
                     'return_all': True }
        )
        self.status.collect = True
        self.collect_output()
        if self.ref_job.server.run_mode.interactive:
            self.ref_job.interactive_close()
        if self["output/convergence"] > 0:
            self.status.finished = True
        else:
            self.status.not_converged = True

    def _update(self, x):
        rerun = False
        if self.input.pressure is not None:
            if not np.allclose(x[:9], self.ref_job.structure.cell.flatten()):
                self.ref_job.structure.set_cell(x[:9].reshape(-3,3), scale_atoms=True)
                rerun = True
            if (not self.input.volume_only
                and not np.allclose(
                    x[9:], self.ref_job.structure.get_scaled_positions().flatten()
                )):
                self.ref_job.structure.set_scaled_positions(x[9:].reshape(-1, 3))
                rerun = True
        elif not np.allclose(x, self.ref_job.structure.positions.flatten()):
            self.ref_job.structure.positions = x.reshape(-1, 3)
            rerun = True
        if rerun:
            self.ref_job.run(delete_existing_job=self._delete_existing_job)

    def check_convergence(self):
        if self.input.ionic_energy_tolerance > 0:
            if len(self.ref_job.output.energy_pot) < 2:
                return False
            elif np.absolute(np.diff(self.ref_job.output.energy_pot)[-1]) > self.input.ionic_energy_tolerance:
                return False
            if self.input.ionic_force_tolerance==0:
                return True
        max_force = np.linalg.norm(self.ref_job.output.forces[-1], axis=-1).max()
        if self.input.pressure is None:
            if max_force > self.input.ionic_force_tolerance:
                return False
        elif self.input.volume_only:
            if np.absolute(self.ref_job.output.pressures[-1]-self.input.pressure).max() > self.input.pressure_tolerance:
                return False
        else:
            if max_force > self.input.ionic_force_tolerance:
                return False
            if np.absolute(self.ref_job.output.pressures[-1]-self.input.pressure).max() > self.input.pressure_tolerance:
                return False
        return True

    def _get_gradient(self, x):
        self._update(x)
        prefactor = 1
        if self.check_convergence():
            prefactor = 0
        if self.input.pressure is not None:
            pressure = -(
                self.ref_job.output.pressures[-1].flatten()-self.input.pressure.flatten()
            )
            if self.input.volume_only:
                return pressure*prefactor*GPa_to_eV_by_A3
            else:
                return np.append(
                    pressure*GPa_to_eV_by_A3,
                    -np.einsum(
                        'ij,ni->nj',
                        np.linalg.inv(self.ref_job.structure.cell),
                        self.ref_job.output.forces[-1]
                    ).flatten()
                ).flatten()*prefactor
        else:
            return -self.ref_job.output.forces[-1].flatten()*prefactor

    def _get_value(self, x):
        self._update(x)
        return self.ref_job.output.energy_pot[-1]

    def collect_output(self):
        self.output.to_hdf(self._hdf5)

    def to_hdf(self, hdf=None, group_name=None):
        super(ScipyMinimizer, self).to_hdf(hdf=hdf, group_name=group_name)
        self.output.to_hdf(self._hdf5)

    def calc_minimize(
        self,
        max_iter=100,
        pressure=None,
        algorithm='CG',
        ionic_energy_tolerance=0,
        ionic_force_tolerance=1.0e-2,
        pressure_tolerance=1.0e-3,
        volume_only=False,
    ):
        """
        Args:
            algorithm (str): scipy algorithm (currently only 'CG' and 'BFGS' run reliably)
            pressure (float/list/numpy.ndarray): target pressures
            max_iter (int): maximum number of iterations
            ionic_energy_tolerance (float): convergence goal in terms of
                                  energy (optional)
            ionic_force_tolerance (float): convergence goal in terms of
                                  forces (optional)
            volume_only (bool): Only pressure minimization
        """
        if pressure is None and volume_only:
            raise ValueError('pressure must be specified if volume_only')
        if pressure is not None and not volume_only:
            warnings.warn(
                'Simultaneous optimization of pressures and positions is a'
                + ' mathematically ill posed problem - there is no guarantee'
                + ' that it converges to the desired structure'
            )
        if pressure is not None and not hasattr(pressure, '__len__'):
            pressure = pressure*np.eye(3)
        self.input.minimizer = algorithm
        self.input.ionic_steps = max_iter
        self.input.pressure = pressure
        self.input.volume_only = volume_only
        self.input.ionic_force_tolerance = ionic_force_tolerance
        self.input.ionic_energy_tolerance = ionic_energy_tolerance
        self.input.pressure_tolerance = pressure_tolerance


class Input(InputList):
    """
    Args:
        minimizer (str): minimizer to use (currently only 'CG' and 'BFGS' run
            reliably)
        ionic_steps (int): max number of steps
        ionic_force_tolerance (float): maximum force tolerance
    """

    def __init__(self, input_file_name=None, table_name="input"):
        self.minimizer = 'CG'
        self.ionic_steps = 100
        self.ionic_force_tolerance = 1.0e-2
        self.pressure = None
        self.volume_only = False
        self.ionic_energy_tolerance = 0
        self.pressure_tolerance = 1.0e-3


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

