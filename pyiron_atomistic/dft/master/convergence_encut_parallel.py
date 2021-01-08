# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron_atomistic.atomistics.master.parallel import AtomisticParallelMaster
from pyiron_base import JobGenerator
import numpy as np

try:
    import pylab as plt
except (ImportError, RuntimeError):
    try:
        import matplotlib.pyplot as plt
    except (ImportError, RuntimeError):
        pass

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


class EncutConvergenceJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        return [
            np.round(encut, 7)
            for encut in np.linspace(
                self._master.input["min"],
                self._master.input["max"],
                self._master.input["num_points"],
            )
        ]

    @staticmethod
    def job_name(parameter):
        return "encut_" + str(parameter).replace(".", "_")

    @staticmethod
    def modify_job(job, parameter):
        job.set_encut(encut=parameter)
        return job


# ToDo: not all abstract methods implemented
class ConvEncutParallel(AtomisticParallelMaster):
    def __init__(self, project, job_name="encut_conv"):
        """

        Args:
            project:
            job_name:
        """
        super(ConvEncutParallel, self).__init__(project, job_name)
        self.__name__ = "ConvEncutParallel"
        self.__version__ = "0.0.1"

        # define default input
        self.input["num_points"] = (11, "number of sample points")
        self.input["min"] = (200, "EnCut Minimum")
        self.input["max"] = (800, "EnCut Maximum")
        self._job_generator = EncutConvergenceJobGenerator(self)

    def collect_output(self):
        eng_lst, encut_lst = [], []
        for job_id in self.child_ids:  # add iter_jobs (should behave like a project)
            ham = self.project_hdf5.inspect(job_id)
            print("job_id: ", job_id, ham.status)
            eng_lst.append(ham["output/generic/energy_tot"][-1])
            encut = ham.job_name.split("_")[1:]
            encut_lst.append(float(encut[0] + "." + encut[1]))
        arg_lst = np.argsort(encut_lst)
        self._output["energy"] = eng_lst[arg_lst]
        self._output["encut"] = encut_lst[arg_lst]

        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val

    def plot(self, plt_show=True):
        df = self.output_to_pandas()
        encut_lst, eng_lst = df["encut"], df["energy"]
        plt.plot(encut_lst, eng_lst, "x-", markersize=20)
        plt.xlabel("EnCut (eV)")
        plt.ylabel("energy (eV)")
        if plt_show:
            plt.show()
