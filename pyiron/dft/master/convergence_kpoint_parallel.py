# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron.atomistics.master.parallel import AtomisticParallelMaster
from pyiron_base.master.parallel import JobGenerator
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


class KpointConvergenceJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        return [
            kpoint
            for kpoint in range(
                self.input["min"],
                self.input["max"] + self.input["steps"],
                self.input["steps"],
            )
        ]

    @staticmethod
    def job_name(parameter):
        return "kpoint_mesh_" + str(parameter)

    @staticmethod
    def modify_job(job, parameter):
        job.set_kpoints(mesh=[parameter, parameter, parameter])
        return job


# ToDo: not all abstract methods implemented
class ConvKpointParallel(AtomisticParallelMaster):
    def __init__(self, project, job_name="encut_conv"):
        """

        Args:
            project:
            job_name:
        """
        super(ConvKpointParallel, self).__init__(project, job_name)
        self.__name__ = "ConvKpointParallel"
        self.__version__ = "0.0.1"

        # define default input
        self.input["steps"] = (2, "increase of kpoints")
        self.input["min"] = (4, "Kpoint Minimum")
        self.input["max"] = (8, "Kpoint Maximum")
        self._job_generator = KpointConvergenceJobGenerator(self)

    def write_input(self):
        self.input["num_points"] = len(
            range(
                self.input["min"],
                self.input["max"] + self.input["steps"],
                self.input["steps"],
            )
        )
        super(ConvKpointParallel, self).write_input()

    def collect_output(self):
        eng_lst, kpoint_lst = [], []
        for job_id in self.child_ids:  # add iter_jobs (should behave like a project)
            ham = self.project_hdf5.inspect(job_id)
            print("job_id: ", job_id, ham.status)
            eng_lst.append(ham["output/generic/energy_tot"][-1])
            kpoint_lst.append(int(ham.job_name.split("_")[-1]))
        arg_lst = np.argsort(kpoint_lst)
        self._output["energy"] = eng_lst[arg_lst]
        self._output["kpoints"] = kpoint_lst[arg_lst]

        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val

    def plot(self, plt_show=True):
        df = self.output_to_pandas()
        kpoint_lst, eng_lst = df["kpoints"], df["energy"]
        plt.plot(kpoint_lst, eng_lst, "x-", markersize=20)
        plt.xlabel("Kpoint Mesh")
        plt.ylabel("energy (eV)")
        if plt_show:
            plt.show()
