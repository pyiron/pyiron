# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron_atomistics.job.parallel import AtomisticParallelMaster
import numpy as np
try:
    import pylab as plt
except (ImportError, RuntimeError):
    try:
        import matplotlib.pyplot as plt
    except (ImportError, RuntimeError):
        pass

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


# ToDo: not all abstract methods implemented
class ConvergenceKpointParallel(AtomisticParallelMaster):
    def __init__(self, project, job_name='encut_conv'):
        """

        :param project:
        :param job_name:
        :return:
        """
        super(ConvergenceKpointParallel, self).__init__(project, job_name)
        self.__name__ = 'ConvergenceKpointParallel'
        self.__version__ = '0.0.1'

        # define default input
        self.input['steps'] = (2, 'increase of kpoints')
        self.input['min'] = (4, 'Kpoint Minimum')
        self.input['max'] = (8, 'Kpoint Maximum')

    def write_input(self):
        self.input['num_points'] = len(range(self.input['min'], self.input['max']+self.input['steps'], self.input['steps']))
        super(ConvergenceKpointParallel, self).write_input()

    def create_jobs(self):
        job_lst = []

        self.submission_status.submitted_jobs = 0
        for kpoint in range(self.input['min'], self.input['max']+self.input['steps'], self.input['steps']):
            if self.job_id and self.project.db.get_item_by_id(self.job_id)['status'] not in ['finished', 'aborted']:
                self._logger.debug("KpointConvergence child project {}".format(self.project.__str__()))
                ham = self._create_child_job("kpoint_mesh_" + str(kpoint))
                if ham.server.run_mode.non_modal and self.get_child_cores() + ham.server.cores > self.server.cores:
                    break
                self._logger.debug('create job: %s %s', ham.job_info_str, ham.master_id)
                ham.set_kpoints(mesh=[kpoint, kpoint, kpoint])
                self._logger.info('KpointConvergence: run job {}'.format(ham.job_name))
                self.submission_status.submit_next()
                if not ham.status.finished:
                    ham.run()
                self._logger.info('KpointConvergence: finished job {}'.format(ham.job_name))
                if ham.server.run_mode.thread:
                    job_lst.append(ham._process)
            else:
                self.refresh_job_status()
        process_lst = [process.communicate() for process in job_lst if process]
        self.status.suspended = True

    def collect_output(self):
        eng_lst, kpoint_lst = [], []
        for job_id in self.child_ids:  # add iter_jobs (should behave like a project)
            ham = self.project_hdf5.inspect(job_id)
            print('job_id: ', job_id, ham.status)
            eng_lst.append(ham["output/generic/energy_tot"][-1])
            kpoint_lst.append(int(ham.job_name.split('_')[-1]))
        arg_lst = np.argsort(kpoint_lst)
        self._output["energy"] = eng_lst[arg_lst]
        self._output["kpoints"] = kpoint_lst[arg_lst]

        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val

    def plot(self, plt_show=True):
        df = self.output_to_pandas()
        kpoint_lst, eng_lst = df["kpoints"], df["energy"]
        plt.plot(kpoint_lst, eng_lst, 'x-', markersize=20)
        plt.xlabel("Kpoint Mesh")
        plt.ylabel("energy (eV)")
        if plt_show:
            plt.show()
