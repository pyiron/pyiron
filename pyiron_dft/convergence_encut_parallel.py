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
class ConvergenceEncutParallel(AtomisticParallelMaster):
    def __init__(self, project, job_name='encut_conv'):
        """

        :param project:
        :param job_name:
        :return:
        """
        super(ConvergenceEncutParallel, self).__init__(project, job_name)
        self.__name__ = 'ConvergenceEncutParallel'
        self.__version__ = '0.0.1'

        # define default input
        self.input['num_points'] = (11, 'number of sample points')
        self.input['min'] = (200, 'EnCut Minimum')
        self.input['max'] = (800, 'EnCut Maximum')

    def create_jobs(self):
        job_lst = []

        self.submission_status.submitted_jobs = 0
        for encut in np.linspace(self.input['min'], self.input['max'], self.input['num_points']):
            if self.job_id and self.project.db.get_item_by_id(self.job_id)['status'] not in ['finished', 'aborted']:
                self._logger.debug("EnCutConvergence child project {}".format(self.project.__str__()))
                ham = self._create_child_job("encut_" + str(encut).replace('.', '_'))
                if ham.server.run_mode.non_modal and self.get_child_cores() + ham.server.cores > self.server.cores:
                    break
                self._logger.debug('create job: %s %s', ham.job_info_str, ham.master_id)
                ham.set_encut(encut=encut)
                self._logger.info('EnCutConvergence: run job {}'.format(ham.job_name))
                self.submission_status.submit_next()
                if not ham.status.finished:
                    ham.run()
                self._logger.info('EnCutConvergence: finished job {}'.format(ham.job_name))
                if ham.server.run_mode.thread:
                    job_lst.append(ham._process)
            else:
                self.refresh_job_status()
        process_lst = [process.communicate() for process in job_lst if process]
        self.status.suspended = True

    def collect_output(self):
        eng_lst, encut_lst = [], []
        for job_id in self.child_ids:  # add iter_jobs (should behave like a project)
            ham = self.project_hdf5.inspect(job_id)
            print('job_id: ', job_id, ham.status)
            eng_lst.append(ham["output/generic/energy_tot"][-1])
            encut = ham.job_name.split('_')[1:]
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
        plt.plot(encut_lst, eng_lst, 'x-', markersize=20)
        plt.xlabel("EnCut (eV)")
        plt.ylabel("energy (eV)")
        if plt_show:
            plt.show()
