# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
import os
from pyiron_atomistic.project import Project


def convergence_goal(self, **qwargs):
    import numpy as np

    eps = 0.2
    if "eps" in qwargs:
        eps = qwargs["eps"]
    erg_lst = self.get_from_childs("output/generic/energy")
    var = 1000 * np.var(erg_lst)
    print(var / len(erg_lst))
    if var / len(erg_lst) < eps:
        return True
    ham_prev = self[-1]
    job_name = self.first_child_name() + "_" + str(len(self))
    ham_next = ham_prev.restart(job_name=job_name)
    return ham_next


class TestSerialMaster(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(
            os.path.join(cls.file_location, "testing_serial_non_modal")
        )
        cls.project.remove_jobs_silently(recursive=True)

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, "testing_serial_non_modal"))
        project.remove(enable=True, enforce=True)

    def test_single_job_non_modal(self):
        ham = self.project.create_job(self.project.job_type.ExampleJob, "job_single_nm")
        ham.server.run_mode.non_modal = True
        job_ser = self.project.create_job(
            self.project.job_type.SerialMasterBase, "sequence_single_nm"
        )
        job_ser.append(ham)
        job_ser.run()
        job_ser.from_hdf()
        self.assertTrue(job_ser.status.finished)
        self.assertTrue(job_ser[0].status.finished)
        self.assertEqual(len(job_ser), 1)
        job_ser.remove()

    def test_single_job_non_modal_all(self):
        ham = self.project.create_job(
            self.project.job_type.ExampleJob, "job_single_nma"
        )
        ham.server.run_mode.non_modal = True
        job_ser = self.project.create_job(
            self.project.job_type.SerialMasterBase, "sequence_single_nma"
        )
        job_ser.server.run_mode.non_modal = True
        job_ser.append(ham)
        job_ser.run()
        self.assertFalse(job_ser.status.finished)
        self.project.wait_for_job(job_ser, interval_in_s=5, max_iterations=50)
        job_ser_reload = self.project.load(job_ser.job_id)
        self.assertTrue(job_ser_reload.status.finished)
        self.assertTrue(job_ser_reload[0].status.finished)
        self.assertEqual(len(job_ser_reload), 1)
        job_ser_reload.remove()

    def test_convergence_goal(self):
        # self.project.set_logging_level("DEBUG")
        ham = self.project.create_job(
            self.project.job_type.ExampleJob, "job_convergence"
        )
        ham.server.run_mode.non_modal = True
        # print (self.project.job_table())
        job_ser = self.project.create_job(
            self.project.job_type.SerialMasterBase, "sequence_convergence"
        )
        job_ser.server.run_mode.non_modal = True
        job_ser.append(ham)
        job_ser.set_goal(convergence_goal=convergence_goal, eps=0.2)
        job_ser.run()
        self.assertFalse(job_ser.status.finished)
        self.project.wait_for_job(job_ser, interval_in_s=5, max_iterations=50)
        self.assertTrue(job_ser.status.finished)
        self.assertTrue(len(job_ser) > 0)
        job_ser.remove()


if __name__ == "__main__":
    unittest.main()
