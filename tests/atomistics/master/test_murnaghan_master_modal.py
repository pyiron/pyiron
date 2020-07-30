# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
from pyiron.atomistics.structure.atoms import CrystalStructure
from pyiron.base.project.generic import Project
import unittest


def convergence_goal(self, **qwargs):
    import numpy as np

    eps = 0.2
    if "eps" in qwargs:
        eps = qwargs["eps"]
    erg_lst = self.get_from_childs("output/generic/energy")
    var = 1000 * np.var(erg_lst)
    # print(var / len(erg_lst))
    if var / len(erg_lst) < eps:
        return True
    ham_prev = self[-1]
    job_name = self.first_child_name() + "_" + str(len(self))
    ham_next = ham_prev.restart(job_name=job_name)
    return ham_next


class TestMurnaghan(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(
            os.path.join(cls.file_location, "testing_murnaghan_master_modal")
        )
        cls.basis = CrystalStructure(
            element="Fe", bravais_basis="bcc", lattice_constant=2.8
        )
        cls.project.remove_jobs(recursive=True)
        # cls.project.remove_jobs(recursive=True)
        # self.project.set_logging_level('INFO')

    @classmethod
    def tearDownClass(cls):
        file_location = os.path.dirname(os.path.abspath(__file__))
        project = Project(os.path.join(file_location, "testing_murnaghan_master_modal"))
        project.remove_jobs(recursive=True)
        project.remove(enable=True, enforce=True)

    def test_run(self):
        # Even though the test is completed successful
        ham = self.project.create_job(
            self.project.job_type.AtomisticExampleJob, "job_test"
        )
        ham.structure = self.basis
        ham.server.run_mode.non_modal = True
        murn = self.project.create_job("Murnaghan", "murnaghan")
        murn.ref_job = ham
        murn.input["num_points"] = 3
        murn.server.run_mode.non_modal = True
        murn.run()
        self.assertFalse(ham.status.finished)
        self.project.wait_for_job(murn, interval_in_s=5, max_iterations=50)
        self.assertTrue(murn.status.finished)
        murn.remove()
        ham.remove()
        self.project.remove(enable=True, enforce=True)


if __name__ == "__main__":
    unittest.main()
