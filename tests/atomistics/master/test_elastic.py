# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import unittest
from pyiron.atomistics.structure.atoms import CrystalStructure
from pyiron_base import Project
from pyiron.atomistics.master.elastic import calc_elastic_tensor, calc_elastic_constants
import numpy as np


class TestMurnaghan(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project = Project(os.path.join(cls.file_location, "test_elast"))
        cls.basis = CrystalStructure(
            element="Fe", bravais_basis="bcc", lattice_constant=2.83
        )
        cls.project.remove_jobs_silently(recursive=True)

    @classmethod
    def tearDownClass(cls):
        cls.file_location = os.path.dirname(os.path.abspath(__file__))
        cls.project.remove(enable=True, enforce=True)

    def test_validate_ready_to_run(self):
        job = self.project.create_job(
            'HessianJob', "job_test"
        )
        job.set_reference_structure(self.basis)
        elast = job.create_job('ElasticTensor', 'elast_job')
        elast.input['min_num_measurements'] = 100
        elast.validate_ready_to_run()
        self.assertEqual(len(elast.input['rotations']), 48)
        self.assertEqual(len(elast.input['strain_matrices']), 100)

    def test_calc_elastic_tensor(self):
        strain = [[[0.005242761019305993, -0.0012053606628952052, -0.0032546722513198236],
                   [-0.0012053606628952052, -0.007127371135828069, 3.041718158272068e-06],
                   [-0.0032546722513198236, 3.041718158272068e-06, 0.006586071210292139]],
                  [[0.0024605481734605596, 0.0016778559303512886, -0.0016403042799169243],
                   [0.0016778559303512886, -8.266496915246613e-05, -0.002503657361580728],
                   [-0.0016403042799169243, -0.002503657361580728, -0.003929702111752978]]]
        stress = [[[1.1937765570713004, -0.2760427230548313, -0.7461461292223824],
                   [-0.2760427230548312, -0.0019061490643724585, -0.00241807297859032],
                   [-0.7461461292223824, -0.0024180729785903197, 1.336480507679297]],
                  [[0.012575219368485248, 0.38833608463518243, -0.3793653414598678],
                   [0.3883360846351825, -0.24167276177709596, -0.5818057844834313],
                   [-0.3793653414598678, -0.5818057844834315, -0.6091472161106105]]]
        C = calc_elastic_tensor(strain=strain, rotations=self.basis.get_symmetry()['rotations'], stress=stress)
        self.assertGreater(np.min(C[:3,:3].diagonal()), 200)
        self.assertLess(np.min(C[:3,:3].diagonal()), 300)
        output = calc_elastic_constants(C)
        self.assertGreater(output['poissons_ratio'], 0)
        self.assertLess(output['poissons_ratio'], 1)
        energy = [-8.02446818, -8.02538938]
        volume = [23.38682061, 23.24219739]
        C = calc_elastic_tensor(strain=strain,
                                rotations=self.basis.get_symmetry()['rotations'],
                                energy=energy,
                                volume=volume)


if __name__ == "__main__":
    unittest.main()
