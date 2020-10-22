# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function

import numpy as np
import scipy.constants
from pyiron.atomistics.master.parallel import AtomisticParallelMaster
from pyiron_base import JobGenerator
from sklearn.linear_model import LinearRegression
from collections import defaultdict

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


eV_div_A3_to_GPa = (
    1e21 / scipy.constants.physical_constants["joule-electron volt relationship"][0]
)

def calc_elastic_tensor(strain, stress=None, energy=None, volume=None, rotations=None):
    """
    Calculate 6x6-elastic tensor from the strain and stress or strain and energy+volume.

    Rotations matrices can be added to take box symmetries into account (unit matrix can
    be added but does not have to be included in the list)

    Args:
        strain (numpy.ndarray): nx3x3 strain tensors
        stress (numpy.ndarray): nx3x3 stress tensors
        energy (numpy.ndarray): n energy values
        volume (numpy.ndarray): n volume values
        rotations (numpy.ndarray): mx3x3 rotation matrices
    """
    if len(strain)==0:
        raise ValueError('Not enough points')
    if rotations is not None:
        rotations = np.append(np.eye(3), rotations).reshape(-1, 3, 3)
        _, indices = np.unique(np.round(rotations, 6), axis=0, return_inverse=True)
        rotations = rotations[np.unique(indices)]
    else:
        rotations = np.eye(3)
    strain = np.einsum('nik,mkl,njl->nmij', rotations, strain, rotations).reshape(-1, 3, 3)
    strain = np.stack((strain[:,0,0], strain[:,1,1], strain[:,2,2], 2*strain[:,1,2], 2*strain[:,0,2], 2*strain[:,0,1]), axis=-1)
    coeff = []
    if stress is not None and len(stress)*len(rotations)==len(strain):
        stress = np.einsum('nik,mkl,njl->nmij', rotations, stress, rotations).reshape(-1, 3, 3)
        stress = np.stack((stress[:,0,0], stress[:,1,1], stress[:,2,2], stress[:,1,2], stress[:,0,2], stress[:,0,1]), axis=-1)
        for ii in range(6):
            reg = LinearRegression().fit(strain, stress[:,ii])
            coeff.append(reg.coef_)
        coeff = np.array(coeff)
    elif energy is not None and volume is not None and len(energy)*len(rotations)==len(strain) and len(energy)==len(volume):
        energy = np.tile(energy, len(rotations))
        volume = np.tile(volume, len(rotations))
        C = np.einsum('n,ni,nj->nij', 0.5*volume, strain, strain)
        C = np.triu(C).reshape(-1, 36)
        C = C[:,np.sum(C, axis=0)!=0]
        C = np.append(np.ones(len(C)).reshape(-1, 1), C, axis=-1)
        reg = LinearRegression().fit(C, energy)
        coeff = np.triu(np.ones((6,6))).flatten()
        coeff[coeff!=0] *= reg.coef_[1:]*eV_div_A3_to_GPa
        coeff = coeff.reshape(6,6)
        coeff = 0.5*(coeff+coeff.T)
    else:
        raise ValueError('Problem with fitting data')
    return coeff

def calc_elastic_constants(elastic_tensor):
    """
    Calculate elastic constants from the elastic tensor.

    For anistropic material (i.e. zener_ratio!=1), the values may or may not make
    sense -> don't trust the results straightforwardly

    Args:
        elastic_tensor (numpy.ndarray): 6x6 tensor
    """
    output = {}
    output['elastic_tensor'] = elastic_tensor
    output['lame_coefficient'] = np.mean(elastic_tensor[:3, :3].diagonal())
    output['shear_modulus'] = np.mean(elastic_tensor[3:, 3:].diagonal())
    output['bulk_modulus'] = np.mean(elastic_tensor[:3,:3])
    output['poissons_ratio'] = 0.5/(1+output['shear_modulus']/output['lame_coefficient'])
    output['youngs_modulus'] = 1/np.mean(np.linalg.inv(elastic_tensor[:3,:3]).diagonal())
    output['poissons_ratio'] = -output['youngs_modulus']*np.sum(np.linalg.inv(elastic_tensor[:3,:3]))/6+0.5
    output['zener_ratio'] = 12*np.mean(elastic_tensor[3:,3:].diagonal())/(3*np.trace(elastic_tensor[:3,:3])-np.sum(elastic_tensor[:3,:3]))
    return output

class ElasticJobGenerator(JobGenerator):
    @property
    def parameter_list(self):
        """

        Returns:
            (list)
        """
        parameter_lst = []
        for ii, epsilon in enumerate(self._job.input['strain_matrices']):
            basis = self._job.ref_job.structure.copy()
            basis.apply_strain(np.array(epsilon))
            parameter_lst.append([ii, basis])
        return parameter_lst

    def job_name(self, parameter):
        return "{}_{}".format(self._job.job_name, parameter[0]).replace(".", "_")

    @staticmethod
    def modify_job(job, parameter):
        job.structure = parameter[1]
        return job

# ToDo: not all abstract methods implemented
class ElasticTensor(AtomisticParallelMaster):
    """
    Class to calculate the elastic tensor and isotropic elastic constants

    Example:

    >>> job = pr.create_job('SomeAtomisticJob', 'atomistic')
    >>> job.structure = pr.create_structure('Fe', 'bcc', 2.83)
    >>> elastic = job.create_job('ElasticTensor', 'elastic')
    >>> elastic.run()

    This class is still under construction and there's no guarantee that it does
    what it looks like it does.
    """
    def __init__(self, project, job_name="elastic"):
        """

        Args:
            project:
            job_name:
        """
        super().__init__(project, job_name)
        self.__name__ = "ElasticTensor"
        self.__version__ = "0.1.0"

        # print ("h5_path: ", self.project_hdf5._h5_path)

        # define default input
        self.input["min_num_measurements"] = (11, "minimum number of samples to be taken")
        self.input["min_num_points"] = (105, "minimum number of data points"
            + "(number of measurements will be min_num_points/len(rotations))")
        self.input["max_strain"] = (
            0.01,
            "relative volume variation around volume defined by ref_ham",
        )
        self.input['strain_matrices'] = []
        self.input['use_symmetry'] = True
        self.input['rotations'] = []
        self.input['use_elements'] = (True, 'whether or not consider chemical elements for '
                                            + 'the symmetry analysis. Could be useful for SQS')
        self._job_generator = ElasticJobGenerator(self)

    @property
    def _number_of_measurements(self):
        return max(
            self.input['min_num_measurements'],
            int(self.input['min_num_points']/max(len(self.input['rotations']), 1))
        )

    def validate_ready_to_run(self):
        super().validate_ready_to_run()
        if self.input['use_symmetry'] and len(self.input['rotations'])==0:
            rotations = self.ref_job.structure.get_symmetry(use_elements=self.input['use_elements'])['rotations']
            _, indices = np.unique(np.round(rotations, 6), axis=0, return_inverse=True)
            rotations = rotations[np.unique(indices)]
            self.input['rotations'] = rotations.tolist()
        if len(self.input['strain_matrices'])==0:
            eps_lst = np.random.random((int(self._number_of_measurements), 3, 3))-0.5
            eps_lst *= self.input['max_strain']
            eps_lst += np.einsum('nij->nji', eps_lst)
            self.input['strain_matrices'] = eps_lst.tolist()

    def list_structures(self):
        if self.ref_job.structure is not None:
            return [parameter[1] for parameter in self._job_generator.parameter_list]
        else:
            return []

    def collect_output(self):
        if self.ref_job.server.run_mode.interactive:
            ham = self.project_hdf5.inspect(self.child_ids[0])
            for key in ['energy_tot', 'pressures', 'volume']:
                if key in ham["output/generic"].list_nodes():
                    self._output[key.split('_')[0]] = ham["output/generic/{}".format(key)]
                else:
                    self._output[key.split('_')[0]] = []
        else:
            output_dict = defaultdict(list)
            for job_id in self.child_ids:
                ham = self.project_hdf5.inspect(job_id)
                print("job_id: ", job_id, ham.status)
                for key in ['energy_tot', 'energy_pot', 'pressures', 'volume']:
                    if key in ham["output/generic"].list_nodes():
                        output_dict[key.split('_')[0]].append(ham["output/generic/{}".format(key)][-1])
                    else:
                        output_dict[key.split('_')[0]] = []
                output_dict['id'].append(job_id)
            for k,v in output_dict.items():
                self._output[k] = np.array(v)
        elastic_tensor = calc_elastic_tensor(
            strain = self.input['strain_matrices'],
            stress = -np.array(self._output['pressures']),
            rotations = self.input['rotations'],
            energy = self._output['energy'],
            volume = self._output['volume'],
        )
        for k,v in calc_elastic_constants(elastic_tensor).items():
            self._output[k] = v
        with self.project_hdf5.open("output") as hdf5_out:
            for key, val in self._output.items():
                hdf5_out[key] = val


