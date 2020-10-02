# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from sklearn.cluster import MeanShift
from scipy.sparse import coo_matrix
from pyiron_base import Settings

__author__ = "Joerg Neugebauer, Sam Waseda"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sam Waseda"
__email__ = "waseda@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


class Neighbors(object):
    """
    Class for storage of the neighbor information for a given atom based on the KDtree algorithm
    """

    def __init__(self, ref_structure, tolerance=2):
        self.distances = None
        self.vecs = None
        self.indices = None
        self._shells = None
        self._tolerance = tolerance
        self._ref_structure = ref_structure
        self._cluster_vecs = None
        self._cluster_dist = None

    @property
    def shells(self):
        """
            Returns the cell numbers of each atom according to the distances
        """
        return self.get_local_shells(tolerance=self._tolerance)

    def update_vectors(self):
        """
            Update vecs and distances with the same indices
        """
        if np.max(np.absolute(self.vecs)) > 0.49*np.min(np.linalg.norm(self._ref_structure.cell, axis=-1)):
            raise AssertionError('Largest distance value is larger than half the box -> rerun get_neighbors')
        myself = np.ones_like(self.indices)*np.arange(len(self.indices))[:, np.newaxis]
        vecs = self._ref_structure.get_distances(
            myself.flatten(), self.indices.flatten(), mic=np.all(self._ref_structure.pbc), vector=True
        )
        self.vecs = vecs.reshape(self.vecs.shape)
        self.distances = np.linalg.norm(self.vecs, axis=-1)

    def get_shell_dict(self, max_shell=2):
        """

        Args:
            max_shell (int): maximum number of shells

        Returns:

        """

        shells = self.shells[0]
        dist = self.distances[0]
        shell_dict = {}
        for i_shell in set(shells):
            if i_shell > max_shell:
                break
            shell_dict[i_shell] = np.mean(dist[shells == i_shell])
            # print ("shells: ", i_shell, shell_dict[i_shell])
        if not (max(shell_dict.keys()) == max_shell):
            raise AssertionError()
        return shell_dict

    def get_local_shells(self, tolerance=2, cluster_by_distances=False, cluster_by_vecs=False):
        """
        Set shell indices based on distances available to each atom. Clustering methods can be used
        at the same time, which will be useful at finite temperature results, but depending on how
        dispersed the atoms are, the algorithm could take some time. If the clustering method(-s)
        have already been launched before this function, it will use the results already available
        and does not execute the clustering method(-s) again.

        Args:
            tolerance (int): decimals in np.round for rounding up distances
            cluster_by_distances (bool): If True, `cluster_by_distances` is called first and the distances obtained
                from the clustered distances are used to calculate the shells. If cluster_by_vecs is True at the same
                time, `cluster_by_distances` will use the clustered vectors for its clustering algorithm. For more,
                see the DocString of `cluster_by_distances`. (default: False)
            cluster_by_vecs (bool): If True, `cluster_by_vectors` is called first and the distances obtained from
                the clustered vectors are used to calculate the shells. (default: False)

        Returns:
            shells (numpy.ndarray): shell indices
        """
        if cluster_by_distances:
            if self._cluster_dist is None:
                self.cluster_by_distances(use_vecs=cluster_by_vecs)
            shells = [np.unique(np.round(dist, decimals=tolerance), return_inverse=True)[1]+1
                         for dist in self._cluster_dist.cluster_centers_[self._cluster_dist.labels_].flatten()
                     ]
            return np.array(shells).reshape(self.indices.shape)
        if cluster_by_vecs:
            if self._cluster_vecs is None:
                self.cluster_by_vecs()
            shells = [np.unique(np.round(dist, decimals=tolerance), return_inverse=True)[1]+1
                         for dist in np.linalg.norm(self._cluster_vecs.cluster_centers_[self._cluster_vecs.labels_], axis=-1)
                     ]
            return np.array(shells).reshape(self.indices.shape)
        if self._shells is None:
            if self.distances is None:
                return None
            self._shells = []
            for dist in self.distances:
                self._shells.append(np.unique(np.round(dist, decimals=tolerance), return_inverse=True)[1]+1)
            if isinstance(self.distances, np.ndarray):
                self._shells = np.array(self._shells)
        return self._shells

    def get_global_shells(self, tolerance=2, cluster_by_distances=False, cluster_by_vecs=False):
        """
        Set shell indices based on all distances available in the system instead of
        setting them according to the local distances (in contrast to shells defined
        as an attribute in this class). Clustering methods can be used at the same time,
        which will be useful at finite temperature results, but depending on how dispersed
        the atoms are, the algorithm could take some time. If the clustering method(-s)
        have already been launched before this function, it will use the results already
        available and does not execute the clustering method(-s) again.

        Args:
            tolerance (int): decimals in np.round for rounding up distances (default: 2)
            cluster_by_distances (bool): If True, `cluster_by_distances` is called first and the distances obtained
                from the clustered distances are used to calculate the shells. If cluster_by_vecs is True at the same
                time, `cluster_by_distances` will use the clustered vectors for its clustering algorithm. For more,
                see the DocString of `cluster_by_distances`. (default: False)
            cluster_by_vecs (bool): If True, `cluster_by_vectors` is called first and the distances obtained from
                the clustered vectors are used to calculate the shells. (default: False)

        Returns:
            shells (numpy.ndarray): shell indices (cf. shells)
        """
        if self.distances is None:
            raise ValueError('neighbors not set')
        distances = self.distances
        if cluster_by_distances:
            if self._cluster_dist is None:
                self.cluster_by_distances(use_vecs=cluster_by_vecs)
            distances = self._cluster_dist.cluster_centers_[self._cluster_dist.labels_].reshape(self.distances.shape)
        elif cluster_by_vecs:
            if self._cluster_vecs is None:
                self.cluster_by_vecs()
            distances = np.linalg.norm(self._cluster_vecs.cluster_centers_[self._cluster_vecs.labels_], axis=-1).reshape(self.distances.shape)
        dist_lst = np.unique(np.round(a=distances, decimals=tolerance))
        shells = distances[:,:,np.newaxis]-dist_lst[np.newaxis,np.newaxis,:]
        shells = np.absolute(shells).argmin(axis=-1)+1
        return shells

    def get_shell_matrix(
        self, chemical_pair=None, cluster_by_distances=False, cluster_by_vecs=False
    ):
        """
        Shell matrices for pairwise interaction. Note: The matrices are always symmetric, meaning if you
        use them as bilinear operators, you have to divide the results by 2.

        Args:
            chemical_pair (list): pair of chemical symbols (e.g. ['Fe', 'Ni'])

        Returns:
            list of sparse matrices for different shells


        Example:
            from pyiron import Project
            structure = Project('.').create_structure('Fe', 'bcc', 2.83).repeat(2)
            J = -0.1 # Ising parameter
            magmoms = 2*np.random.random((len(structure)), 3)-1 # Random magnetic moments between -1 and 1
            neigh = structure.get_neighbors(num_neighbors=8) # Iron first shell
            shell_matrices = neigh.get_shell_matrix()
            print('Energy =', 0.5*J*magmoms.dot(shell_matrices[0].dot(matmoms)))
        """

        pairs = np.stack((self.indices,
            np.ones_like(self.indices)*np.arange(len(self.indices))[:,np.newaxis],
            self.get_global_shells(cluster_by_distances=cluster_by_distances, cluster_by_vecs=cluster_by_vecs)-1),
            axis=-1
        ).reshape(-1, 3)
        shell_max = np.max(pairs[:,-1])+1
        if chemical_pair is not None:
            c = self._ref_structure.get_chemical_symbols()
            pairs = pairs[np.all(np.sort(c[pairs[:,:2]], axis=-1)==np.sort(chemical_pair), axis=-1)]
        shell_matrix = []
        for ind in np.arange(shell_max):
            indices = pairs[ind==pairs[:,-1]]
            if len(indices)>0:
                ind_tmp = np.unique(indices[:,:-1], axis=0, return_counts=True)
                shell_matrix.append(coo_matrix((ind_tmp[1], (ind_tmp[0][:,0], ind_tmp[0][:,1])),
                    shape=(len(self._ref_structure), len(self._ref_structure))
                ))
            else:
                shell_matrix.append(coo_matrix((len(self._ref_structure), len(self._ref_structure))))
        return shell_matrix

    def find_neighbors_by_vector(self, vector, deviation=False):
        """
        Args:
            vector (list/np.ndarray): vector by which positions are translated (and neighbors are searched)
            deviation (bool): whether to return distance between the expect positions and real positions

        Returns:
            np.ndarray: list of id's for the specified translation

        Example:
            a_0 = 2.832
            structure = pr.create_structure('Fe', 'bcc', a_0)
            id_list = structure.find_neighbors_by_vector([0, 0, a_0])
            # In this example, you get a list of neighbor atom id's at z+=a_0 for each atom.
            # This is particularly powerful for SSA when the magnetic structure has to be translated
            # in each direction.
        """

        dist = np.linalg.norm(self.vecs-np.array(vector), axis=-1)
        if deviation:
            return self.indices[np.arange(len(dist)), np.argmin(dist, axis=-1)], np.min(dist, axis=-1)
        return self.indices[np.arange(len(dist)), np.argmin(dist, axis=-1)]

    def cluster_by_vecs(self, bandwidth=None, n_jobs=None, max_iter=300):
        """
        Method to group vectors which have similar values. This method should be used as a part of
        neigh.get_global_shells(cluster_by_vecs=True) or neigh.get_local_shells(cluster_by_vecs=True).
        However, in order to specify certain arguments (such as n_jobs or max_iter), it might help to
        have run this function before calling parent functions, as the data obtained with this function
        will be stored in the variable `_cluster_vecs`

        Args:
            bandwidth (float): Resolution (cf. sklearn.cluster.MeanShift)
            n_jobs (int): Number of cores (cf. sklearn.cluster.MeanShift)
            max_iter (int): Number of maximum iterations (cf. sklearn.cluster.MeanShift)
        """
        if bandwidth is None:
            bandwidth = 0.2*np.min(self.distances)
        dr = self.vecs.copy().reshape(-1, 3)
        self._cluster_vecs = MeanShift(bandwidth=bandwidth, n_jobs=n_jobs, max_iter=max_iter).fit(dr)
        self._cluster_vecs.labels_ = self._cluster_vecs.labels_.reshape(self.indices.shape)

    def cluster_by_distances(self, bandwidth=None, use_vecs=False, n_jobs=None, max_iter=300):
        """
        Method to group vectors which have similar values. This method should be used as a part of
        neigh.get_global_shells(cluster_by_vecs=True) or neigh.get_local_shells(cluster_by_distances=True).
        However, in order to specify certain arguments (such as n_jobs or max_iter), it might help to
        have run this function before calling parent functions, as the data obtained with this function
        will be stored in the variable `_cluster_distances`

        Args:
            bandwidth (float): Resolution (cf. sklearn.cluster.MeanShift)
            use_vecs (bool): Whether to form clusters for vecs beforehand. If true, the distances obtained
                from the clustered vectors is used for the distance clustering.  Otherwise neigh.distances
                is used.
            n_jobs (int): Number of cores (cf. sklearn.cluster.MeanShift)
            max_iter (int): Number of maximum iterations (cf. sklearn.cluster.MeanShift)
        """
        if bandwidth is None:
            bandwidth = 0.05*np.min(self.distances)
        dr = self.distances
        if use_vecs:
            if self._cluster_vecs is None:
                self.cluster_by_vecs()
            dr = np.linalg.norm(self._cluster_vecs.cluster_centers_[self._cluster_vecs.labels_], axis=-1)
        self._cluster_dist = MeanShift(bandwidth=bandwidth, n_jobs=n_jobs, max_iter=max_iter).fit(dr.reshape(-1, 1))
        self._cluster_dist.labels_ = self._cluster_dist.labels_.reshape(self.indices.shape)

    def reset_clusters(self, vecs=True, distances=True):
        """
        Method to reset clusters.

        Args:
            vecs (bool): Reset `_cluster_vecs` (cf. `cluster_by_vecs`)
            distances (bool): Reset `_cluster_distances` (cf. `cluster_by_distances`)
        """
        if vecs:
            self._cluster_vecs = None
        if distances:
            self._cluster_distances = None

