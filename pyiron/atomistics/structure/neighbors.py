# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.sparse import coo_matrix
from pyiron_base import Settings
from pyiron.atomistics.structure.analyse import get_average_of_unique_labels
import warnings

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

class Tree:
    """
    Class to get tree structure for the neighborhood information.

    Main attributes (do not modify them):

    - distances (numpy.ndarray/list): Distances to the neighbors of given positions
    - indices (numpy.ndarray/list): Indices of the neighbors of given positions
    - vecs (numpy.ndarray/list): Vectors to the neighbors of given positions

    Auxiliary attributes:

    - allow_ragged (bool): Whether to allow a ragged list of numpy arrays or not. This is relevant
        only if the cutoff length was specified, in which case the number of atoms for each
        position could differ.
    - wrap_positions (bool): Whether to wrap back the positions entered by user in get_neighborhood
        etc. Since the information outside the original box is limited to a few layer,
        wrap_positions=False might miss some points without issuing an error.

    Furthermore, you can re-employ the original tree structure to get neighborhood information via
    get_indices, get_vectors, get_distances and get_neighborhood. The information delivered by
    get_neighborhood is the same as what can be obtained through the other three getters (except
    for the fact that get_neighborhood returns a Tree instance, while the others return numpy
    arrays), but since get_vectors anyway has to call get_distances and get_indices internally, the
    computational cost of get_neighborhood and get_vectors is the same.

    """
    def __init__(self, ref_structure):
        self.distances = None
        self.vecs = None
        self.indices = None
        self._extended_positions = None
        self._wrapped_indices = None
        self._allow_ragged = False
        self._cell = None
        self._extended_indices = None
        self._ref_structure = ref_structure.copy()
        self.wrap_positions = False
        self._tree = None
        self.num_neighbors = None
        self.cutoff_radius = np.inf

    def __repr__(self):
        """
            Returns: __repr__
        """
        to_return = (
            "Main attributes:\n"
            + "- distances : Distances to the neighbors of given positions\n"
            + "- indices : Indices of the neighbors of given positions\n"
        )
        if self.vecs is not None:
            to_return += "- vecs : Vectors to the neighbors of given positions\n"
        return to_return

    def copy(self):
        new_neigh = Tree(self._ref_structure)
        new_neigh.distances = self.distances.copy()
        new_neigh.vecs = self.vecs.copy()
        new_neigh.indices = self.indices.copy()
        new_neigh._extended_positions = self._extended_positions
        new_neigh._wrapped_indices = self._wrapped_indices
        new_neigh._allow_ragged = self._allow_ragged
        new_neigh._cell = self._cell
        new_neigh._extended_indices = self._extended_indices
        new_neigh.wrap_positions = self.wrap_positions
        new_neigh._tree = self._tree
        new_neigh.num_neighbors = self.num_neighbors
        new_neigh.cutoff_radius = self.cutoff_radius
        return new_neigh

    def _get_max_length(self, ref_vector=None):
        if ref_vector is None:
            ref_vector = self.distances
        if (ref_vector is None
            or len(ref_vector)==0
            or not hasattr(ref_vector[0], '__len__')):
            return None
        return max(len(dd[dd<np.inf]) for dd in ref_vector)

    def _fill(self, value, filler=np.inf):
        max_length = self._get_max_length()
        if max_length is None:
            return value
        arr = np.zeros((len(value), max_length)+value[0].shape[1:], dtype=type(filler))
        arr.fill(filler)
        for ii, vv in enumerate(value):
            arr[ii,:len(vv)] = vv
        return arr

    def _contract(self, value, ref_vector=None):
        if self._get_max_length(ref_vector=ref_vector) is None:
            return value
        return [vv[:np.sum(dist<np.inf)] for vv, dist in zip(value, self.distances)]

    @property
    def allow_ragged(self):
        """
        Whether to allow ragged list of distancs/vectors/indices or fill empty slots with numpy.inf
        to get rectangular arrays
        """
        return self._allow_ragged

    @allow_ragged.setter
    def allow_ragged(self, new_bool):
        if not isinstance(new_bool, bool):
            raise ValueError('allow_ragged must be a boolean')
        if self._allow_ragged == new_bool:
            return
        self._allow_ragged = new_bool
        if new_bool:
            self.distances = self._contract(self.distances)
            if self.vecs is not None:
                self.vecs = self._contract(self.vecs)
            self.indices = self._contract(self.indices)
        else:
            self.distances = self._fill(self.distances)
            self.indices = self._fill(self.indices, filler=len(self._ref_structure))
            if self.vecs is not None:
                self.vecs = self._fill(self.vecs)

    def _get_extended_positions(self):
        if self._extended_positions is None:
            return self._ref_structure.positions
        return self._extended_positions

    def _get_wrapped_indices(self):
        if self._wrapped_indices is None:
            return np.arange(len(self._ref_structure.positions))
        return self._wrapped_indices

    def _get_wrapped_positions(self, positions):
        if not self.wrap_positions:
            return np.asarray(positions)
        x = np.array(positions).copy()
        cell = self._ref_structure.cell
        x_scale = np.dot(x, np.linalg.inv(cell))+1.0e-12
        x[...,self._ref_structure.pbc] -= np.dot(np.floor(x_scale),
                                                 cell)[...,self._ref_structure.pbc]
        return x

    def _get_distances_and_indices(
        self,
        positions=None,
        allow_ragged=None,
        num_neighbors=None,
        cutoff_radius=np.inf,
        width_buffer=1.2,
    ):
        if positions is None:
            if allow_ragged is None:
                return self.distances, self.indices
            if allow_ragged:
                return (self._contract(self.distances),
                        self._contract(self.indices))
            return (self._fill(self.distances),
                    self._fill(self.indices, filler=len(self._ref_structure)))
        num_neighbors = self._estimate_num_neighbors(
            num_neighbors=num_neighbors,
            cutoff_radius=cutoff_radius,
            width_buffer=width_buffer,
        )
        if len(self._get_extended_positions()) < num_neighbors and cutoff_radius==np.inf:
            raise ValueError(
                'num_neighbors too large - make width_buffer larger and/or make '
                + 'num_neighbors smaller'
            )
        distances, indices = self._tree.query(
            self._get_wrapped_positions(positions),
            k=num_neighbors,
            distance_upper_bound=cutoff_radius
        )
        if cutoff_radius<np.inf and np.any(distances.T[-1]<np.inf):
            warnings.warn(
                'Number of neighbors found within the cutoff_radius is equal to (estimated) '
                + 'num_neighbors. Increase num_neighbors (or set it to None) or '
                + 'width_buffer to find all neighbors within cutoff_radius.'
            )
        self._extended_indices = indices.copy()
        indices[distances<np.inf] = self._get_wrapped_indices()[indices[distances<np.inf]]
        if allow_ragged is None:
            allow_ragged = self.allow_ragged
        if allow_ragged:
            return self._contract(distances, distances), self._contract(indices, distances)
        return distances, indices

    def get_indices(
        self,
        positions=None,
        allow_ragged=None,
        num_neighbors=None,
        cutoff_radius=np.inf,
        width_buffer=1.2,
    ):
        """
        Get current indices or neighbor indices for given positions using the same Tree structure
        used for the instantiation of the Neighbor class. This function should not be used if the
        structure changed in the meantime. If `positions` is None and `allow_ragged` is None, this
        function returns the same content as `indices`.

        Args:
            positions (list/numpy.ndarray/None): Positions around which neighborhood vectors
                are to be computed (None to get current vectors)
            allow_ragged (bool): Whether to allow ragged list of arrays or rectangular
                numpy.ndarray filled with np.inf for values outside cutoff_radius
            num_neighbors (int/None): Number of neighboring atoms to calculate vectors for.
                Ignored if `positions` is None.
            cutoff_radius (float): cutoff radius. Ignored if `positions` is None.
            width_buffer (float): Buffer length for the estimation of num_neighbors. Ignored if
                `positions` is None.

        Returns:
            (list/numpy.ndarray) list (if allow_ragged=True) or numpy.ndarray (otherwise) of
                neighbor indices

        """
        return self._get_distances_and_indices(
            positions=positions,
            allow_ragged=allow_ragged,
            num_neighbors=num_neighbors,
            cutoff_radius=cutoff_radius,
            width_buffer=width_buffer,
        )[1]

    def get_distances(
        self,
        positions=None,
        allow_ragged=None,
        num_neighbors=None,
        cutoff_radius=np.inf,
        width_buffer=1.2,
    ):
        """
        Get current positions or neighbor positions for given positions using the same Tree
        structure used for the instantiation of the Neighbor class. This function should not be
        used if the structure changed in the meantime. If `positions` is None and `allow_ragged` is
        None, this function returns the same content as `positions`.

        Args:
            positions (list/numpy.ndarray/None): Positions around which neighborhood vectors
                are to be computed (None to get current vectors)
            allow_ragged (bool): Whether to allow ragged list of arrays or rectangular
                numpy.ndarray filled with np.inf for values outside cutoff_radius
            num_neighbors (int/None): Number of neighboring atoms to calculate vectors for.
                Ignored if `positions` is None.
            cutoff_radius (float): cutoff radius. Ignored if `positions` is None.
            width_buffer (float): Buffer length for the estimation of num_neighbors. Ignored if
                `positions` is None.

        Returns:
            (list/numpy.ndarray) list (if allow_ragged=True) or numpy.ndarray (otherwise) of
                neighbor distances
        """
        return self._get_distances_and_indices(
            positions=positions,
            allow_ragged=allow_ragged,
            num_neighbors=num_neighbors,
            cutoff_radius=cutoff_radius,
            width_buffer=width_buffer,
        )[0]

    def get_vectors(
        self,
        positions=None,
        allow_ragged=False,
        num_neighbors=None,
        cutoff_radius=np.inf,
        width_buffer=1.2,
    ):
        """
        Get current vectors or neighbor vectors for given positions using the same Tree structure
        used for the instantiation of the Neighbor class. This function should not be used if the
        structure changed in the meantime. If `positions` is None and `allow_ragged` is None, this
        function returns the same content as `vecs`.

        Args:
            positions (list/numpy.ndarray/None): Positions around which neighborhood vectors
                are to be computed (None to get current vectors)
            allow_ragged (bool): Whether to allow ragged list of arrays or rectangular
                numpy.ndarray filled with np.inf for values outside cutoff_radius
            num_neighbors (int/None): Number of neighboring atoms to calculate vectors for.
                Ignored if `positions` is None.
            cutoff_radius (float): cutoff radius. Ignored if `positions` is None.
            width_buffer (float): Buffer length for the estimation of num_neighbors. Ignored if
                `positions` is None.


        Returns:
            (list/numpy.ndarray) list (if allow_ragged=True) or numpy.ndarray (otherwise) of
                neighbor vectors
        """
        return self._get_vectors(
            positions=positions,
            allow_ragged=allow_ragged,
            num_neighbors=num_neighbors,
            cutoff_radius=cutoff_radius,
            width_buffer=width_buffer,
        )

    def _get_vectors(
        self,
        positions=None,
        allow_ragged=None,
        num_neighbors=None,
        cutoff_radius=np.inf,
        distances=None,
        indices=None,
        width_buffer=1.2,
    ):
        if positions is not None:
            if distances is None or indices is None:
                distances, indices = self._get_distances_and_indices(
                    positions=positions,
                    allow_ragged=False,
                    num_neighbors=num_neighbors,
                    cutoff_radius=cutoff_radius,
                    width_buffer=width_buffer,
                )
            vectors = np.zeros(distances.shape+(3,))
            vectors -= self._get_wrapped_positions(positions).reshape(distances.shape[:-1]+(1, 3))
            vectors[distances<np.inf] += self._get_extended_positions()[
                self._extended_indices[distances<np.inf]
            ]
            vectors[distances==np.inf] = np.array(3*[np.inf])
            if self._cell is not None:
                vectors[distances<np.inf] -= self._cell*np.rint(vectors[distances<np.inf]/self._cell)
        elif self.vecs is not None:
            vectors = self.vecs
        else:
            raise AssertionError(
                'vectors not created yet -> put positions or reinitialize with t_vec=True'
            )
        if allow_ragged==self.allow_ragged:
            return vectors
        if allow_ragged:
            return self._contract(vectors, distances)
        return self._fill(vectors)

    def _estimate_num_neighbors(self, num_neighbors=None, cutoff_radius=np.inf, width_buffer=1.2):
        """

        Args:
            num_neighbors (int): number of neighbors
            width_buffer (float): width of the layer to be added to account for pbc.
            cutoff_radius (float): self-explanatory

        Returns:
            Number of atoms required for a given cutoff

        """
        if num_neighbors is None and cutoff_radius==np.inf and self.num_neighbors is None:
            raise ValueError('Specify num_neighbors or cutoff_radius')
        elif num_neighbors is None and self.num_neighbors is None:
            volume = self._ref_structure.get_volume(per_atom=True)
            num_neighbors = max(14, int((1+width_buffer)*4./3.*np.pi*cutoff_radius**3/volume))
        elif num_neighbors is None:
            num_neighbors = self.num_neighbors
        if self.num_neighbors is None:
            self.num_neighbors = num_neighbors
            self.cutoff_radius = cutoff_radius
        if num_neighbors > self.num_neighbors:
            warnings.warn("Taking a larger search area after initialization has the risk of "
                          + "missing neighborhood atoms")
        return num_neighbors

    def _estimate_width(self, num_neighbors=None, cutoff_radius=np.inf, width_buffer=1.2):
        """

        Args:
            num_neighbors (int): number of neighbors
            width_buffer (float): width of the layer to be added to account for pbc.
            cutoff_radius (float): cutoff radius

        Returns:
            Width of layer required for the given number of atoms

        """
        if num_neighbors is None and cutoff_radius==np.inf:
            raise ValueError('Define either num_neighbors or cutoff_radius')
        if all(self._ref_structure.pbc==False):
            return 0
        elif cutoff_radius!=np.inf:
            return cutoff_radius
        pbc = self._ref_structure.pbc
        prefactor = [1, 1/np.pi, 4/(3*np.pi)]
        prefactor = prefactor[sum(pbc)-1]
        width = np.prod(
            (np.linalg.norm(self._ref_structure.cell, axis=-1)-np.ones(3))*pbc+np.ones(3)
        )
        width *= prefactor*np.max([num_neighbors, 8])/len(self._ref_structure)
        cutoff_radius = width_buffer*width**(1/np.sum(pbc))
        return cutoff_radius

    def get_neighborhood(
        self,
        positions,
        num_neighbors=None,
        t_vec=True,
        cutoff_radius=np.inf,
        width_buffer=1.2,
    ):
        """
        Get neighborhood information of `positions`. What it returns is in principle the same as
        `get_neighborhood` in `Atoms`. The only one difference is the reuse of the same Tree
        structure, which makes the algorithm more efficient, but will fail if the reference
        structure changed in the meantime.

        Args:
            position: Position in a box whose neighborhood information is analysed
            num_neighbors (int): Number of nearest neighbors
            t_vec (bool): True: compute distance vectors (pbc are taken into account)
            cutoff_radius (float): Upper bound of the distance to which the search is to be done
            width_buffer (float): Width of the layer to be added to account for pbc.

        Returns:

            pyiron.atomistics.structure.atoms.Tree: Neighbors instances with the neighbor indices,
                distances and vectors

        """
        new_neigh = self.copy()
        return new_neigh._get_neighborhood(
            positions=positions,
            num_neighbors=num_neighbors,
            t_vec=t_vec,
            cutoff_radius=cutoff_radius,
            exclude_self=False,
            width_buffer=width_buffer,
        )

    def _get_neighborhood(
        self,
        positions,
        num_neighbors=12,
        t_vec=True,
        cutoff_radius=np.inf,
        exclude_self=False,
        width_buffer=1.2,
    ):
        start_column = 0
        if exclude_self:
            start_column = 1
            if num_neighbors is not None:
                num_neighbors += 1
        distances, indices = self._get_distances_and_indices(
            positions,
            num_neighbors=num_neighbors,
            cutoff_radius=cutoff_radius,
            width_buffer=width_buffer,
        )
        if num_neighbors is not None:
            self.num_neighbors -= 1
        max_column = np.sum(distances<np.inf, axis=-1).max()
        self.distances = distances[...,start_column:max_column]
        self.indices = indices[...,start_column:max_column]
        self._extended_indices = self._extended_indices[...,start_column:max_column]
        if t_vec:
            self.vecs = self._get_vectors(
                positions=positions, distances=self.distances, indices=self._extended_indices
            )
        return self

    def _check_width(self, width, pbc):
        if any(pbc) and np.prod(self.distances.shape)>0 and self.vecs is not None:
            if np.linalg.norm(
                self._fill(self._contract(self.vecs), filler=0.0)[...,pbc], axis=-1
            ).max() > width:
                return True
        return False


class Neighbors(Tree):
    """
    Class for storage of the neighbor information for a given atom based on the KDtree algorithm

    Main attributes (do not modify them):

    - distances (numpy.ndarray/list): Distances to the neighbors of all atoms
    - indices (numpy.ndarray/list): Indices of the neighbors of all atoms
    - vecs (numpy.ndarray/list): Vectors to the neighbors of all atoms

    Auxiliary attributes:

    - allow_ragged (bool): Whether to allow a ragged list of numpy arrays or not. This is relevant
        only if the cutoff length was specified, in which case the number of neighbor atoms for each
        atom could differ. allow_ragged = False is computationally more efficient in most of the
        methods.
    - wrap_positions (bool): Whether to wrap back the positions entered by user in get_neighborhood
        etc. Since the information outside the original box is limited to a few layer,
        wrap_positions=False might miss some points without issuing an error.

    Furthermore, you can re-employ the original tree structure to get neighborhood information via
    get_indices, get_vectors, get_distances and get_neighborhood. The information delivered by
    get_neighborhood is the same as what can be obtained through the other three getters (except
    for the fact that get_neighborhood returns a Tree instance, while the others return numpy
    arrays), but since get_vectors anyway has to call get_distances and get_indices internally, the
    computational cost of get_neighborhood and get_vectors is the same.
    """

    def __init__(self, ref_structure, tolerance=2):
        super().__init__(ref_structure=ref_structure)
        self._shells = None
        self._tolerance = tolerance
        self._cluster_vecs = None
        self._cluster_dist = None

    def __repr__(self):
        """
            Returns: __repr__
        """
        to_return = super().__repr__()
        return to_return.replace('given positions', 'each atom')

    @property
    def shells(self):
        """
            Returns the cell numbers of each atom according to the distances
        """
        if self._shells is None:
            self._shells = self.get_local_shells()
        if self.allow_ragged:
            return self._contract(self._shells)
        return self._shells

    def get_local_shells(self, tolerance=None, cluster_by_distances=False, cluster_by_vecs=False):
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
        allow_ragged = self.allow_ragged
        if allow_ragged:
            self.allow_ragged = False
        if tolerance is None:
            tolerance = self._tolerance
        if cluster_by_distances:
            if self._cluster_dist is None:
                self.cluster_by_distances(use_vecs=cluster_by_vecs)
            shells = np.array([np.unique(np.round(dist, decimals=tolerance), return_inverse=True)[1]+1
                         for dist in self._cluster_dist.cluster_centers_[self._cluster_dist.labels_]
                     ])
            shells[self._cluster_dist.labels_<0] = -1
            shells = shells.reshape(self.indices.shape)
        elif cluster_by_vecs:
            if self._cluster_vecs is None:
                self.cluster_by_vecs()
            shells = np.array([np.unique(np.round(dist, decimals=tolerance), return_inverse=True)[1]+1
                         for dist in np.linalg.norm(self._cluster_vecs.cluster_centers_[self._cluster_vecs.labels_], axis=-1)
                     ])
            shells[self._cluster_vecs.labels_<0] = -1
            shells = shells.reshape(self.indices.shape)
        else:
            shells = []
            for dist in self.distances:
                shells.append(np.unique(np.round(dist[dist<np.inf], decimals=tolerance), return_inverse=True)[1]+1)
            shells = self._fill(shells, filler=-1)
        self.allow_ragged = allow_ragged 
        if allow_ragged:
            return self._contract(shells)
        return shells

    def get_global_shells(self, tolerance=None, cluster_by_distances=False, cluster_by_vecs=False):
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
        allow_ragged = self.allow_ragged
        if allow_ragged:
            self.allow_ragged = False
        if tolerance is None:
            tolerance = self._tolerance
        if self.distances is None:
            raise ValueError('neighbors not set')
        distances = self.distances
        if cluster_by_distances:
            if self._cluster_dist is None:
                self.cluster_by_distances(use_vecs=cluster_by_vecs)
            distances = self._cluster_dist.cluster_centers_[
                self._cluster_dist.labels_
            ].reshape(self.distances.shape)
            distances[self._cluster_dist.labels_<0] = np.inf
        elif cluster_by_vecs:
            if self._cluster_vecs is None:
                self.cluster_by_vecs()
            distances = np.linalg.norm(
                self._cluster_vecs.cluster_centers_[self._cluster_vecs.labels_], axis=-1
            ).reshape(self.distances.shape)
            distances[self._cluster_vecs.labels_<0] = np.inf
        dist_lst = np.unique(np.round(a=distances, decimals=tolerance))
        shells = -np.ones_like(self.indices).astype(int)
        shells[distances<np.inf] = np.absolute(
            distances[distances<np.inf, np.newaxis]-dist_lst[np.newaxis, dist_lst<np.inf]
        ).argmin(axis=-1)+1
        self.allow_ragged = allow_ragged 
        if allow_ragged:
            return self._contract(shells)
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

    def find_neighbors_by_vector(self, vector, return_deviation=False):
        """
        Args:
            vector (list/np.ndarray): vector by which positions are translated (and neighbors are searched)
            return_deviation (bool): whether to return distance between the expect positions and real positions

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

        z = np.zeros(len(self._ref_structure)*3).reshape(-1, 3)
        v = np.append(z[:,np.newaxis,:], self.vecs, axis=1)
        dist = np.linalg.norm(v-np.array(vector), axis=-1)
        indices = np.append(np.arange(len(self._ref_structure))[:,np.newaxis], self.indices, axis=1)
        if return_deviation:
            return indices[np.arange(len(dist)), np.argmin(dist, axis=-1)], np.min(dist, axis=-1)
        return indices[np.arange(len(dist)), np.argmin(dist, axis=-1)]

    def cluster_by_vecs(
        self, distance_threshold=None, n_clusters=None, linkage='complete', affinity='euclidean'
    ):
        """
        Method to group vectors which have similar values. This method should be used as a part of
        neigh.get_global_shells(cluster_by_vecs=True) or neigh.get_local_shells(cluster_by_vecs=True).
        However, in order to specify certain arguments (such as n_jobs or max_iter), it might help to
        have run this function before calling parent functions, as the data obtained with this function
        will be stored in the variable `_cluster_vecs`

        Args:
            distance_threshold (float/None): The linkage distance threshold above which, clusters
                will not be merged. (cf. sklearn.cluster.AgglomerativeClustering)
            n_clusters (int/None): The number of clusters to find.
                (cf. sklearn.cluster.AgglomerativeClustering)
            linkage (str): Which linkage criterion to use. The linkage criterion determines which
                distance to use between sets of observation. The algorithm will merge the pairs of
                cluster that minimize this criterion. (cf. sklearn.cluster.AgglomerativeClustering)
            affinity (str/callable): Metric used to compute the linkage. Can be `euclidean`, `l1`,
                `l2`, `manhattan`, `cosine`, or `precomputed`. If linkage is `ward`, only
                `euclidean` is accepted.

        """
        allow_ragged = self.allow_ragged
        if allow_ragged:
            self.allow_ragged = False
        if distance_threshold is None and n_clusters is None:
            distance_threshold = np.min(self.distances)
        dr = self.vecs[self.distances<np.inf]
        self._cluster_vecs = AgglomerativeClustering(
            distance_threshold=distance_threshold,
            n_clusters=n_clusters,
            linkage=linkage,
            affinity=affinity,
        ).fit(dr)
        self._cluster_vecs.cluster_centers_ = get_average_of_unique_labels(
            self._cluster_vecs.labels_, dr
        )
        new_labels = -np.ones_like(self.indices).astype(int)
        new_labels[self.distances<np.inf] = self._cluster_vecs.labels_
        self._cluster_vecs.labels_ = new_labels
        self.allow_ragged = allow_ragged

    def cluster_by_distances(
        self,
        distance_threshold=None,
        n_clusters=None,
        linkage='complete',
        affinity='euclidean',
        use_vecs=False,
    ):
        """
        Method to group vectors which have similar lengths. This method should be used as a part of
        neigh.get_global_shells(cluster_by_vecs=True) or
        neigh.get_local_shells(cluster_by_distances=True).  However, in order to specify certain
        arguments (such as n_jobs or max_iter), it might help to have run this function before
        calling parent functions, as the data obtained with this function will be stored in the
        variable `_cluster_distances`

        Args:
            distance_threshold (float/None): The linkage distance threshold above which, clusters
                will not be merged. (cf. sklearn.cluster.AgglomerativeClustering)
            n_clusters (int/None): The number of clusters to find.
                (cf. sklearn.cluster.AgglomerativeClustering)
            linkage (str): Which linkage criterion to use. The linkage criterion determines which
                distance to use between sets of observation. The algorithm will merge the pairs of
                cluster that minimize this criterion. (cf. sklearn.cluster.AgglomerativeClustering)
            affinity (str/callable): Metric used to compute the linkage. Can be `euclidean`, `l1`,
                `l2`, `manhattan`, `cosine`, or `precomputed`. If linkage is `ward`, only
                `euclidean` is accepted.
            use_vecs (bool): Whether to form clusters for vecs beforehand. If true, the distances
                obtained from the clustered vectors is used for the distance clustering. Otherwise
                neigh.distances is used.
        """
        allow_ragged = self.allow_ragged
        if allow_ragged:
            self.allow_ragged = False
        if distance_threshold is None:
            distance_threshold = 0.1*np.min(self.distances)
        dr = self.distances[self.distances<np.inf]
        if use_vecs:
            if self._cluster_vecs is None:
                self.cluster_by_vecs()
            labels_to_consider = self._cluster_vecs.labels_[self._cluster_vecs.labels_>=0]
            dr = np.linalg.norm(self._cluster_vecs.cluster_centers_[labels_to_consider], axis=-1)
        self._cluster_dist = AgglomerativeClustering(
            distance_threshold=distance_threshold,
            n_clusters=n_clusters,
            linkage=linkage,
            affinity=affinity,
        ).fit(dr.reshape(-1, 1))
        self._cluster_dist.cluster_centers_ = get_average_of_unique_labels(
            self._cluster_dist.labels_, dr
        )
        new_labels = -np.ones_like(self.indices).astype(int)
        new_labels[self.distances<np.inf] = self._cluster_dist.labels_
        self._cluster_dist.labels_ = new_labels
        self.allow_ragged = allow_ragged

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

    def cluster_analysis(
        self, id_list, return_cluster_sizes=False
    ):
        """

        Args:
            id_list:
            return_cluster_sizes:

        Returns:

        """
        self._cluster = [0] * len(self._ref_structure)
        c_count = 1
        # element_list = self.get_atomic_numbers()
        for ia in id_list:
            # el0 = element_list[ia]
            nbrs = self.indices[ia]
            # print ("nbrs: ", ia, nbrs)
            if self._cluster[ia] == 0:
                self._cluster[ia] = c_count
                self.__probe_cluster(c_count, nbrs, id_list)
                c_count += 1

        cluster = np.array(self._cluster)
        cluster_dict = {
            i_c: np.where(cluster == i_c)[0].tolist() for i_c in range(1, c_count)
        }
        if return_cluster_sizes:
            sizes = [self._cluster.count(i_c + 1) for i_c in range(c_count - 1)]
            return cluster_dict, sizes

        return cluster_dict  # sizes

    def __probe_cluster(self, c_count, neighbors, id_list):
        """

        Args:
            c_count:
            neighbors:
            id_list:

        Returns:

        """
        for nbr_id in neighbors:
            if self._cluster[nbr_id] == 0:
                if nbr_id in id_list:  # TODO: check also for ordered structures
                    self._cluster[nbr_id] = c_count
                    nbrs = self.indices[nbr_id]
                    self.__probe_cluster(c_count, nbrs, id_list)

    # TODO: combine with corresponding routine in plot3d
    def get_bonds(self, radius=np.inf, max_shells=None, prec=0.1):
        """

        Args:
            radius:
            max_shells:
            prec: minimum distance between any two clusters (if smaller considered to be single cluster)

        Returns:

        """

        def get_cluster(dist_vec, ind_vec, prec=prec):
            ind_where = np.where(np.diff(dist_vec) > prec)[0] + 1
            ind_vec_cl = [np.sort(group) for group in np.split(ind_vec, ind_where)]
            return ind_vec_cl

        dist = self.distances
        ind = self.indices
        el_list = self._ref_structure.get_chemical_symbols()

        ind_shell = []
        for d, i in zip(dist, ind):
            id_list = get_cluster(d[d < radius], i[d < radius])
            # print ("id: ", d[d<radius], id_list, dist_lst)
            ia_shells_dict = {}
            for i_shell_list in id_list:
                ia_shell_dict = {}
                for i_s in i_shell_list:
                    el = el_list[i_s]
                    if el not in ia_shell_dict:
                        ia_shell_dict[el] = []
                    ia_shell_dict[el].append(i_s)
                for el, ia_lst in ia_shell_dict.items():
                    if el not in ia_shells_dict:
                        ia_shells_dict[el] = []
                    if max_shells is not None:
                        if len(ia_shells_dict[el]) + 1 > max_shells:
                            continue
                    ia_shells_dict[el].append(ia_lst)
            ind_shell.append(ia_shells_dict)
        return ind_shell

