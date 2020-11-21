# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron_base import Settings
from sklearn.cluster import AgglomerativeClustering
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

class Analyse:
    """
    Class to analyse atom structure.
    """
    def __init__(self, ref_structure):
        self._ref_structure = ref_structure

    def __repr__(self):
        """
            Returns: __repr__
        """
        return ''

    def get_layers(self, distance_threshold=0.01, id_list=None):
        """
        Get an array of layer numbers.

        Args:
            distance_threshold (float): Distance below which two points are
                considered to belong to the same layer. For detailed
                description: sklearn.cluster.AgglomerativeClustering
            id_list (list/numpy.ndarray): List of atoms for which the layers
                should be considered.

        Returns: Array of layer numbers (same shape as structure.positions)

        Example I - how to get the number of layers in each direction:

        >>> structure = Project('.').create_structure('Fe', 'bcc', 2.83).repeat(5)
        >>> print('Numbers of layers:', np.max(structure.analyse.get_layers(), axis=0)+1)

        Example II - get layers of only one species:

        >>> print('Iron layers:', structure.analyse.get_layers(
        ...     id_list=structure.select_index('Fe'))
        ... )
        """
        if id_list is None:
            id_list = np.arange(len(self._ref_structure))
        layers = []
        for x in self._ref_structure.positions[np.array(id_list)].T:
            cluster = AgglomerativeClustering(
                linkage='complete',
                n_clusters=None,
                distance_threshold=distance_threshold
            ).fit(x.reshape(-1,1))
            args = np.unique(cluster.labels_, return_index=True)[1]
            layers.append(x[args].argsort().argsort()[cluster.labels_])
        return np.vstack(layers).T
