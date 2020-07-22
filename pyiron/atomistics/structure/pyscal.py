# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from sklearn import cluster
from pyiron.lammps.structure import UnfoldingPrism
from pyiron.atomistics.structure.atoms import pyiron_to_ase

try:
    import pyscal.core as pc
except ImportError:
    pass  # pyscal is currently not available on windows

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Nov 6, 2019"


def get_steinhardt_parameter_job(job, neighbor_method="cutoff", cutoff=0, n_clusters=2, 
                                 q=[4, 6], averaged=False, clustering=True):
    """
    Calculate Steinhardts parameters

    Args:
        job (job): pyiron job
        neighbor_method (str) : can be ['cutoff', 'voronoi']
        cutoff (float) : can be 0 for adaptive cutoff or any other value
        n_clusters (int) : number of clusters for K means clustering
        q (list) : can be from 2-12, the required q values to be calculated
        averaged (bool) : If True, calculates the averaged versions of the parameter
        clustering (bool) : If True, cluster based on the q values

    Returns:
        q (list) : calculated q parameters
    
    """
    return get_steinhardt_parameter_structure(
        structure=job.structure, 
        neighbor_method=neighbor_method,
        cutoff=cutoff, 
        n_clusters=n_clusters, 
        q=q,
        averaged=averaged,
        clustering=cluster
    )


def get_steinhardt_parameter_structure(structure, neighbor_method="cutoff", cutoff=0, n_clusters=2, 
                                       q=[4, 6], averaged=False, clustering=True):
    """
    Calculate Steinhardts parameters

    Args:
        job (job): pyiron job
        neighbor_method (str) : can be ['cutoff', 'voronoi']
        cutoff (float) : can be 0 for adaptive cutoff or any other value
        n_clusters (int) : number of clusters for K means clustering
        q (list) : can be from 2-12, the required q values to be calculated
        averaged (bool) : If True, calculates the averaged versions of the parameter
        clustering (bool) : If True, cluster based on the q values

    Returns:
        q (list) : calculated q parameters
    
    """
    sys = pc.System()
    sys.read_inputfile(
        pyiron_to_ase(structure), 
        format='ase', 
        is_triclinic=not UnfoldingPrism(structure.cell, digits=15).is_skewed()
    )

    sys.find_neighbors(
        method=neighbor_method,
        cutoff=cutoff 
    )

    sys.calculate_q(
        q, 
        averaged=averaged
    )

    sysq = sys.get_qvals(
        q, 
        averaged=averaged
    )

    if clustering:
        cl = cluster.KMeans(
            n_clusters=n_clusters
        )

        ind = cl.fit(list(zip(*sysq))).labels_ == 0
        return sysq, ind
    else:
        return sysq

def analyse_pyscal_centro_symmetry(atoms, num):
    pass