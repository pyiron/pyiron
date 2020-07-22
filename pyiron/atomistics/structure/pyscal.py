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


def get_steinhardt_parameter_job(job, cutoff=0, n_clusters=2, q=[4, 6]):
    return get_steinhardt_parameter_structure(
        structure=job.get_stucture(), 
        cutoff=cutoff, 
        n_clusters=n_clusters, 
        q=q
    )


def get_steinhardt_parameter_structure(structure, cutoff=0, n_clusters=2, q=[4, 6]):
    sys = pc.System()
    sys.read_inputfile(
        pyiron_to_ase(structure), 
        format='ase', 
        is_triclinic=not UnfoldingPrism(structure.cell, digits=15).is_skewed()
    )
    sys.find_neighbors(
        method='cutoff', 
        cutoff=cutoff
    )
    sys.calculate_q(
        q, 
        averaged=True
    )
    sysq = sys.get_qvals(
        q, 
        averaged=True
    )
    cl = cluster.KMeans(
        n_clusters=n_clusters
    )
    ind = cl.fit(list(zip(*sysq))).labels_ == 0
    return sysq, ind
