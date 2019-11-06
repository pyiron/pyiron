# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import pyscal.core as pc
from sklearn import cluster
from pyiron.lammps.structure import UnfoldingPrism

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Nov 6, 2019"


def _get_steinhardt_parameter(cell, positions, cutoff=3.50, n_clusters=2):
    sys = pc.System()
    prism = UnfoldingPrism(cell, digits=15)
    xhi, yhi, zhi, xy, xz, yz = prism.get_lammps_prism_str()
    coords = [prism.pos_to_lammps(position) for position in positions]
    sys.box = [[0.0, float(xhi)], [0.0, float(yhi)], [0.0, float(zhi)]]
    sys.atoms = [pc.Atom(pos=p, id=i) for i, p in enumerate(coords)]
    sys.find_neighbors(method='cutoff', cutoff=cutoff)
    sys.calculate_q([4, 6], averaged=True)
    sysq = sys.get_qvals([4, 6], averaged=True)
    cl = cluster.KMeans(n_clusters=n_clusters)
    ind = cl.fit(list(zip(*sysq))).labels_ == 0
    return sysq, ind


def get_steinhardt_parameter_job(job, cutoff=3.50, n_clusters=2):
    cell = job['output/generic/cells'][-1]
    positions = job['output/generic/positions'][-1]
    return _get_steinhardt_parameter(cell=cell, positions=positions, cutoff=cutoff, n_clusters=n_clusters)


def get_steinhardt_parameter_structure(structure, cutoff=3.50, n_clusters=2):
    cell = structure.cell.copy()
    positions = structure.positions.copy()
    return _get_steinhardt_parameter(cell=cell, positions=positions, cutoff=cutoff, n_clusters=n_clusters)
