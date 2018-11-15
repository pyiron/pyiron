# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from ovito import ObjectNode
from ovito.data import Bonds, DataCollection
from ovito.modifiers import CreateBondsModifier, CommonNeighborAnalysisModifier, CentroSymmetryModifier, VoronoiAnalysisModifier

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


def analyse_ovito_cna_adaptive(atoms, mode='total'):
    """
    Args:
        mode (str): ['total', 'numeric', 'str']
    """
    if not mode in ['total', 'numeric', 'str']:
        raise ValueError('Unsupported mode')
    data = DataCollection.create_from_ase_atoms(atoms.copy())
    node = ObjectNode()
    node.source = data
    node.modifiers.append(CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff))
    output = node.compute()
    if mode == 'total':
        return output.attributes
    else:
        atoms_output = output.to_ase_atoms()
        if mode == 'numeric':
            return atoms_output.get_array("Structure Type")
        elif mode == 'str':
            cna_property = output.particle_properties.structure_type
            return np.array([cna_property.get_type_by_id(cnatype).name for cnatype in atoms_output.get_array("Structure Type")])

def analyse_ovito_centro_symmetry(atoms, num_neighbors=12):
    """
    Args:
        mode (str): ['total', 'numeric', 'str']
    """
    data = DataCollection.create_from_ase_atoms(atoms.copy())
    node = ObjectNode()
    node.source = data
    node.modifiers.append(CentroSymmetryModifier(num_neighbors = num_neighbors))
    output = node.compute()
    return output.particle_properties['Centrosymmetry'].array

def analyse_ovito_voronoi_volume(atoms):
    """
    Args:
        mode (str): ['total', 'numeric', 'str']
    """
    data = DataCollection.create_from_ase_atoms(atoms.copy())
    node = ObjectNode()
    node.source = data
    node.modifiers.append(VoronoiAnalysisModifier())
    output = node.compute()
    return output.particle_properties['Atomic Volume'].array
