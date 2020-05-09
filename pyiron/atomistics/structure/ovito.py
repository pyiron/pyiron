# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron.base.settings.generic import Settings
from ovito import ObjectNode
from ovito.data import Bonds, DataCollection
from ovito.modifiers import (
    CreateBondsModifier,
    CommonNeighborAnalysisModifier,
    CentroSymmetryModifier,
    VoronoiAnalysisModifier,
)

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"

s = Settings()


def analyse_ovito_cna_adaptive(atoms, mode="total"):
    """
    Use Ovito's common neighbor analysis binding.

    Args:
        atoms (pyrion.structure.atoms.Atoms): The structure to analyze.
        mode ("total"/"numeric"/"str"): Controls the style and level of detail of the output. (Default is "total", only
            return a summary of the values in the structure.)

    Returns:
        (depends on `mode`)
    """
    s.publication_add(publication())
    if not mode in ["total", "numeric", "str"]:
        raise ValueError("Unsupported mode")
    data = DataCollection.create_from_ase_atoms(atoms.copy())
    node = ObjectNode()
    node.source = data
    node.modifiers.append(
        CommonNeighborAnalysisModifier(
            mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff
        )
    )
    output = node.compute()
    if mode == "total":
        return output.attributes
    else:
        atoms_output = output.to_ase_atoms()
        if mode == "numeric":
            return atoms_output.get_array("Structure Type")
        elif mode == "str":
            cna_property = output.particle_properties.structure_type
            return np.array(
                [
                    cna_property.get_type_by_id(cnatype).name
                    for cnatype in atoms_output.get_array("Structure Type")
                ]
            )


def analyse_ovito_centro_symmetry(atoms, num_neighbors=12):
    """

    Args:
        atoms:
        num_neighbors:

    Returns:

    """
    s.publication_add(publication())
    data = DataCollection.create_from_ase_atoms(atoms.copy())
    node = ObjectNode()
    node.source = data
    node.modifiers.append(CentroSymmetryModifier(num_neighbors=num_neighbors))
    output = node.compute()
    return output.particle_properties["Centrosymmetry"].array


def analyse_ovito_voronoi_volume(atoms):
    """

    Args:
        atoms:

    Returns:

    """
    s.publication_add(publication())
    data = DataCollection.create_from_ase_atoms(atoms.copy())
    node = ObjectNode()
    node.source = data
    node.modifiers.append(VoronoiAnalysisModifier())
    output = node.compute()
    return output.particle_properties["Atomic Volume"].array


def publication():
    return {
        "ovito": {
            "Stukowski_2009": {
                "doi": "10.1088/0965-0393/18/1/015012",
                "url": "https://doi.org/10.1088%2F0965-0393%2F18%2F1%2F015012",
                "year": "2009",
                "month": "dec",
                "publisher": "{IOP} Publishing",
                "volume": "18",
                "number": "1",
                "pages": "015012",
                "author": ["Alexander Stukowski"],
                "title": "Visualization and analysis of atomistic simulation data with "
                "{OVITO}{\textendash}the Open Visualization Tool",
                "journal": "Modelling and Simulation in Materials Science and Engineering",
                "abstract": "The Open Visualization Tool (OVITO) is a new 3D visualization "
                "software designed for post-processing atomistic data obtained "
                "from molecular dynamics or Monte Carlo simulations. Unique "
                "analysis, editing and animations functions are integrated into "
                "its easy-to-use graphical user interface. The software is "
                "written in object-oriented C++, controllable via Python scripts "
                "and easily extendable through a plug-in interface. It is "
                "distributed as open-source software and can be downloaded "
                "from the website http://ovito.sourceforge.net/.",
            }
        }
    }
