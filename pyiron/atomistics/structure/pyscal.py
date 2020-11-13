# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
from pyiron_base import Settings
from pyiron.atomistics.structure.atoms import pyiron_to_ase
import pyscal.core as pc
from sklearn import cluster

__author__ = "Sarath Menon, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sarath Menon"
__email__ = "sarath.menon@rub.de"
__status__ = "development"
__date__ = "Nov 6, 2019"

s = Settings()


def get_steinhardt_parameter_structure(atoms, neighbor_method="cutoff", cutoff=0, n_clusters=2,
                                       q=(4, 6), averaged=False, clustering=True):
    """
    Calculate Steinhardts parameters

    Args:
        atoms: Atoms object
        neighbor_method (str) : can be ['cutoff', 'voronoi']
        cutoff (float) : can be 0 for adaptive cutoff or any other value
        n_clusters (int) : number of clusters for K means clustering
        q (list) : can be from 2-12, the required q values to be calculated
        averaged (bool) : If True, calculates the averaged versions of the parameter
        clustering (bool) : If True, cluster based on the q values

    Returns:
        q (list) : calculated q parameters

    """
    s.publication_add(publication())
    sys = pc.System()
    sys.read_inputfile(
        pyiron_to_ase(atoms),
        format='ase'
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


def analyse_centro_symmetry(atoms, num_neighbors=12):
    """
    Analyse centrosymmetry parameter

    Args:
        atoms: Atoms object
        num_neighbors (int) : number of neighbors

    Returns:
        csm (list) : list of centrosymmetry parameter
    """
    s.publication_add(publication())
    sys = pc.System()
    sys.read_inputfile(atoms, format="ase")
    sys.calculate_centrosymmetry(nmax=num_neighbors)
    atoms = sys.atoms
    return np.array([atom.centrosymmetry for atom in atoms])


def analyse_diamond_structure(atoms, mode="total", ovito_compatibility=False):
    """
    Analyse diamond structure

    Args:
        atoms: Atoms object
        mode ("total"/"numeric"/"str"): Controls the style and level
        of detail of the output.
            - total : return number of atoms belonging to each structure
            - numeric : return a per atom list of numbers- 0 for unknown,
                1 fcc, 2 hcp, 3 bcc and 4 icosa
            - str : return a per atom string of sructures
        ovito_compatibility(bool): use ovito compatiblity mode

    Returns:
        (depends on `mode`)
    """
    s.publication_add(publication())
    sys = pc.System()
    sys.read_inputfile(atoms, format="ase")
    diamond_dict = sys.identify_diamond()

    ovito_identifiers = [
        'IdentifyDiamond.counts.CUBIC_DIAMOND',
        'IdentifyDiamond.counts.CUBIC_DIAMOND_FIRST_NEIGHBOR',
        'IdentifyDiamond.counts.CUBIC_DIAMOND_SECOND_NEIGHBOR',
        'IdentifyDiamond.counts.HEX_DIAMOND',
        'IdentifyDiamond.counts.HEX_DIAMOND_FIRST_NEIGHBOR',
        'IdentifyDiamond.counts.HEX_DIAMOND_SECOND_NEIGHBOR',
        'IdentifyDiamond.counts.OTHER'
    ]
    pyscal_identifiers = [
        'others', 'fcc', 'hcp', 'bcc', 'ico', 'cubic diamond',
        'cubic diamond 1NN', 'cubic diamond 2NN',
        'hex diamond', 'hex diamond 1NN', 'hex diamond 2NN'
    ]
    convert_to_ovito = {
        0: 6, 1: 6, 2: 6, 3: 6, 4: 6,
        5: 0, 6: 1, 7: 2, 8: 3, 9: 4, 10: 5
    }

    if mode == "total":
        if not ovito_compatibility:
            return diamond_dict
        else:
            return {
                'IdentifyDiamond.counts.CUBIC_DIAMOND': diamond_dict['cubic diamond'],
                'IdentifyDiamond.counts.CUBIC_DIAMOND_FIRST_NEIGHBOR': diamond_dict['cubic diamond 1NN'],
                'IdentifyDiamond.counts.CUBIC_DIAMOND_SECOND_NEIGHBOR': diamond_dict['cubic diamond 2NN'],
                'IdentifyDiamond.counts.HEX_DIAMOND': diamond_dict['hex diamond'],
                'IdentifyDiamond.counts.HEX_DIAMOND_FIRST_NEIGHBOR': diamond_dict['hex diamond 1NN'],
                'IdentifyDiamond.counts.HEX_DIAMOND_SECOND_NEIGHBOR': diamond_dict['hex diamond 2NN'],
                'IdentifyDiamond.counts.OTHER':
                    diamond_dict['others'] +
                    diamond_dict['fcc'] +
                    diamond_dict['hcp'] +
                    diamond_dict['bcc'] +
                    diamond_dict['ico']
            }
    elif mode == "numeric":
        if not ovito_compatibility:
            return [atom.structure for atom in sys.atoms]
        else:
            return [convert_to_ovito[atom.structure] for atom in sys.atoms]
    elif mode == "str":
        if not ovito_compatibility:
            return [pyscal_identifiers[atom.structure] for atom in sys.atoms]
        else:
            return [ovito_identifiers[convert_to_ovito[atom.structure]] for atom in sys.atoms]
    else:
        raise ValueError("Only total, str and numeric mode is imported for analyse_diamond_structure()")


def analyse_cna_adaptive(atoms, mode="total", ovito_compatibility=False):
    """
    Use common neighbor analysis

    Args:
        atoms (pyiron.structure.atoms.Atoms): The structure to analyze.
        mode ("total"/"numeric"/"str"): Controls the style and level
            of detail of the output.
            - total : return number of atoms belonging to each structure
            - numeric : return a per atom list of numbers- 0 for unknown,
                1 fcc, 2 hcp, 3 bcc and 4 icosa
            - str : return a per atom string of sructures
        ovito_compatibility(bool): use ovito compatiblity mode

    Returns:
        (depends on `mode`)
    """
    s.publication_add(publication())
    if mode not in ["total", "numeric", "str"]:
        raise ValueError("Unsupported mode")

    pyscal_parameter = ['others', 'fcc', 'hcp', 'bcc', 'ico']
    ovito_parameter = [
        'CommonNeighborAnalysis.counts.OTHER',
        'CommonNeighborAnalysis.counts.FCC',
        'CommonNeighborAnalysis.counts.HCP',
        'CommonNeighborAnalysis.counts.BCC',
        'CommonNeighborAnalysis.counts.ICO'
    ]

    sys = pc.System()
    sys.read_inputfile(atoms, format="ase")
    cna = sys.calculate_cna()

    if mode == "total":
        if not ovito_compatibility:
            return cna
        else:
            return {o: cna[p] for o, p in zip(
                ovito_parameter,
                pyscal_parameter
            )}
    else:
        atoms = sys.atoms
        cnalist = ([atom.structure for atom in atoms])
        if mode == "numeric":
            return cnalist
        elif mode == "str":
            if not ovito_compatibility:
                dd = ['others', 'fcc', 'hcp', 'bcc', 'ico']
                return [dd[int(x)] for x in cnalist]
            else:
                dd = ["OTHER", "FCC", "HCP", "BCC", "ICO"]
                return [dd[int(x)] for x in cnalist]
        else:
            raise ValueError("Only total, str and numeric mode is imported for analyse_cna_adaptive()")


def analyse_voronoi_volume(atoms):
    """
    Calculate the Voronoi volume of atoms

    Args:
        atoms : (pyiron.structure.atoms.Atoms): The structure to analyze.
    """
    s.publication_add(publication())
    sys = pc.System()
    sys.read_inputfile(atoms, format="ase")
    sys.find_neighbors(method="voronoi")
    atoms = sys.atoms
    vols = np.array([atom.volume for atom in atoms])
    return vols


def publication():
    return {
        "pyscal": {
            "Menon2019": {
                "doi": "10.21105/joss.01824",
                "url": "https://doi.org/10.21105/joss.01824",
                "year": "2019",
                "volume": "4",
                "number": "43",
                "pages": "1824",
                "author": ["Sarath Menon", "Grisell Diaz Leines", "Jutta Rogal"],
                "title": "pyscal: A python module for structural analysis of atomic environments",
                "journal": "Journal of Open Source Software",
            }
        }
    }
