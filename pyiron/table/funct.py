# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import json
import numpy as np
import warnings

from pyiron.atomistics.structure.atoms import Atoms, pyiron_to_ase


def _get_value_from_incar(job, key):
    return eval(
        job["input/incar/data_dict"]["Value"][
            job["input/incar/data_dict"]["Parameter"].index(key)
        ]
    )


def get_majority(lst, minority=False):
    elements_dict = {name: lst.count(name) for name in set(lst)}
    max_value = np.max(list(elements_dict.values()))
    majority_element = [
        key for key, value in elements_dict.items() if value == max_value
    ][0]
    if minority:
        minority_lst = list(elements_dict.keys())
        del minority_lst[minority_lst.index(majority_element)]
        return majority_element, minority_lst
    else:
        return majority_element


def get_incar(job):
    data_dict = job["input/incar/data_dict"]
    return {
        key: value for key, value in zip(data_dict["Parameter"], data_dict["Value"])
    }


def get_sigma(job):
    return {"sigma": _get_value_from_incar(job=job, key="SIGMA")}


def get_ismear(job):
    return {"ismear": _get_value_from_incar(job=job, key="ISMEAR")}


def get_encut(job):
    return {"encut": _get_value_from_incar(job=job, key="ENCUT")}


def get_n_kpts(job):
    return {"n_kpts": eval(job["input/kpoints/data_dict"]["Value"][3].split()[0])}


def get_n_equ_kpts(job):
    return {"n_equ_kpts": len(job['output/generic/dft/bands/k_points'])}


def get_total_number_of_atoms(job):
    return {"Number_of_atoms": len(job["input/structure/indices"])}


def get_average_waves(job):
    _, weights, planewaves = job["output/outcar/irreducible_kpoints"]
    return {"avg. plane waves": sum(weights * planewaves) / sum(weights)}


def get_plane_waves(job):
    _, weights, planewaves = job["output/outcar/irreducible_kpoints"]
    return {"plane waves": sum(weights * planewaves)}


def get_ekin_error(job):
    return {"energy_tot_wo_kin_corr": job["output/outcar/kin_energy_error"] + job["output/generic/energy_tot"][-1]}


def get_volume(job):
    return {"volume": job["output/generic/volume"][-1]}


def get_volume_per_atom(job):
    return {"volume": job["output/generic/volume"][-1] / get_total_number_of_atoms(job=job)["Number_of_atoms"]}


def get_elements(job):
    species = job["input/structure/species"]
    indices_lst = job["input/structure/indices"]
    indices_set = set(indices_lst)
    count_lst = [indices_lst.tolist().count(ind) for ind in indices_set]
    main_element = species[count_lst.index(np.max(count_lst))]
    return {species[ind]: count_lst[ind] for ind in indices_set}


def get_convergence_check(job):
    try:
        conv = job.project.load(job.job_id).convergence_check()
    except:
        conv = None
    return {"Convergence": conv}


def get_number_of_species(job):
    return {"Number_of_species": len(job["output/structure/species"])}


def get_number_of_ionic_steps(job):
    return {"Number_of_ionic_steps": len(job["output/generic/energy_tot"])}


def get_number_of_final_electronic_steps(job):
    el_steps = job["output/generic/scf_energies"]
    if len(el_steps) != 0:
        return {"Number_of_final_electronic_steps": len(el_steps[-1])}
    else:
        return {"Number_of_final_electronic_steps": None}


def get_majority_species(job):
    indices_lst = job["input/structure/indices"].tolist()
    element_lst = job["input/structure/species"]
    majority_element, minority_lst = get_majority(
        [element_lst[ind] for ind in indices_lst], minority=True
    )
    return {"majority_element": majority_element, "minority_element_list": minority_lst}


def get_majority_crystal_structure(job):
    basis = Atoms().from_hdf(job["input"])
    majority_element = basis.get_majority_species()["symbol"]
    majority_index = [
        ind for ind, el in enumerate(basis) if el.symbol == majority_element
    ]
    type_list = list(basis[majority_index].analyse_pyscal_cna_adaptive(
        mode="str",
        ovito_compatibility=True
    ))
    return {"crystal_structure": get_majority(type_list, minority=False)}


def get_job_name(job):
    return {"job_name": job.job_name}


def get_job_id(job):
    return {"job_id": job.job_id}


def get_energy_tot_per_atom(job):
    return {"energy_tot": job["output/generic/energy_tot"][-1] / get_total_number_of_atoms(job=job)["Number_of_atoms"]}


def get_energy_tot(job):
    return {"energy_tot": job["output/generic/energy_tot"][-1]}


def get_energy_pot_per_atom(job):
    return {"energy_pot": job["output/generic/energy_pot"][-1] / get_total_number_of_atoms(job=job)["Number_of_atoms"]}


def get_energy_pot(job):
    return {"energy_pot": job["output/generic/energy_pot"][-1]}


def get_energy_free_per_atom(job):
    return {"energy_free": job["output/generic/dft/energy_free"][-1] / get_total_number_of_atoms(job=job)["Number_of_atoms"]}


def get_energy_free(job):
    return {"energy_free": job["output/generic/dft/energy_free"][-1]}


def get_energy_int_per_atom(job):
    return {"energy_int": job["output/generic/dft/energy_int"][-1] / get_total_number_of_atoms(job=job)["Number_of_atoms"]}


def get_energy_int(job):
    return {"energy_int": job["output/generic/dft/energy_int"][-1]}


def get_f_states(job):
    if "occ_matrix" in job["output/electronic_structure"].list_nodes():
        return {
            "f_states": job["output/electronic_structure/occ_matrix"].flatten().tolist()
        }
    elif "occupancy_matrix" in job["output/electronic_structure"].list_nodes():
        return {
            "f_states": job["output/electronic_structure/occupancy_matrix"]
            .flatten()
            .tolist()
        }
    else:
        print("get_f_states(): ", job.job_name, job.status)
        return {"f_states": [0.0]}


def get_e_band(job):
    if "occ_matrix" in job["output/electronic_structure"].list_nodes():
        f_occ = job["output/electronic_structure/occ_matrix"].flatten()
        ev_mat = job["output/electronic_structure/eig_matrix"].flatten()
    elif "occupancy_matrix" in job["output/electronic_structure"].list_nodes():
        f_occ = job["output/electronic_structure/occupancy_matrix"].flatten()
        ev_mat = job["output/electronic_structure/eigenvalue_matrix"].flatten()
    else:
        print("get_e_band(): ", job.job_name, job.status)
        f_occ = np.array([0.0])
        ev_mat = np.array([0.0])
    return {"e_band": np.sum(ev_mat * f_occ)}


def get_equilibrium_parameters(job):
    return {
        key: job["output/" + key]
        for key in [
            "equilibrium_energy",
            "equilibrium_b_prime",
            "equilibrium_bulk_modulus",
            "equilibrium_volume",
        ]
    }


def get_structure(job):
    atoms = pyiron_to_ase(job.load_object().get_structure())
    atoms_dict = {
        "symbols": atoms.get_chemical_symbols(),
        "positions": atoms.get_positions().tolist(),
        "cell": atoms.get_cell().tolist(),
        "pbc": atoms.get_pbc().tolist(),
        "celldisp": atoms.get_celldisp().tolist(),
    }
    if atoms.has("tags"):
        atoms_dict["tags"] = atoms.get_tags().tolist()
    if atoms.has("masses"):
        atoms_dict["masses"] = atoms.get_masses().tolist()
    if atoms.has("momenta"):
        atoms_dict["momenta"] = atoms.get_momenta().tolist()
    if atoms.has("initial_magmoms"):
        atoms_dict["magmoms"] = atoms.get_initial_magnetic_moments().tolist()
    if atoms.has("initial_charges"):
        atoms_dict["charges"] = atoms.get_initial_charges().tolist()
    if not atoms.__dict__["_calc"] == None:
        warnings.warn("Found calculator: " + str(atoms.__dict__["_calc"]))
    if not atoms.__dict__["_constraints"] == []:
        warnings.warn("Found constraint: " + str(atoms.__dict__["_constraints"]))
    return {"structure": json.dumps(atoms_dict)}


def get_forces(job):
    return {"forces": json.dumps(job["output/generic/forces"][-1].tolist())}


def get_magnetic_structure(job):
    basis = Atoms().from_hdf(job["input"])
    magmons = basis.get_initial_magnetic_moments()
    if all(magmons == None):
        return {"magnetic_structure": "non magnetic"}
    else:
        abs_sum_mag = sum(np.abs(magmons))
        sum_mag = sum(magmons)
        if abs_sum_mag == 0 and sum_mag == 0:
            return {"magnetic_structure": "non magnetic"}
        elif abs_sum_mag == np.abs(sum_mag):
            return {"magnetic_structure": "ferro-magnetic"}
        elif abs_sum_mag > 0 and sum_mag == 0:
            return {"magnetic_structure": "para-magnetic"}
        else:
            return {"magnetic_structure": "unknown"}


def get_e_conv_level(job):
    return {'el_conv': np.max(np.abs(
        job['output/generic/dft/scf_energy_free'][0] -
        job['output/generic/dft/scf_energy_free'][0][-1]
    )[-10:])}
