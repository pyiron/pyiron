# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from ase import Atoms
from ase.constraints import dict2constraint
import copy
import importlib
import numpy as np
from pyiron.atomistics.job.interactive import GenericInteractive
from pyiron.atomistics.structure.atoms import pyiron_to_ase, Atoms as PAtoms

try:
    from ase.cell import Cell
except ImportError:
    Cell = None

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


def ase_structure_todict(structure):
    atoms_dict = {
        "symbols": structure.get_chemical_symbols(),
        "positions": structure.get_positions(),
        "pbc": structure.get_pbc(),
        "celldisp": structure.get_celldisp(),
        "constraint": [c.todict() for c in structure.constraints],
        "info": copy.deepcopy(structure.info),
    }
    if Cell is not None:
        atoms_dict["cell"] = structure.get_cell().todict()
    else:
        atoms_dict["cell"] = structure.get_cell()
    if structure.has("tags"):
        atoms_dict["tags"] = structure.get_tags()
    if structure.has("masses"):
        atoms_dict["masses"] = structure.get_masses()
    if structure.has("momenta"):
        atoms_dict["momenta"] = structure.get_momenta()
    if structure.has("initial_magmoms"):
        atoms_dict["magmoms"] = structure.get_initial_magnetic_moments()
    if structure.has("initial_charges"):
        atoms_dict["charges"] = structure.get_initial_charges()
    if structure.calc is not None:
        calculator_dict = structure.calc.todict()
        calculator_dict["calculator_class"] = (
            str(structure.calc.__class__).replace("'", " ").split()[1]
        )
        calculator_dict["label"] = structure.calc.label
        atoms_dict["calculator"] = calculator_dict
    return atoms_dict


def ase_calculator_fromdict(class_path, class_dict):
    module_loaded = importlib.import_module(".".join(class_path.split(".")[:-1]))
    module_class = getattr(module_loaded, class_path.split(".")[-1])
    return module_class(**class_dict)


def ase_structure_fromdict(atoms_dict):
    def cell_fromdict(celldict):
        celldict.pop("pbc", None)
        if Cell is not None:
            return Cell(**celldict)
        else:
            return celldict

    atoms_dict_copy = copy.deepcopy(atoms_dict)
    if "calculator" in atoms_dict_copy.keys():
        calculator_dict = atoms_dict_copy["calculator"]
        calculator_class = calculator_dict["calculator_class"]
        del calculator_dict["calculator_class"]
        atoms_dict_copy["calculator"] = ase_calculator_fromdict(
            calculator_class, calculator_dict
        )
    if "constraint" in atoms_dict_copy.keys():
        atoms_dict_copy["constraint"] = [
            dict2constraint(const_dict)
            for const_dict in atoms_dict_copy["constraint"]
        ]
    atoms_dict_copy["cell"] = cell_fromdict(celldict=atoms_dict_copy["cell"])
    atoms = Atoms(**atoms_dict_copy)
    if atoms.calc is not None:
        atoms.calc.read(atoms.calc.label)
    return atoms


class AseJob(GenericInteractive):
    def __init__(self, project, job_name):
        super(AseJob, self).__init__(project, job_name)
        self.__name__ = "AseJob"
        self.__version__ = (
            None
        )  # Reset the version number to the executable is set automatically

    @property
    def structure(self):
        return GenericInteractive.structure.fget(self)

    @structure.setter
    def structure(self, structure):
        if isinstance(structure, PAtoms):
            structure = pyiron_to_ase(structure)
        GenericInteractive.structure.fset(self, structure)

    def to_hdf(self, hdf=None, group_name=None):
        super(AseJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf_input:
            hdf_input["structure"] = ase_structure_todict(self._structure)

    def from_hdf(self, hdf=None, group_name=None):
        super(AseJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf_input:
            self.structure = ase_structure_fromdict(hdf_input["structure"])

    def run_static(self):
        pre_run_mode = self.server.run_mode
        self.server.run_mode.interactive = True
        self.run_if_interactive()
        self.interactive_close()
        self.server.run_mode = pre_run_mode

    def run_if_interactive(self):
        if self.structure.calc is None:
            self.set_calculator()
        super(AseJob, self).run_if_interactive()
        self.interactive_collect()

    def set_calculator(self):
        raise NotImplementedError(
            "The _set_calculator function is not implemented for this code."
        )

    def interactive_structure_setter(self, structure):
        self.structure.calc.calculate(structure)

    def interactive_initialize_interface(self):
        self.status.running = True
        self._structure.calc.set_label(self.working_directory + "/")
        self._interactive_library = True

    def interactive_close(self):
        if self.interactive_is_activated():
            super(AseJob, self).interactive_close()
            with self.project_hdf5.open("output") as h5:
                if "interactive" in h5.list_groups():
                    for key in h5["interactive"].list_nodes():
                        h5["generic/" + key] = h5["interactive/" + key]

    def interactive_forces_getter(self):
        return self.structure.get_forces()

    def interactive_pressures_getter(self):
        return -self.structure.get_stress(voigt=False)

    def interactive_energy_pot_getter(self):
        return self.structure.get_potential_energy()

    def interactive_energy_tot_getter(self):
        return self.structure.get_potential_energy()

    def interactive_indices_getter(self):
        element_lst = sorted(list(set(self.structure.get_chemical_symbols())))
        return np.array(
            [element_lst.index(el) for el in self.structure.get_chemical_symbols()]
        )

    def interactive_positions_getter(self):
        return self.structure.positions.copy()

    def interactive_steps_getter(self):
        return len(self.interactive_cache[list(self.interactive_cache.keys())[0]])

    def interactive_time_getter(self):
        return self.interactive_steps_getter()

    def interactive_volume_getter(self):
        return self.structure.get_volume()

    def interactive_cells_getter(self):
        return self.structure.cell.copy()

    def write_input(self):
        pass

    def collect_output(self):
        pass

    def run_if_scheduler(self):
        self._create_working_directory()
        super(AseJob, self).run_if_scheduler()

    def interactive_index_organizer(self):
        index_merge_lst = self._interactive_species_lst.tolist() + list(
            np.unique(self._structure_current.get_chemical_symbols())
        )
        el_lst = sorted(set(index_merge_lst), key=index_merge_lst.index)
        current_structure_index = [
            el_lst.index(el)
            for el in self._structure_current.get_chemical_symbols()
        ]
        previous_structure_index = [
            el_lst.index(el)
            for el in self._structure_previous.get_chemical_symbols()
        ]
        if not np.array_equal(
            np.array(current_structure_index),
            np.array(previous_structure_index),
        ):
            self._logger.debug("Generic library: indices changed!")
            self.interactive_indices_setter(self._structure_current.indices)

    def get_structure(self, iteration_step=-1, wrap_atoms=True):
        """
        Gets the structure from a given iteration step of the simulation (MD/ionic relaxation). For static calculations
        there is only one ionic iteration step
        Args:
            iteration_step (int): Step for which the structure is requested
            wrap_atoms (bool):

        Returns:
            atomistics.structure.atoms.Atoms object
        """
        if (
            self.server.run_mode.interactive
            or self.server.run_mode.interactive_non_modal
        ):
            # Warning: We only copy symbols, positions and cell information - no tags.
            if self.output.indices is not None and len(self.output.indices) != 0:
                indices = self.output.indices[iteration_step]
            else:
                return None
            if len(self._interactive_species_lst) == 0:
                el_lst = list(
                    np.unique(self._structure_current.get_chemical_symbols())
                )
            else:
                el_lst = self._interactive_species_lst.tolist()
            if indices is not None:
                if wrap_atoms:
                    positions = self.output.positions[iteration_step]
                else:
                    if len(self.output.unwrapped_positions) > max([iteration_step, 0]):
                        positions = self.output.unwrapped_positions[iteration_step]
                    else:
                        positions = (
                            self.output.positions[iteration_step]
                            + self.output.total_displacements[iteration_step]
                        )
                atoms = Atoms(
                    symbols=np.array([el_lst[el] for el in indices]),
                    positions=positions,
                    cell=self.output.cells[iteration_step],
                    pbc=self.structure.pbc,
                )
                # Update indicies to match the indicies in the cache.
                atoms.set_species([self._periodic_table.element(el) for el in el_lst])
                atoms.indices = indices
                if wrap_atoms:
                    atoms = atoms.center_coordinates_in_unit_cell()
                return atoms
            else:
                return None
        else:
            if (
                self.get("output/generic/cells") is not None
                and len(self.get("output/generic/cells")) != 0
            ):
                return super(GenericInteractive, self).get_structure(
                    iteration_step=iteration_step, wrap_atoms=wrap_atoms
                )
            else:
                return None


class AseAdapter(object):
    def __init__(self, ham, fast_mode=False):
        self._ham = ham
        self._fast_mode = fast_mode
        if self._ham.server.run_mode.interactive and fast_mode:
            self.interactive_cache = {
                "velocities": [],
                "energy_kin": [],
                "momenta": [],
                "positions": [],
                "energy_tot": [],
                "energy_pot": []
            }
            self._ham.run()
            self._ham.interactive_cache = {}
        elif self._ham.server.run_mode.interactive:
            self.interactive_cache = {
                "velocities": [],
                "energy_kin": [],
                "momenta": []
            }
        self.constraints = []
        try:
            self.arrays = {
                "positions": self._ham.structure.positions.copy(),
                "numbers": self._ham.structure.numbers,
            }
        except AttributeError:
            self.arrays = {
                "positions": self._ham.structure.positions.copy(),
                "numbers": self._ham.structure.get_atomic_numbers(),
            }

    @property
    def communicator(self):
        return None

    def get_masses(self):
        return np.array(self._ham.structure.get_masses())

    def get_positions(self):
        return self.arrays["positions"]

    def set_positions(self, positions):
        self.arrays["positions"] = positions

    def get_forces(self, md=True):
        if self._fast_mode:
            self._ham.interactive_positions_setter(self.arrays["positions"])
            self.interactive_cache["positions"].append(self.arrays["positions"])
            self._ham.interactive_execute()
            self.interactive_cache["energy_pot"].append(self._ham.interactive_energy_pot_getter())
            return np.array(self._ham.interactive_forces_getter())
        else:
            self._ham.structure.positions = self.arrays["positions"]
            if self._ham.server.run_mode.interactive:
                self._ham.run()
            else:
                self._ham.run(delete_existing_job=True)
            return self._ham.output.forces[-1]

    def interactive_close(self):
        self._ham.interactive_store_in_cache(
            "velocities", self.interactive_cache["velocities"]
        )
        self._ham.interactive_store_in_cache(
            "energy_kin", self.interactive_cache["energy_kin"]
        )
        if self._fast_mode:
            self._ham.interactive_store_in_cache(
                "positions", self.interactive_cache["positions"]
            )
            self._ham.interactive_store_in_cache(
                "energy_pot", self.interactive_cache["energy_pot"][::2]
            )
            self._ham.interactive_store_in_cache(
                "energy_tot",
                (
                    np.array(self.interactive_cache["energy_pot"][::2])
                    + np.array(self.interactive_cache["energy_kin"])
                ).tolist(),
            )
        else:
            self._ham.interactive_store_in_cache(
                "energy_tot",
                (
                    np.array(self._ham.output.energy_pot)[::2]
                    + np.array(self.interactive_cache["energy_kin"])
                ).tolist(),
            )
        self._ham.interactive_close()

    def get_number_of_atoms(self):
        return self._ham.structure.get_number_of_atoms()

    # ASE functions
    def get_kinetic_energy(self):
        """Get the kinetic energy."""
        momenta = self.arrays.get("momenta")
        if momenta is None:
            return 0.0
        return 0.5 * np.vdot(momenta, self.get_velocities())

    def set_momenta(self, momenta, apply_constraint=True):
        """Set momenta."""
        if apply_constraint and len(self.constraints) > 0 and momenta is not None:
            momenta = np.array(momenta)  # modify a copy
            for constraint in self.constraints:
                if hasattr(constraint, "adjust_momenta"):
                    constraint.adjust_momenta(self, momenta)
        self.set_array("momenta", momenta, float, (3,))
        self.interactive_cache["velocities"].append(self.get_velocities())
        self.interactive_cache["energy_kin"].append(self.get_kinetic_energy())

    def set_velocities(self, velocities):
        """Set the momenta by specifying the velocities."""
        self.set_momenta(self.get_masses()[:, np.newaxis] * velocities)

    def get_momenta(self):
        """Get array of momenta."""
        if "momenta" in self.arrays:
            return self.arrays["momenta"].copy()
        else:
            return np.zeros((len(self), 3))

    def set_array(self, name, a, dtype=None, shape=None):
        """Update array.

        If *shape* is not *None*, the shape of *a* will be checked.
        If *a* is *None*, then the array is deleted."""

        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name, a, dtype, shape)
        else:
            if a is None:
                del self.arrays[name]
            else:
                a = np.asarray(a)
                if a.shape != b.shape:
                    raise ValueError(
                        "Array has wrong shape %s != %s." % (a.shape, b.shape)
                    )
                b[:] = a

    def get_angular_momentum(self):
        """Get total angular momentum with respect to the center of mass."""
        com = self.get_center_of_mass()
        positions = self.get_positions()
        positions -= com  # translate center of mass to origin
        return np.cross(positions, self.get_momenta()).sum(0)

    def new_array(self, name, a, dtype=None, shape=None):
        """Add new array.

        If *shape* is not *None*, the shape of *a* will be checked."""

        if dtype is not None:
            a = np.array(a, dtype, order="C")
            if len(a) == 0 and shape is not None:
                a.shape = (-1,) + shape
        else:
            if not a.flags["C_CONTIGUOUS"]:
                a = np.ascontiguousarray(a)
            else:
                a = a.copy()

        if name in self.arrays:
            raise RuntimeError

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError("Array has wrong length: %d != %d." % (len(a), len(b)))
            break

        if shape is not None and a.shape[1:] != shape:
            raise ValueError(
                "Array has wrong shape %s != %s." % (a.shape, (a.shape[0:1] + shape))
            )

        self.arrays[name] = a

    def has(self, name):
        """Check for existence of array.

        name must be one of: 'tags', 'momenta', 'masses', 'initial_magmoms',
        'initial_charges'."""
        # XXX extend has to calculator properties
        return name in self.arrays

    def get_center_of_mass(self, scaled=False):
        """Get the center of mass.

        If scaled=True the center of mass in scaled coordinates
        is returned."""
        m = self.get_masses()
        com = np.dot(m, self.arrays["positions"]) / m.sum()
        if scaled:
            if self._fast_mode:
                return np.linalg.solve(self._ham.structure.cells[-1].T, com)
            else:
                return np.linalg.solve(self._ham.output.cells[-1].T, com)
        else:
            return com

    def get_velocities(self):
        """Get array of velocities."""
        momenta = self.arrays.get("momenta")
        if momenta is None:
            return None
        m = self.get_masses()
        # m = self.arrays.get('masses')
        # if m is None:
        #     m = atomic_masses[self.arrays['numbers']]
        return momenta / m.reshape(-1, 1)

    def __len__(self):
        return len(self._ham.structure)
