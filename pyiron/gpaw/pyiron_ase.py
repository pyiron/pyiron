# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from ase import Atoms
from ase.constraints import dict2constraint
import copy
import importlib
import numpy as np
from pyiron_base import InteractiveBase
from pyiron.atomistics.job.interactive import GenericInteractiveOutput

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


class AseJob(InteractiveBase):
    def __init__(self, project, job_name):
        super(AseJob, self).__init__(project, job_name)
        self.__name__ = "AseJob"
        self.__version__ = (
            None
        )  # Reset the version number to the executable is set automatically
        self._structure = None
        self.output = GenericInteractiveOutput(job=self)
        self._python_only_job = True

    @staticmethod
    def _cellfromdict(celldict):
        if Cell is not None:
            return Cell(**celldict)
        else:
            return celldict

    def _todict(self):
        atoms_dict = {
            "symbols": self._structure.get_chemical_symbols(),
            "positions": self._structure.get_positions(),
            "pbc": self._structure.get_pbc(),
            "celldisp": self._structure.get_celldisp(),
            "constraint": [c.todict() for c in self._structure.constraints],
            "info": copy.deepcopy(self._structure.info),
        }
        if Cell is not None:
            atoms_dict["cell"] = self._structure.get_cell().todict()
        else:
            atoms_dict["cell"] = self._structure.get_cell()
        if self._structure.has("tags"):
            atoms_dict["tags"] = self._structure.get_tags()
        if self._structure.has("masses"):
            atoms_dict["masses"] = self._structure.get_masses()
        if self._structure.has("momenta"):
            atoms_dict["momenta"] = self._structure.get_momenta()
        if self._structure.has("initial_magmoms"):
            atoms_dict["magmoms"] = self._structure.get_initial_magnetic_moments()
        if self._structure.has("initial_charges"):
            atoms_dict["charges"] = self._structure.get_initial_charges()
        if self._structure._calc is not None:
            calculator_dict = self._structure._calc.todict()
            calculator_dict["calculator_class"] = (
                str(self._structure._calc.__class__).replace("'", " ").split()[1]
            )
            calculator_dict["label"] = self._structure._calc.label
            atoms_dict["calculator"] = calculator_dict
        return atoms_dict

    def _fromdict(self, atoms_dict):
        atoms_dict_copy = copy.deepcopy(atoms_dict)
        if "calculator" in atoms_dict_copy.keys():
            calculator_dict = atoms_dict_copy["calculator"]
            calculator_class = calculator_dict["calculator_class"]
            del calculator_dict["calculator_class"]
            atoms_dict_copy["calculator"] = self._dict2calculator(
                calculator_class, calculator_dict
            )
        if "constraint" in atoms_dict_copy.keys():
            atoms_dict_copy["constraint"] = [
                dict2constraint(const_dict)
                for const_dict in atoms_dict_copy["constraint"]
            ]
        atoms_dict_copy["cell"] = self._cellfromdict(celldict=atoms_dict_copy["cell"])
        atoms = Atoms(**atoms_dict_copy)
        if atoms.calc is not None:
            atoms.calc.read(atoms.calc.label)
        return atoms

    @staticmethod
    def _dict2calculator(class_path, class_dict):
        module_loaded = importlib.import_module(".".join(class_path.split(".")[:-1]))
        module_class = getattr(module_loaded, class_path.split(".")[-1])
        return module_class(**class_dict)

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, basis):
        self._structure = basis

    def to_hdf(self, hdf=None, group_name=None):
        super(AseJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf_input:
            hdf_input["structure"] = self._todict()

    def from_hdf(self, hdf=None, group_name=None):
        super(AseJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf_input:
            self.structure = self._fromdict(hdf_input["structure"])

    def run_static(self):
        raise ValueError("ASEJob requires interactive mode.")

    def run_if_interactive(self):
        self.status.running = True
        self._structure.calc.set_label(self.working_directory + "/")
        self.structure.calc.calculate(self.structure)
        self.interactive_collect()

    def interactive_close(self):
        self.status.collect = True
        if (
            len(list(self.interactive_cache.keys())) > 0
            and len(self.interactive_cache[list(self.interactive_cache.keys())[0]]) != 0
        ):
            self.interactive_flush(path="interactive")
        with self.project_hdf5.open("output") as h5:
            if "interactive" in h5.list_groups():
                for key in h5["interactive"].list_nodes():
                    h5["generic/" + key] = h5["interactive/" + key]
        self.status.finished = True

    def interactive_forces_getter(self):
        return self._structure.get_forces()

    def interactive_energy_pot_getter(self):
        return self._structure.get_potential_energy()

    def interactive_energy_tot_getter(self):
        return self._structure.get_potential_energy()

    def interactive_indices_getter(self):
        element_lst = sorted(list(set(self._structure.get_chemical_symbols())))
        return np.array(
            [element_lst.index(el) for el in self._structure.get_chemical_symbols()]
        )

    def interactive_positions_getter(self):
        return self._structure.positions.copy()

    def interactive_steps_getter(self):
        return len(self.interactive_cache[list(self.interactive_cache.keys())[0]])

    def interactive_time_getter(self):
        return self.interactive_steps_getter()

    def interactive_volume_getter(self):
        return self._structure.get_volume()

    def interactive_cells_getter(self):
        return self._structure.cell.copy()


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
