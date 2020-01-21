# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
import hashlib
import numpy as np
import warnings
from pyiron.base.generic.parameters import GenericParameters
import scipy.constants as spc

__author__ = "Joerg Neugebauer, Sudarsan Surendralal, Jan Janssen"
__copyright__ = (
    "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Sudarsan Surendralal"
__email__ = "surendralal@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


# Conversion factors for transfroming pyiron units to Lammps units (alphabetical)
AMU_TO_G = spc.atomic_mass * spc.kilo
AMU_TO_KG = spc.atomic_mass
ANG_PER_FS_TO_ANG_PER_PS = spc.pico / spc.femto
ANG_PER_FS_TO_BOHR_PER_FS = spc.angstrom / spc.physical_constants['Bohr radius'][0]
ANG_PER_FS_TO_CM_PER_S = (spc.angstrom / spc.femto) / spc.centi
ANG_PER_FS_TO_M_PER_S = spc.angstrom / spc.femto
ANG_TO_BOHR = spc.angstrom / spc.physical_constants['Bohr radius'][0]
ANG_TO_CM = spc.angstrom / spc.centi
ANG_TO_M = spc.angstrom
EL_TO_COUL = spc.elementary_charge
EV_PER_ANG_TO_DYNE = (spc.electron_volt / spc.angstrom) / spc.dyne
EV_PER_ANG_TO_HA_PER_BOHR = spc.physical_constants["electron volt-hartree relationship"][0] * \
                            spc.physical_constants['Bohr radius'][0] / spc.angstrom
EV_PER_ANG_TO_KCAL_PER_MOL_ANG = spc.eV / (spc.kilo * spc.calorie / spc.N_A)
EV_PER_ANG_TO_N = spc.electron_volt / spc.angstrom
EV_TO_ERG = spc.electron_volt / spc.erg
EV_TO_HA = spc.physical_constants["electron volt-hartree relationship"][0]
EV_TO_J = spc.electron_volt
EV_TO_KCAL_PER_MOL = spc.eV / (spc.kilo * spc.calorie / spc.N_A)
FS_TO_PS = spc.femto / spc.pico
FS_TO_S = spc.femto
GPA_TO_ATM = spc.giga / spc.atm
GPA_TO_BAR = spc.giga / spc.bar
GPA_TO_BARYE = spc.giga / (spc.micro * spc.bar)  # "barye" = 1e-6 bar
GPA_TO_PA = spc.giga

# Conversions for most of the Lammps units to Pyiron units
# Lammps units source doc: https://lammps.sandia.gov/doc/units.html
# Pyrion units source doc: https://pyiron.github.io/source/faq.html
# At time of writing, not all these conversion factors are used, but may be helpful later.
LAMMPS_UNIT_CONVERSIONS = {
    "metal": {
        "mass": 1.,
        "distance": 1.,
        "time": FS_TO_PS,
        "energy": 1.,
        "velocity": ANG_PER_FS_TO_ANG_PER_PS,
        "force": 1.,
        "temperature": 1.,
        "pressure": GPA_TO_BAR,
        "charge": 1.
    },
    "si": {
        "mass": AMU_TO_KG,
        "distance": ANG_TO_M,
        "time": FS_TO_S,
        "energy": EV_TO_J,
        "velocity": ANG_PER_FS_TO_M_PER_S,
        "force": EV_PER_ANG_TO_N,
        "temperature": 1.,
        "pressure": GPA_TO_PA,
        "charge": EL_TO_COUL
    },
    "cgs": {
        "mass": AMU_TO_G,
        "distance": ANG_TO_CM,
        "time": FS_TO_S,
        "energy": EV_TO_ERG,
        "velocity": ANG_PER_FS_TO_CM_PER_S,
        "force": EV_PER_ANG_TO_DYNE,
        "temperature": 1.,
        "pressure": GPA_TO_BARYE,
        "charge": 4.8032044e-10  # In statCoulombs, but these are deprecated and thus not in scipt.constants
    },
    "real": {
        "mass": 1.,
        "distance": 1.,
        "time": 1.,
        "energy": EV_TO_KCAL_PER_MOL,
        "velocity": 1.,
        "force": EV_PER_ANG_TO_KCAL_PER_MOL_ANG,
        "temperature": 1.,
        "pressure": GPA_TO_ATM,
        "charge": 1.
    },
    "electron": {
        "mass": 1.,
        "distance": ANG_TO_BOHR,
        "time": 1.,
        "energy": EV_TO_HA,
        "velocity": ANG_PER_FS_TO_BOHR_PER_FS,
        "force": EV_PER_ANG_TO_HA_PER_BOHR,
        "temperature": 1.,
        "pressure": GPA_TO_PA,
        "charge": 1.
    },
}


class LammpsControl(GenericParameters):
    def __init__(self, input_file_name=None, **qwargs):
        super(LammpsControl, self).__init__(
            input_file_name=input_file_name, table_name="control_inp", comment_char="#"
        )
        self._init_block_dict()

    @property
    def dataset(self):
        return self._dataset

    def _init_block_dict(self):
        block_dict = OrderedDict()
        block_dict["read_restart"] = "read_restart"
        block_dict["structure"] = (
            "units",
            "dimension",
            "boundary",
            "atom_style",
            "read_data",
            "include",
        )
        block_dict["potential"] = (
            "group",
            "set",
            "pair_style",
            "pair_coeff",
            "bond_style",
            "bond_coeff",
            "angle_style",
            "angle_coeff",
            "kspace_style",
            "neighbor",
        )
        block_dict["compute"] = (
            "compute",
            "fix",
            "variable",
            "fix_modify",
        )
        block_dict["setup"] = (
            "min_style",
            "min_modify",
            "neigh_modify",
            "velocity",
            "timestep",
            "dielectric",
            "reset_timestep",
        )
        block_dict["control"] = (
            "dump",
            "dump_modify",
            "thermo_style",
            "thermo_modify",
            "thermo",
        )
        block_dict["run"] = ("run", "minimize")
        block_dict["write_restart"] = "write_restart"
        self._block_dict = block_dict

    def load_default(self, file_content=None):
        if file_content is None:
            file_content = (
                "units                 metal\n"
                + "dimension             3\n"
                + "boundary              p p p\n"
                + "atom_style            atomic\n"
                + "read_data             structure.inp\n"
                + "include               potential.inp\n"
                + "fix___ensemble        all nve\n"
                + "variable___dumptime   equal 100\n"
                + "variable___thermotime equal 100\n"
                + "dump___1              all custom ${dumptime} dump.out id type xsu ysu zsu fx fy fz vx vy vz\n"
                + "dump_modify___1       sort id format line \"%d %d %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\"\n"
                + "thermo_style          custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol\n"
                + "thermo_modify         format float %20.15g\n"
                + "thermo                ${thermotime}\n"
                + "run                   0\n"
            )
        self.load_string(file_content)

    def calc_minimize(
        self,
        e_tol=0.0,
        f_tol=1e-4,
        max_iter=100000,
        pressure=None,
        n_print=100,
        style='cg'
    ):
        """
        Sets parameters required for minimization.

        Args:
            e_tol (float): If the magnitude of difference between energies of two consecutive steps is lower than or
                equal to `e_tol`, the minimisation terminates. (Default is 0.0 eV.)
            f_tol (float): If the magnitude of the global force vector at a step is lower than or equal to `f_tol`, the
                minimisation terminates. (Default is 1e-4 eV/angstrom.)
            max_iter (int): Maximum number of minimisation steps to carry out. If the minimisation converges before
                `max_iter` steps, terminate at the converged step. If the minimisation does not converge up to
                `max_iter` steps, terminate at the `max_iter` step. (Default is 100000.)
            pressure (float/list/tuple/numpy.ndarray): Pressure in GPa at which minimisation is to be carried out. If
                None, isochoric (constant volume) condition will be used. If a float, cell shape changes are allowed.
                in each of the three primary axes, but no shearing occurs. If array-like, (up to) the first six
                elements will be interpreted as the x, y, z, xy, xz, and yz components of the pressure tensor,
                respectively. If this component-wise mode is used, cell changes corresponding to any one element of the
                stress tensor can be selectively disabled by setting this element to None. (Default is None, run
                isochorically.)
            n_print (int): Write (dump or print) to the output file every n steps (Default: 100)
            style ('cg'/'sd'/other values from Lammps docs): The style of the numeric minimization, either conjugate
                gradient, steepest descent, or other keys permissible from the Lammps docs on 'min_style'. (Default
                is 'cg' -- conjugate gradient.)
        """
        # This docstring is a source for the calc_minimize method in pyiron.lammps.base.LammpsBase.calc_minimize and
        # pyiron.lammps.interactive.LammpsInteractive.calc_minimize -- Please ensure that changes to signature or
        # defaults stay consistent!

        max_evaluations = 100 * max_iter

        if self["units"] not in LAMMPS_UNIT_CONVERSIONS.keys():
            raise NotImplementedError
        energy_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["energy"]
        force_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["force"]
        pressure_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["pressure"]

        e_tol *= energy_units
        f_tol *= force_units

        if pressure is not None:
            if type(pressure) == float or type(pressure) == int:
                pressure = pressure * np.ones(3)
            str_press = ""
            for press, str_axis in zip(pressure, [" x ", " y ", " z ", " xy ", " xz ", " yz "][:len(pressure)]):
                if press is not None:
                    str_press += str_axis + str(press * pressure_units)
            if len(str_press) == 0:
                raise ValueError("Pressure values cannot all be None")
            elif len(str_press) > 1:
                str_press += " couple none"
            self.set(fix___ensemble=r"all box/relax" + str_press)
        self.remove_keys(["fix___nve"])
        self.set(min_style=style)
        self.set(
            minimize=str(e_tol)
            + " "
            + str(f_tol)
            + " "
            + str(int(max_iter))
            + " "
            + str(int(max_evaluations))
        )
        self.remove_keys(["run", "velocity"])
        self.modify(variable___dumptime="equal " + str(n_print), thermo=n_print)
        self.modify(variable___thermotime="equal " + str(n_print), thermo=n_print)

    def calc_static(self):
        self.set(run="0")
        self.remove_keys(["minimize", "velocity"])

    def set_initial_velocity(
        self,
        temperature,
        seed=None,
        gaussian=False,
        append_value=False,
        zero_lin_momentum=True,
        zero_rot_momentum=True,
        job_name="",
    ):
        """
        Create initial velocities via velocity all create. More information can be found on LAMMPS website:
        https://lammps.sandia.gov/doc/velocity.html

        Args:
            temperature: (int or float)
            seed: (int) Seed for the initial random number generator
            gaussian: (True/False) Create velocity according to the Gaussian distribution (otherwise uniform)
            append_value: (True/False) Add the velocity values to the current velocities (probably not functional now)
            zero_lin_momentum: (True/False) Cancel the total linear momentum
            zero_rot_momentum: (True/False) Cancel the total angular momentum
            job_name: (str) job name to generate seed
        """

        if seed is None:
            seed = self.generate_seed_from_job(job_name=job_name, seed=1)
        arg = ""
        if gaussian:
            arg = " dist gaussian"
        if append_value:
            arg += " sum yes"
        if not zero_lin_momentum:
            arg += " mom no"
        if not zero_rot_momentum:
            arg += " rot no"
        self.modify(
            velocity="all create " + str(temperature) + " " + str(seed) + arg,
            append_if_not_present=True,
        )

    @staticmethod
    def generate_seed_from_job(job_name="", seed=0):
        """
        Generate a unique seed from the job name.

        Args:
            job_name (str): job_name of the current job to generate the seed
            seed (int): integer to access part of the seed

        Returns:
            int: random seed generated based on the hash
        """
        return int(
            str(int(hashlib.sha256(job_name.encode()).hexdigest(), 16))[
                5 * seed : 5 * seed + 5
            ]
        )

    def calc_md(
        self,
        temperature=None,
        pressure=None,
        n_ionic_steps=1000,
        time_step=1.0,
        n_print=100,
        temperature_damping_timescale=100.0,
        pressure_damping_timescale=1000.0,
        seed=None,
        tloop=None,
        initial_temperature=None,
        langevin=False,
        delta_temp=None,
        delta_press=None,
        job_name="",
    ):
        """
        Set an MD calculation within LAMMPS. Nosé Hoover is used by default.

        Args:
            temperature (None/float): Target temperature. If set to None, an NVE calculation is performed.
                                      It is required when the pressure is set or langevin is set
            pressure (None/float/numpy.ndarray/list): Target pressure. If set to None, an NVE or an NVT calculation is
                performed. A length-3 list or array may be given to specify x-, y- and z-components individually. In
                this case, floats and `None` may be mixed to allow relaxation only in particular directions.
            n_ionic_steps (int): Number of ionic steps
            time_step (float): Step size in fs between two steps.
            n_print (int):  Print frequency
            temperature_damping_timescale (float): The time associated with the thermostat adjusting the temperature.
                                                   (In fs. After rescaling to appropriate time units, is equivalent to
                                                   Lammps' `Tdamp`.)
            pressure_damping_timescale (float): The time associated with the barostat adjusting the temperature.
                                                (In fs. After rescaling to appropriate time units, is equivalent to
                                                Lammps' `Pdamp`.)
            seed (int):  Seed for the random number generation (required for the velocity creation)
            tloop:
            initial_temperature (None/float):  Initial temperature according to which the initial velocity field
                                               is created. If None, the initial temperature will be twice the target
                                               temperature (which would go immediately down to the target temperature
                                               as described in equipartition theorem). If 0, the velocity field is not
                                               initialized (in which case  the initial velocity given in structure will
                                               be used). If any other number is given, this value is going to be used
                                               for the initial temperature.
            langevin (bool): (True or False) Activate Langevin dynamics
            delta_temp (float): Thermostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
            delta_press (float): Barostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
            job_name (str): Job name of the job to generate a unique random seed.
        """
        if self["units"] not in LAMMPS_UNIT_CONVERSIONS.keys():
            raise NotImplementedError
        time_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["time"]
        temperature_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["temperature"]
        pressure_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["pressure"]

        # Transform time
        if time_step is not None:
            try:
                self["timestep"] = time_step * time_units
            except KeyError:
                raise NotImplementedError()

        # Transform thermostat strength (time)
        if delta_temp is not None:
            warnings.warn(
                "WARNING: `delta_temp` is deprecated, please use `temperature_damping_timescale`."
            )
            temperature_damping_timescale = delta_temp
        else:
            temperature_damping_timescale *= time_units

        # Transform barostat strength (time)
        if delta_press is not None:
            warnings.warn(
                "WARNING: `delta_press` is deprecated, please use `pressure_damping_timescale`."
            )
            pressure_damping_timescale = delta_press
        else:
            pressure_damping_timescale *= time_units

        # Transform temperature
        if temperature is not None:
            temperature *= temperature_units

        # Apply initial overheating (default uses the theorem of equipartition of energy between KE and PE)
        if initial_temperature is None and temperature is not None:
            initial_temperature = 2 * temperature

        if seed is None:
            seed = self.generate_seed_from_job(job_name=job_name)

        # Set thermodynamic ensemble
        if pressure is not None:  # NPT
            if not hasattr(pressure, "__len__"):
                pressure = pressure * np.ones(3)
            else:
                pressure = np.array(pressure, dtype=float)

            not_none_mask = [p is not None for p in pressure]
            if not np.any(not_none_mask):
                raise ValueError("Pressure cannot be three times None")

            if len(pressure) > 6:
                raise ValueError("Pressure must be a float or a vector with length <= 6")

            if temperature is None or temperature == 0.0:
                raise ValueError("Target temperature for fix nvt/npt/nph cannot be 0")

            pressure[not_none_mask] *= pressure_units

            pressure_string = ""
            for coord, value in zip(["x", "y", "z", "xy", "xz", "yz"][:len(pressure)], pressure):
                if value is not None:
                    pressure_string += " {0} {1} {1} {2}".format(
                        coord, str(value), str(pressure_damping_timescale)
                    )

            if langevin:  # NPT(Langevin)
                fix_ensemble_str = "all nph" + pressure_string
                self.modify(
                    fix___langevin="all langevin {0} {1} {2} {3} zero yes".format(
                        str(temperature),
                        str(temperature),
                        str(temperature_damping_timescale),
                        str(seed),
                    ),
                    append_if_not_present=True,
                )
            else:  # NPT(Nose-Hoover)
                fix_ensemble_str = "all npt temp {0} {1} {2}".format(
                    str(temperature),
                    str(temperature),
                    str(temperature_damping_timescale),
                )
                fix_ensemble_str += pressure_string
        elif temperature is not None:  # NVT
            if temperature == 0.0:
                raise ValueError("Target temperature for fix nvt/npt/nph cannot be 0.0")

            if langevin:  # NVT(Langevin)
                fix_ensemble_str = "all nve"
                self.modify(
                    fix___langevin="all langevin {0} {1} {2} {3} zero yes".format(
                        str(temperature),
                        str(temperature),
                        str(temperature_damping_timescale),
                        str(seed),
                    ),
                    append_if_not_present=True,
                )
            else:  # NVT(Nose-Hoover)
                fix_ensemble_str = "all nvt temp {0} {1} {2}".format(
                    str(temperature),
                    str(temperature),
                    str(temperature_damping_timescale),
                )
        else:  # NVE
            if langevin:
                warnings.warn("Temperature not set; Langevin ignored.")
            fix_ensemble_str = "all nve"
            initial_temperature = 0

        if tloop is not None:
            fix_ensemble_str += " tloop " + str(tloop)

        self.remove_keys(["minimize"])
        self.modify(
            fix___ensemble=fix_ensemble_str,
            variable___dumptime=" equal {} ".format(n_print),
            thermo=int(n_print),
            run=int(n_ionic_steps),
            append_if_not_present=True,
        )
        self.modify(
            fix___ensemble=fix_ensemble_str,
            variable___thermotime=" equal {} ".format(n_print),
            thermo=int(n_print),
            run=int(n_ionic_steps),
            append_if_not_present=True,
        )

        if initial_temperature > 0:
            self.set_initial_velocity(
                temperature=initial_temperature,
                seed=seed,
                gaussian=True,
                job_name=job_name
            )

    def calc_vcsgc(
        self,
        mu,
        ordered_element_list,
        target_concentration=None,
        kappa=1000.,
        mc_step_interval=100,
        swap_fraction=0.1,
        temperature_mc=None,
        window_size=None,
        window_moves=None,
        temperature=None,
        pressure=None,
        n_ionic_steps=1000,
        time_step=1.0,
        n_print=100,
        temperature_damping_timescale=100.0,
        pressure_damping_timescale=1000.0,
        seed=None,
        initial_temperature=None,
        langevin=False,
        job_name="",
    ):
        """
        Run variance-constrained semi-grand-canonical MD/MC for a binary system. In addition to VC-SGC arguments, all
        arguments for a regular MD calculation are also accepted.

        https://vcsgc-lammps.materialsmodeling.org

        Note:
            For easy visualization later (with `get_structure`), it is highly recommended that the initial structure
            contain at least one atom of each species.

        Warning:
            Assumes the units are metal, otherwise units for the constraints may be off.

        Args:
            mu (dict): A dictionary of chemical potentials, one for each element the potential treats, where the
                dictionary keys are just the chemical symbol. Note that only the *relative* chemical potentials are used
                here, such that the swap acceptance probability is influenced by the chemical potential difference
                between the two species (a more negative value increases the odds of swapping *to* that element.)
            ordered_element_list (list): A list of the chemical species symbols in the order they appear in the
                definition of the potential in the Lammps' input file.
            target_concentration: A dictionary of target simulation domain concentrations for each species *in the
                potential*. Dictionary keys should be the chemical symbol of the corresponding species, and the sum of
                all concentrations must be 1. (Default is None, which runs regular semi-grand-canonical MD/MC without
                any variance constraint.)
            kappa: Variance constraint for the MC. Larger value means a tighter adherence to the target concentrations.
                (Default is 1000.)
            mc_step_interval (int): How many steps of MD between each set of MC moves. (Default is 100.) Must divide the
                number of ionic steps evenly.
            swap_fraction (float): The fraction of atoms whose species is swapped at each MC phase. (Default is 0.1.)
            temperature_mc (float): The temperature for accepting MC steps. (Default is None, which uses the MD
                temperature.)
            window_size (float): The size of the sampling window for parallel calculations as a fraction of something
                unspecified in the VC-SGC docs, but it must lie between 0.5 and 1. (Default is None, window is
                determined automatically.)
            window_moves (int): The number of times the sampling window is moved during one MC cycle. (Default is None,
                number of moves is determined automatically.)
        """
        self.calc_md(
            temperature=temperature,
            pressure=pressure,
            n_ionic_steps=n_ionic_steps,
            time_step=time_step,
            n_print=n_print,
            temperature_damping_timescale=temperature_damping_timescale,
            pressure_damping_timescale=pressure_damping_timescale,
            seed=seed,
            tloop=None,
            initial_temperature=initial_temperature,
            langevin=langevin,
            job_name=job_name,
        )

        if self["units"] not in LAMMPS_UNIT_CONVERSIONS.keys():
            raise NotImplementedError
        temperature_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["temperature"]
        energy_units = LAMMPS_UNIT_CONVERSIONS[self["units"]]["energy"]

        if temperature_mc is None:
            if temperature is None:
                raise ValueError("If temperature is not given, temperature_mc must be.")
            temperature_mc = temperature * temperature_units

        if seed is None:
            seed = self.generate_seed_from_job(job_name=job_name)

        if set(mu.keys()) != set(ordered_element_list):
            raise ValueError("Exactly one chemical potential must be given for each element treated by the potential.")

        calibrating_el = ordered_element_list[0]
        # Apply the actual SGC string
        fix_vcsgc_str = "all sgcmc {0} {1} {2} {3} randseed {4}".format(
            str(mc_step_interval),
            str(swap_fraction),
            str(temperature_mc),
            str(" ".join(str((mu[el] - mu[calibrating_el]) * energy_units) for el in ordered_element_list[1:])),
            str(seed),
        )

        # Add VC to the SGC if target concentrations were provided
        if target_concentration is not None:
            if set(target_concentration.keys()) != set(ordered_element_list):
                raise ValueError("Exactly one target concentration must be given for each element treated by the "
                                 "potential.")

            if not np.isclose(np.sum([v for _, v in target_concentration.items()]), 1):
                raise ValueError("Target concentrations must sum to 1 but were {}.".format(
                    np.sum([v for _, v in target_concentration.items()])
                ))

            fix_vcsgc_str += " variance {0} {1}".format(
                str(kappa),
                str(" ".join([str(target_concentration[el]) for el in ordered_element_list[1:]]))
            )

        # Set optional windowing parameters
        if window_moves is not None:
            if int(window_moves) != window_moves or window_moves < 0:
                raise ValueError("Window moves must be a non-negative integer.")
            fix_vcsgc_str += " window_moves {0}".format(window_moves)
        if window_size is not None:
            if not 0.5 <= window_size <= 1.0:
                raise ValueError("Window size must be a fraction between 0.5 and 1")
            fix_vcsgc_str += " window_size {0}".format(window_size)

        self.modify(
            fix___vcsgc=fix_vcsgc_str,
            append_if_not_present=True
        )

    def measure_mean_value(self, key, every=1, repeat=None, name=None, atom=False):
        """
            Args:
                key (str): property to take an average value of (e.g. 'energy_pot' v.i.)
                every (int): number of steps there should be between two measurements
                repeat (int): number of measurements for each output (default: n_print/every)
                name (str): name to give in the output string (ignored if a pyiron predefined tag is used)

            Comments:
                Currently available keys: 'energy_pot', 'energy_tot', 'temperature', 'volume',
                                          'pressures', 'positions', 'forces, 'velocities'
                Future keys: 'cells'
        """

        if every<=0:
            raise AssertionError('every must be a positive integer')
        if repeat is None:
            self['variable___mean_repeat_times'] = 'equal round(${thermotime}/'+str(every)+')'
        else:
            self['variable___mean_repeat_times'] = 'equal {}'.format(repeat)
        if key=='energy_pot':
            self._measure_mean_value('energy_pot', 'pe', every)
        elif key=='energy_tot':
            self._measure_mean_value('energy_tot', 'etotal', every)
        elif key=='temperature':
            self._measure_mean_value('temperature', 'temp', every)
        elif key=='volume':
            self._measure_mean_value('volume', 'vol', every)
        elif key=='pressures':
            self._measure_mean_value('pressure', ['pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz'], every)
        elif key=='positions':
            self['compute___unwrap'] = 'all property/atom xu yu zu'
            self['fix___mean_positions'] = ('all ave/atom '
                                            +str(every)
                                            +' ${mean_repeat_times} ${thermotime} c_unwrap[*]'
                                           )
            self['dump___1'] = self['dump___1']+' '+' '.join(['f_mean_positions[{}]'.format(ii+1) for ii in range(3)])
            self['dump_modify___1'] = self['dump_modify___1'][:-1]+' '+' '.join(['%20.15g']*3)+'"'
        elif key=='forces':
            self._measure_mean_value('forces', ['fx', 'fy', 'fz'], every, atom=True)
        elif key=='velocities':
            self._measure_mean_value('velocities', ['vx', 'vy', 'vz'], every, atom=True)
        elif name is not None:
            if '**' in key:
                warnings.warn('** is replaced by ^ (as it is understood by LAMMPS)')
                key = key.replace('**', '^')
            self._measure_mean_value(name, key, every, atom)
        else:
            raise NotImplementedError(key+' is not implemented')

    def _measure_mean_value(self, key_pyiron, key_lmp, every, atom=False):
        """
            Args:
                key (str): property to take an average value of (e.g. 'energy_pot' v.i.)
                every (int): number of steps there should be between two measurements
                repeat (int): number of measurements for each output (default: n_print/every)
                freq (int): output frequency (default: n_print)

            Comments:
                Currently available keys: 'energy_pot', 'energy_tot', 'temperature', 'volume',
                                          'pressures', 'positions', 'forces, 'velocities'
                Future keys: 'cells'
        """
        if isinstance(key_lmp, str):
            self['variable___{}'.format(key_pyiron)] = 'equal {}'.format(key_lmp)
            self['fix___mean_{}'.format(key_pyiron)] = 'all ave/time '+str(every)+' ${mean_repeat_times} ${thermotime} v_'+str(key_pyiron)
            self['thermo_style'] = self['thermo_style']+' f_mean_{}'.format(key_pyiron)
        else:
            for ii, kk in enumerate(key_lmp):
                if atom is True:
                    self['variable___{}_{}'.format(key_pyiron, ii)] = 'atom {}'.format(key_lmp[ii])
                else:
                    self['variable___{}_{}'.format(key_pyiron, ii)] = 'equal {}'.format(key_lmp[ii])
            if atom is True:
                self['fix___mean_{}'.format(key_pyiron)] = ('all ave/atom '
                                                            +str(every)+
                                                            ' ${mean_repeat_times} ${thermotime} '
                                                            +str(' '.join(['v_{}_{}'.format(key_pyiron, ii) for ii in range(len(key_lmp))]))
                                                           )
                self['dump___1'] = self['dump___1']+' '+' '.join(['f_mean_{}[{}]'.format(key_pyiron, ii+1) for ii in range(len(key_lmp))])
                self['dump_modify___1'] = self['dump_modify___1'][:-1]+' '+' '.join(['%20.15g']*len(key_lmp))+'"'
            else:
                self['fix___mean_{}'.format(key_pyiron)] = 'all ave/time '+str(every)+' ${mean_repeat_times} ${thermotime} '+str(' '.join(['v_{}_{}'.format(key_pyiron, ii) for ii in range(len(key_lmp))]))
                self['thermo_style'] = self['thermo_style']+' '+' '.join(['f_mean_{}[{}]'.format(key_pyiron, ii+1) for ii in range(len(key_lmp))])

