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
        block_dict["compute"] = "compute"
        block_dict["control"] = (
            "fix___ensemble",
            "fix___langevin",
            "fix___vcsgc",
            "fix",
            "min_style",
            "min_modify",
            "neigh_modify",
            "velocity",
            "dump",
            "dump_modify",
            "thermo_style",
            "thermo_modify",
            "thermo",
            "timestep",
            "dielectric",
            "fix_modify",
            "reset_timestep",
        )
        block_dict["run"] = ("run", "minimize")
        block_dict["write_restart"] = "write_restart"
        self._block_dict = block_dict

    def load_default(self, file_content=None):
        if file_content is None:
            file_content = (
                "units               metal\n"
                + "dimension           3\n"
                + "boundary            p p p\n"
                + "atom_style          atomic\n"
                + "read_data           structure.inp\n"
                + "include             potential.inp\n"
                + "fix___ensemble      all nve\n"
                + "variable___dumptime equal 100\n"
                + "dump___1            all custom ${dumptime} dump.out id type xsu ysu zsu fx fy fz\n"
                + 'dump_modify___1     sort id format line "%d %d %20.15g %20.15g %20.15g %20.15g %20.15g %20.15g"\n'
                + "thermo_style        custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol\n"
                + "thermo_modify       format float %20.15g\n"
                + "thermo              100\n"
                + "run                 0\n"
            )
        self.load_string(file_content)

    def calc_minimize(
        self, e_tol=0.0, f_tol=1e-8, max_iter=100000, pressure=None, n_print=100
    ):
        max_evaluations = 10 * max_iter
        if pressure is not None:
            if type(pressure) == float or type(pressure) == int:
                pressure = pressure * np.ones(3)
            str_press = ""
            for press, str_axis in zip(pressure, [" x ", " y ", " z "]):
                if press is not None:
                    str_press += str_axis + str(press * 1.0e4)
            if len(str_press) == 0:
                raise ValueError("Pressure values cannot be three times None")
            elif len(str_press) > 1:
                str_press += " couple none"
            self.set(fix___ensemble=r"all box/relax" + str_press)
        else:
            self.remove_keys(["fix___nve"])
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
            time_step (float): Step size between two steps. In fs if units==metal
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
        # Conversion factors for transfroming pyiron units to Lammps units
        fs_to_ps = spc.femto / spc.pico
        fs_to_s = spc.femto / 1.0
        GPa_to_bar = spc.giga * 1.0 / spc.bar
        GPa_to_Pa = spc.giga
        GPa_to_barye = spc.giga * 1.0 / (1.0e-6 * spc.bar)  # Lammps is in "barye"
        GPa_to_atm = spc.giga * 1.0 / spc.atm
        lammps_unit_conversions = {
            "metal": {"time": fs_to_ps, "pressure": GPa_to_bar},
            "si": {"time": fs_to_s, "pressure": GPa_to_Pa},
            "cgs": {"time": fs_to_s, "pressure": GPa_to_barye},
            "real": {"time": 1, "pressure": GPa_to_atm},
            "electron": {"time": 1, "pressure": GPa_to_Pa},
        }
        time_units = lammps_unit_conversions[self["units"]]["time"]
        pressure_units = lammps_unit_conversions[self["units"]]["pressure"]
        # No need for temperature conversion; pyiron and all available Lammps units are both in Kelvin
        # (well, except unitless Lennard-Jones units...)
        if self["units"] == "lj":
            raise NotImplementedError

        # Transform time
        if time_step is not None:
            try:
                self["timestep"] = time_step * time_units
            except KeyError:
                raise NotImplementedError()

        # Transform thermostat strength
        if delta_temp is not None:
            warnings.warn(
                "WARNING: `delta_temp` is deprecated, please use `temperature_damping_timescale`."
            )
            temperature_damping_timescale = delta_temp
        else:
            temperature_damping_timescale *= time_units

        # Transform barostat strength
        if delta_press is not None:
            warnings.warn(
                "WARNING: `delta_press` is deprecated, please use `pressure_damping_timescale`."
            )
            pressure_damping_timescale = delta_press
        else:
            pressure_damping_timescale *= time_units

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
                pressure = np.array(pressure)

            if sum(pressure != None) == 0:
                raise ValueError("Pressure cannot be three times None")

            if len(pressure) != 3:
                raise ValueError("Pressure must be a float or a 3d vector")

            if temperature is None or temperature == 0.0:
                raise ValueError("Target temperature for fix nvt/npt/nph cannot be 0")

            pressure[pressure != None] *= pressure_units

            pressure_string = ""
            for coord, value in zip(["x", "y", "z"], pressure):
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

        if initial_temperature > 0:
            self.set_initial_velocity(
                temperature=initial_temperature,
                seed=seed,
                gaussian=True,
                job_name=job_name
            )

    def calc_vcsgc_binary(
        self,
        mc_step_interval=100,
        swap_fraction=0.1,
        temperature_mc=None,
        delta_mu=0.,
        kappa=1000.,
        target_concentration=0.,
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
        Run variance-constrained semi-grand-canonical MD/MC for a binary system.

        https://vcsgc-lammps.materialsmodeling.org

        Note:
            For easy visualization later (with `get_structure`), it is highly recommended that the initial structure
            contain at least one atom of each species.

        Warning:
            Assumes the units are metal, otherwise units for the constraints may be off.

        Args:
            mc_step_interval: How many steps of MD between each set of MC moves. (Default is 100.) Must divide the
                number of ionic steps evenly.
            swap_fraction: The fraction of atoms whose species is swapped at each MC phase. (Default is 0.1.)
            temperature_mc: The temperature for accepting MC steps. (Default is None, which uses the MD temperature.)
            deltamu: The chemical potential difference to replace the first species by the second in eV. (Default is
                0 eV.)
            kappa: Variance constraint for the MC. (Default is 1000.)
            target_concentration: Variance constraint for the MC. (Default is 1000.)
        """
        assert(n_ionic_steps % mc_step_interval == 0)

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
            job_name=job_name
        )

        if temperature_mc is None:
            temperature_mc = temperature
        if seed is None:
            self.generate_seed_from_job(job_name=job_name)

        fix_vcsgc_str = "all sgcmc {0} {1} {2} {3} randseed {4} variance {5} {6}".format(
            str(mc_step_interval),
            str(swap_fraction),
            str(temperature_mc),
            str(delta_mu),
            str(seed),
            str(kappa),
            str(target_concentration)
        )
        self.modify(
            fix___vcsgc=fix_vcsgc_str,
            append_if_not_present=True
        )

    def calc_vcsgc(
        self,
        delta_mu,
        ordered_element_list,
        target_concentration=None,
        kappa=1000.,
        mc_step_interval=100,
        swap_fraction=0.1,
        temperature_mc=None,
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
            delta_mu (dict): A dictionary of N-1 chemical potential differences, where N is the number of species *in
                the potential*. Dictionary keys must be the chemical symbols of the two species the chemical potential
                difference is for, separated by an underscore. E.g. for an X-Y-Z alloy, {'X_Y': -0.4, 'Y_Z': 0.6},
                where -0.4 is the change in free energy to replace an X atom with a Y atom.
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

        if temperature_mc is None:
            temperature_mc = temperature

        if seed is None:
            self.generate_seed_from_job(job_name=job_name)

        if len(delta_mu.keys()) != len(ordered_element_list) - 1:
            raise ValueError("Expected `delta_mu` to have {} items, but got {}".format(
                len(ordered_element_list) - 1,
                len(delta_mu.keys())
            ))

        if set(np.unique(np.array([k.split("_") for k in delta_mu.keys()]).flatten())) != set(ordered_element_list):
            raise ValueError("Chemical potential differences must include all possible species, and may include no "
                             "species not treated by the potential.")

        if not np.all([len(np.unique(k.split("_"))) == 2 for k in delta_mu.keys()]):
            raise ValueError("Chemical potential differences must be between exactly two unique species.")

        # Re-order the given chemical potential differences so they're all relative to the initial species in the potl
        n = len(ordered_element_list)
        mat = np.zeros((n, n))
        vec = np.empty(n)
        order_dict = {}
        for n, el in enumerate(ordered_element_list):
            order_dict[el] = n
        for n, (k, v) in enumerate(delta_mu.items()):
            el1, el2 = k.split('_')
            i, j = order_dict[el1], order_dict[el2]
            mat[n, i] = 1
            mat[n, j] = -1
            vec[n] = v
        mat[-1, 0] = 1
        vec[-1] = 0.  # The initial species is our zero-energy reference
        sol = np.linalg.solve(mat, vec)

        # Apply the actual SGC string
        fix_vcsgc_str = "all sgcmc {0} {1} {2} {3} randseed {4}".format(
            str(mc_step_interval),
            str(swap_fraction),
            str(temperature_mc),
            str(" ".join(str(dmu) for dmu in sol[1:])),
            str(seed),
        )

        # Add VC to the SGC if target concentrations were provided
        if target_concentration is not None:
            if set(target_concentration.keys()) != set(ordered_element_list):
                raise ValueError("Exactly one target concentration must be given for each element treated by the "
                                 "potential.")

            if not np.isclose(np.sum([v for _, v in target_concentration.items()]), 1):
                raise ValueError("Target concentrations must sum to 1.")

            fix_vcsgc_str += " variance {0} {1}".format(
                str(kappa),
                str(" ".join([str(target_concentration[el]) for el in ordered_element_list[1:]]))
            )

        self.modify(
            fix___vcsgc=fix_vcsgc_str,
            append_if_not_present=True
        )
