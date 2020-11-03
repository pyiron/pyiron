# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import subprocess
import os
import time
import posixpath
import warnings
from pyiron_base import Settings, GenericParameters, Executable
from pyiron.atomistics.job.interactivewrapper import (
    InteractiveWrapper,
    ReferenceJobOutput,
)
from pyiron.atomistics.job.interactive import InteractiveInterface
from pyiron.sphinx.base import InputWriter

__author__ = "Jan Janssen, Osamu Waseda"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


s = Settings()


class SxExtOpt(InteractiveInterface):
    def __init__(
        self,
        structure,
        working_directory=None,
        maxDist=5,
        ionic_steps=1000,
        ionic_energy=None,
        ionic_forces=None,
        ionic_energy_tolerance=1.0e-3,
        ionic_force_tolerance=1.0e-2,
        max_step_length=1.0e-1,
        soft_mode_damping=1,
        executable=None,
        ssa=False,
    ):
        if ionic_forces is not None:
            warnings.warn(
                (
                        'ionic_forces is deprecated as of vers. 0.3.0.' +
                        'It is not guaranteed to be in service in vers. 0.4.0.' +
                        'Use ionic_force_tolerance instead'
                ), DeprecationWarning
            )
            ionic_force_tolerance = ionic_forces
        if ionic_energy is not None:
            warnings.warn(
                (
                        'ionic_energy is deprecated as of vers. 0.3.0.' +
                        'It is not guaranteed to be in service in vers. 0.4.0.' +
                        'Use ionic_energy_tolerance instead'
                ), DeprecationWarning
            )
            ionic_energy_tolerance = ionic_energy
        super().__init__()
        self.__name__ = "SxExtOpt"
        if working_directory is None:
            warnings.warn("WARNING: working_directory not set; current folder is used")
            working_directory = os.getcwd()
        self._interactive_library = None
        self._interactive_library_read = None
        self.working_directory = working_directory
        if executable is None:
            executable = Executable(
                path_binary_codes=s.resource_paths,
                codename="SxExtOptInteractive",
                module=self.__module__.split(".")[1],
                overwrite_nt_flag=False,
            ).executable_path
        self._start_process(
            structure=structure,
            executable=executable,
            maxDist=maxDist,
            ionic_steps=ionic_steps,
            ionic_energy_tolerance=ionic_energy_tolerance,
            ionic_force_tolerance=ionic_force_tolerance,
            max_step_length=max_step_length,
            soft_mode_damping=soft_mode_damping,
            selective_dynamics="selective_dynamics" in structure._tag_list.keys(),
            ssa=ssa,
        )
        self._cell = structure.cell
        if ssa:
            self._elements = structure.get_parent_symbols()
        else:
            magmom = structure.get_initial_magnetic_moments()
            magmom[magmom!=None] = np.round(magmom[magmom!=None], decimals=1)
            magmom = np.char.mod('%s', magmom)
            self._elements = np.char.add(structure.get_parent_symbols(), magmom)
            self._elements = np.char.replace(self._elements, '-', 'm')
            self._elements = np.char.replace(self._elements, '.', 'p')
        self._positions = structure.positions
        self._converged = False

    def _start_process(
        self,
        structure,
        executable,
        maxDist=5,
        ionic_steps=1000,
        ionic_energy_tolerance=1.0e-3,
        ionic_force_tolerance=1.0e-2,
        max_step_length=1.0e-1,
        soft_mode_damping=1,
        selective_dynamics=False,
        ssa=False,
    ):
        if selective_dynamics:
            input_writer_obj = InputWriter()
            input_writer_obj.structure = structure
            if ssa:
                input_writer_obj.structure.set_initial_magnetic_moments(len(structure)*[None])
            input_writer_obj.write_structure(
                file_name="structure.sx",
                cwd=self.working_directory,
                structure_str=None,
                symmetry_enabled=True,
                keep_angstrom=True,
            )
        self._write_input(
            working_directory=self.working_directory,
            maxDist=maxDist,
            ionic_steps=ionic_steps,
            ionic_energy_tolerance=ionic_energy_tolerance,
            ionic_force_tolerance=ionic_force_tolerance,
            max_step_length=max_step_length,
            soft_mode_damping=soft_mode_damping,
            selective_dynamics=selective_dynamics,
        )

        shell = os.name == "nt"
        try:
            with open(
                posixpath.join(self.working_directory, "out.txt"), mode="w"
            ) as f_out:
                with open(
                    posixpath.join(self.working_directory, "error.txt"), mode="w"
                ) as f_err:
                    self._process = subprocess.Popen(
                        [executable],
                        cwd=self.working_directory,
                        shell=shell,
                        stdout=f_out,
                        stderr=f_err,
                        universal_newlines=True,
                    )
        except subprocess.CalledProcessError as e:
            raise ValueError("run_job.py crashed")
        while not self._interactive_pipes_initialized(self.working_directory):
            time.sleep(1)
        self._interactive_initialize_interface()

    @staticmethod
    def _write_input(
        working_directory,
        maxDist=5,
        ionic_steps=1000,
        ionic_energy_tolerance=1.0e-3,
        ionic_force_tolerance=1.0e-2,
        max_step_length=1.0e-1,
        soft_mode_damping=1,
        selective_dynamics=False,
    ):
        with open(os.path.join(working_directory, "optim.sx"), "w") as f:
            content = (
                "main { ricQN { ric { maxDist = %f; withAngles; } maxSteps = %i; dEnergy = %f; dF = %f; maxStepLength = %f; softModeDamping = %f;}}"
                % (
                    maxDist,
                    ionic_steps,
                    ionic_energy_tolerance,
                    ionic_force_tolerance,
                    max_step_length,
                    soft_mode_damping,
                )
            )
            if selective_dynamics:
                content += "structure { include <structure.sx>; }"
            f.write(content)

    @staticmethod
    def _interactive_pipes_initialized(working_directory):
        return os.path.exists(
            os.path.join(working_directory, "control")
        ) and os.path.exists(os.path.join(working_directory, "response"))

    def _interactive_write_line(self, line):
        self._interactive_library.write("%s\n" % line)
        self._interactive_library.flush()

    def _interactive_initialize_interface(self):
        self._interactive_library_read = open(
            os.path.join(self.working_directory, "control"), "r"
        )
        self._interactive_library = open(
            os.path.join(self.working_directory, "response"), "w"
        )

    def interactive_close(self):
        if self.interactive_is_activated():
            self._interactive_library.close()
            self._interactive_library_read.close()
            self._delete_named_pipes(working_directory=self.working_directory)

    @staticmethod 
    def _delete_named_pipes(working_directory):
        for file in ["control", "response"]:
            file_path = posixpath.join(working_directory, file)
            if os.path.exists(file_path):
                os.remove(file_path)

    def interactive_is_activated(self):
        if self._interactive_library is None:
            return False
        else:
            return True

    def _write_cell(self, cell):
        for c in cell:
            self._interactive_write_line("%.16f %.16f %.16f" % (c[0], c[1], c[2]))

    def _write_number_of_atoms(self, count):
        self._interactive_write_line("%s" % count)

    def _write_positions(self, positions, elements):
        for pos, el in zip(positions, elements):
            self._interactive_write_line(
                "%.16f %.16f %.16f %s" % (pos[0], pos[1], pos[2], str(el))
            )

    def _write_energy(self):
        self._interactive_write_line("0")

    def _write_forces(self, forces):
        for f in forces:
            self._interactive_write_line("%.16f %.16f %.16f" % (f[0], f[1], f[2]))

    def _read_positions(self, count):
        return [
            [float(c) for c in self._interactive_library_read.readline().split()]
            for _ in range(count)
        ]

    def set_forces(self, forces):
        line = self._interactive_library_read.readline().split()
        if len(line) == 0 or line[0] == "end":
            print("Ending calculation")
            self._converged = True
        elif line[0] == "get":
            if line[1] == "cell":
                self._write_cell(cell=self._cell)
            elif line[1] == "natoms":
                self._write_number_of_atoms(count=len(self._positions))
            elif line[1] == "positions":
                self._write_positions(
                    positions=self._positions, elements=self._elements
                )
            elif line[1] == "energy":
                self._write_energy()
            elif line[1] == "forces":
                self._write_forces(forces=forces)
            else:
                raise ValueError("Unknown command:", line)
            self.set_forces(forces)
        elif line[0] == "run":
            self.set_forces(forces)
        elif line[0] == "set":
            if line[1] == "positions":
                self._positions = np.array(self._read_positions(count=len(forces)))
            else:
                raise ValueError("Unknown command:", line)
        else:
            raise ValueError("Unknown command:", line)

    def get_positions(self):
        return self._positions

    def end(self):
        while not self.converged:
            self.set_forces(np.zeros_like(self._positions))

    @property
    def converged(self):
        return self._converged

    def __del__(self):
        self.end()
        if self.interactive_is_activated():
            self.interactive_close()
        else:
            self._delete_named_pipes(working_directory=self.working_directory)


class SxExtOptInteractive(InteractiveWrapper):
    def __init__(self, project, job_name):
        super(SxExtOptInteractive, self).__init__(project, job_name)
        self.__name__ = "SxExtOptInteractive"
        self.__version__ = (
            None
        )  # Reset the version number to the executable is set automatically
        self._executable_activate()
        self.input = Input()
        self.output = SxExtOptOutput(job=self)
        self._interactive_interface = None
        self._interactive_number_of_steps = 0

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        super(SxExtOptInteractive, self).set_input_to_read_only()
        self.input.read_only = True

    def write_input(self):
        pass

    def run_static(self):
        """
        The run if modal function is called by run to execute the simulation, while waiting for the output. For this we
        use subprocess.check_output()
        """
        self._create_working_directory()
        self._interactive_interface = SxExtOpt(
            structure=self.ref_job.structure,
            working_directory=self.working_directory,
            maxDist=int(self.input["maxDist"]),
            ionic_steps=int(self.input["ionic_steps"]),
            ionic_energy_tolerance=float(self.input["ionic_energy_tolerance"]),
            ionic_force_tolerance=float(self.input["ionic_force_tolerance"]),
            max_step_length=float(self.input["max_step_length"]),
            soft_mode_damping=float(self.input["soft_mode_damping"]),
            executable=self.executable.executable_path,
            ssa=self.input['ssa'],
        )
        self.status.running = True
        self._logger.info("job status: %s", self.status)
        new_positions = self.ref_job.structure.positions
        self.ref_job_initialize()
        while (
            self._interactive_number_of_steps < self.input["ionic_steps"]
            and not self._interactive_interface.converged
        ):
            str_temp = self.ref_job.structure
            str_temp.positions = new_positions
            self.ref_job.structure = str_temp
            if self.ref_job.server.run_mode.interactive:
                self._logger.debug("SxExtOpt: step start!")
                self.ref_job.run()
                self._logger.debug("SxExtOpt: step finished!")
            else:
                self.ref_job.run(delete_existing_job=True)
            self._interactive_interface.set_forces(forces=self.get_forces())
            new_positions = self._interactive_interface.get_positions()
            self._interactive_number_of_steps += 1
        self.status.collect = True
        if self.ref_job.server.run_mode.interactive:
            self.ref_job.interactive_close()
        self._interactive_interface.interactive_close()
        self.run()

    def get_forces(self):
        ff = np.array(self.ref_job.output.forces[-1])
        if hasattr(self.ref_job.structure, "selective_dynamics"):
            ff[np.array(self.ref_job.structure.selective_dynamics) == False] = 0
            return ff
        return ff

    def convergence_check(self):
        """
        Validate the convergence of the calculation.

        Returns:
             (bool): If the calculation is converged
        """
        if self._interactive_number_of_steps < self.input["ionic_steps"]:
            return True
        else:
            return False


class Input(GenericParameters):
    """
    class to control the generic input for a Sphinx calculation.

    Args:
        input_file_name (str): name of the input file
        table_name (str): name of the GenericParameters table
    """

    def __init__(self, input_file_name=None, table_name="input"):
        super(Input, self).__init__(
            input_file_name=input_file_name,
            table_name=table_name,
            comment_char="//",
            separator_char="=",
            end_value_char=";",
        )

    def load_default(self):
        """
        Loads the default file content
        """
        file_content = (
            "ionic_steps = 1000 // maximum number of ionic steps\n"
            "ionic_energy_tolerance = 1.0e-3\n"
            "ionic_force_tolerance = 1.0e-2\n"
            "maxDist = 5 // maximum possible distance for considering neighbors\n"
            "max_step_length = 1.0e-1 // maximum displacement at each step\n"
            "ssa = False // ignore different magnetic moment values when internal symmetries are considered\n"
            "soft_mode_damping = 1.0 // Tikhonov damper\n"
        )
        self.load_string(file_content)


class SxExtOptOutput(ReferenceJobOutput):
    def __init__(self, job):
        super(SxExtOptOutput, self).__init__(job=job)
