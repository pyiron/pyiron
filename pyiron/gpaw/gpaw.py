import numpy as np
import warnings
from pyiron.gpaw.pyiron_ase import AseJob
from pyiron.atomistics.structure.atoms import pyiron_to_ase, Atoms as PAtoms
from pyiron_base import GenericParameters, Settings
from ase import Atoms

s = Settings()

try:
    from gpaw import GPAW, PW, MethfesselPaxton
except ImportError:
    pass


class Gpaw(AseJob):
    def __init__(self, project, job_name):
        super(Gpaw, self).__init__(project, job_name)
        self.__name__ = "GpawJob"
        self.input = GpawInput()

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, basis):
        if isinstance(basis, PAtoms):
            basis = pyiron_to_ase(basis)
        self._structure = basis

    @property
    def encut(self):
        return self.plane_wave_cutoff

    @encut.setter
    def encut(self, val):
        self.plane_wave_cutoff = val

    @property
    def plane_wave_cutoff(self):
        return self.input["encut"]

    @plane_wave_cutoff.setter
    def plane_wave_cutoff(self, val):
        self.input["encut"] = val

    def get_k_mesh_by_cell(self, kpoints_per_reciprocal_angstrom, cell=None):
        if cell is None:
            if self.structure is None:
                raise AssertionError('structure not set')
            cell = self.structure.cell
        latlens = np.linalg.norm(cell, axis=-1)
        kmesh = np.rint( 2 * np.pi / latlens * kpoints_per_reciprocal_angstrom)
        if kmesh.min() <= 0:
            raise AssertionError("kpoint per angstrom too low")
        return [int(k) for k in kmesh]

    def set_kpoints(
        self,
        mesh=None,
        scheme="MP",
        center_shift=None,
        symmetry_reduction=True,
        manual_kpoints=None,
        weights=None,
        reciprocal=True,
        kpoints_per_reciprocal_angstrom=None,
    ):
        """
        Function to setup the k-points

        Args:
            mesh (list): Size of the mesh (ignored if scheme is not set to 'MP' or kpoints_per_reciprocal_angstrom is set)
            scheme (str): Type of k-point generation scheme (MP/GP(gamma point)/Manual/Line)
            center_shift (list): Shifts the center of the mesh from the gamma point by the given vector in relative
                                 coordinates
            symmetry_reduction (boolean): Tells if the symmetry reduction is to be applied to the k-points
            manual_kpoints (list/numpy.ndarray): Manual list of k-points
            weights(list/numpy.ndarray): Manually supplied weights to each k-point in case of the manual mode
            reciprocal (bool): Tells if the supplied values are in reciprocal (direct) or cartesian coordinates (in
            reciprocal space)
            kpoints_per_reciprocal_angstrom (float): Number of kpoint per angstrom in each direction
        """
        if kpoints_per_reciprocal_angstrom is not None:
            if mesh is not None:
                warnings.warn("mesh value is overwritten by kpoints_per_reciprocal_angstrom")
            mesh = self.get_k_mesh_by_cell(kpoints_per_reciprocal_angstrom=kpoints_per_reciprocal_angstrom)
        if mesh is not None:
            if np.min(mesh) <= 0:
                raise ValueError("mesh values must be larger than 0")
        if center_shift is not None:
            if np.min(center_shift) < 0 or np.max(center_shift) > 1:
                warnings.warn("center_shift is given in relative coordinates")
        self._set_kpoints(
            mesh=mesh,
            scheme=scheme,
            center_shift=center_shift,
            symmetry_reduction=symmetry_reduction,
            manual_kpoints=manual_kpoints,
            weights=weights,
            reciprocal=reciprocal,
        )

    # Backward compatibility
    def get_encut(self):
        return self.encut

    def set_encut(self, encut):
        """
        Sets the plane-wave energy cutoff
        Args:
            encut (float): The energy cutoff in eV
        """
        self.plane_wave_cutoff = encut

    def _set_kpoints(
        self,
        mesh=None,
        scheme="MP",
        center_shift=None,
        symmetry_reduction=True,
        manual_kpoints=None,
        weights=None,
        reciprocal=True,
    ):
        if scheme != "MP":
            raise ValueError("Currently only MP is supported in the pyiron wrapper.")
        if center_shift is not None:
            raise ValueError("centershift is not implemented in the pyiron wrapper.")
        if not symmetry_reduction:
            raise ValueError(
                "symmetry_reduction is not implemented in the pyiron wrapper."
            )
        if manual_kpoints is not None:
            raise ValueError(
                "manual_kpoints are not implemented in the pyiron wrapper."
            )
        if weights is not None:
            raise ValueError("weights are not implemented in the pyiron wrapper.")
        if not reciprocal:
            raise ValueError("reciprocal is not implemented in the pyiron wrapper.")
        self.input["kpoints"] = mesh

    def write_input(self):
        pass

    def collect_output(self):
        pass

    def run_static(self):
        pre_run_mode = self.server.run_mode
        self.server.run_mode.interactive = True
        self.run_if_interactive()
        self.interactive_close()
        self.server.run_mode = pre_run_mode

    def run_if_scheduler(self):
        self._create_working_directory()
        super(GpawJob, self).run_if_scheduler()

    def run_if_interactive(self):
        if self.structure.calc is None:
            kpoints = self.input["kpoints"]
            if isinstance(kpoints, str):
                kpoints = (
                    self.input["kpoints"].replace("[", "").replace("]", "").split()
                )
            self._create_working_directory()
            calc = GPAW(
                mode=PW(float(self.input["encut"])),
                xc=self.input["potential"],
                occupations=MethfesselPaxton(width=float(self.input["sigma"])),
                kpts=kpoints,
                txt=self.working_directory + "/" + self.job_name + ".txt",
            )
            self.structure.set_calculator(calc)
        self.status.running = True
        self.structure.calc.calculate(self.structure)
        self.interactive_collect()

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(GpawJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ExampleJob object in the HDF5 File

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(GpawJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            self.input.from_hdf(hdf5_input)


class GpawInput(GenericParameters):
    """
    Input class for the ExampleJob based on the GenericParameters class.

    Args:
        input_file_name (str): Name of the input file - optional
    """

    def __init__(self, input_file_name=None):
        super(GpawInput, self).__init__(
            input_file_name=input_file_name, table_name="input", comment_char="#"
        )

    def load_default(self):
        """
        Loading the default settings for the input file.
        """
        input_str = """\
kpoints [6,6,6]
sigma 0.1
encut 350
potential PBE
"""
        self.load_string(input_str)
