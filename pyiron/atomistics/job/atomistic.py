# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from ase.io import write as ase_write
import copy

import numpy as np

from pyiron.atomistics.structure.atoms import Atoms
from pyiron.base.objects.generic.parameters import GenericParameters
from pyiron.base.objects.job.generic import GenericJob as GenericJobCore
from pyiron.base.objects.job.master import GenericMaster

try:
    from pyiron.base.core.project.gui import ProjectGUI
except (ImportError, TypeError, AttributeError):
    pass

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class AtomisticGenericJob(GenericJobCore):
    def __init__(self, project, job_name):
        super(AtomisticGenericJob, self).__init__(project, job_name)
        self.__name__ = "AtomisticGenericJob"
        self.__version__ = "0.1"
        self._structure = None
        self._generic_input = GenericInput()
        self.output = GenericOutput(job=self)

    @property
    def structure(self):
        """

        Returns:

        """
        return self._structure

    @structure.setter
    def structure(self, basis):
        """

        Args:
            basis:

        Returns:

        """
        self._generic_input['structure'] = 'atoms'
        self._structure = basis

    def copy_to(self, project=None, new_job_name=None, input_only=False, new_database_entry=True):
        """

        Args:
            destination:
            new_job_name:
            input_only:
            new_database_entry:

        Returns:

        """
        new_generic_job = super(AtomisticGenericJob, self).copy_to(project=project,
                                                                   new_job_name=new_job_name,
                                                                   input_only=input_only,
                                                                   new_database_entry=new_database_entry)
        if not new_generic_job._structure:
            new_generic_job._structure = copy.copy(self._structure)
        return new_generic_job

    def calc_minimize(self, e_tol=1e-8, f_tol=1e-8, max_iter=1000, pressure=None, n_print=1):
        """

        Args:
            e_tol:
            f_tol:
            max_iter:
            pressure:
            n_print:

        Returns:

        """
        self._generic_input['calc_mode'] = 'minimize'
        self._generic_input['max_iter'] = max_iter
        self._generic_input['pressure'] = pressure
        self._generic_input.remove_keys(['temperature', 'n_ionic_steps', 'n_print', 'velocity'])

    def calc_static(self):
        """

        Returns:

        """
        self._generic_input['calc_mode'] = 'static'
        self._generic_input.remove_keys(['max_iter', 'pressure', 'temperature', 'n_ionic_steps', 'n_print', 'velocity'])

    def calc_md(self, temperature=None, pressure=None, n_ionic_steps=1000, time_step=None, n_print=100, delta_temp=1.0,
                delta_press=None, seed=None, tloop=None, rescale_velocity=True, langevin=False):
        self._generic_input['calc_mode'] = 'md'
        self._generic_input['temperature'] = temperature
        self._generic_input['n_ionic_steps'] = n_ionic_steps
        self._generic_input['n_print'] = n_print
        self._generic_input.remove_keys(['max_iter', 'pressure'])

    def from_hdf(self, hdf=None, group_name=None):
        """
        Recreates instance from the hdf5 file
        Args:
            hdf (str): Path to the hdf5 file
            group_name (str): Name of the group which contains the object
        """
        super(AtomisticGenericJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self._hdf5.open("input") as hdf5_input:
            try:
                self._generic_input.from_hdf(hdf5_input)
            except ValueError:
                pass

    def to_hdf(self, hdf=None, group_name=None):
        """
        Stores the instance attributes into the hdf5 file
        """
        super(AtomisticGenericJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self._hdf5.open("input") as hdf5_input:
            self._generic_input.to_hdf(hdf5_input)

    def store_structure(self):
        """

        Returns:

        """
        if self.structure is not None:
            structure_container = self.create_job(self.project.job_type.StructureContainer)
            structure_container.structure = self.structure
            self.parent_id = structure_container.job_id
        else:
            ValueError('There is no structure attached to the current Job.')

    def animate_structure(self, spacefill=True, show_cell=True, stride=1, center_of_mass=False, particle_size=0.5):
        """
        Animates the job if a trajectory is present

        Args:
            spacefill (bool):
            show_cell (bool):
            stride (int): show animation every stride [::stride]
                          use value >1 to make animation faster
                           default=1
            center_of_mass (bool):

        Returns:
            animation: nglview IPython widget

        """
        if not self.status.finished:
            raise ValueError("This job can't be animated until it is finished")
        try:
            import nglview
        except ImportError:
            raise ImportError("The animate() function requires the package nglview to be installed")

        animation = nglview.show_asetraj(self.trajectory(stride=stride, center_of_mass=center_of_mass))
        if spacefill:
            animation.add_spacefill(radius_type='vdw', scale=0.5, radius=particle_size)
            animation.remove_ball_and_stick()
        else:
            animation.add_ball_and_stick()
        if show_cell:
            if self.structure.cell is not None:
                animation.add_unitcell()
        return animation

    def view_structure(self, snapshot=-1, spacefill=True, show_cell=True):
        """

        Args:
            snapshot (int): Snapshot of the trajectory one wants
            spacefill (bool):
            show_cell (bool):

        Returns:
            view: nglview IPython widget

        """
        import nglview
        atoms = self.get_structure(snapshot)
        picture = nglview.show_ase(atoms)
        if spacefill:
            picture.add_spacefill(radius_type='vdw', scale=0.5)
            picture.remove_ball_and_stick()
        else:
            picture.add_ball_and_stick()
        if show_cell:
            if atoms.cell is not None:
                picture.add_unitcell()
        return picture

    def validate_ready_to_run(self):
        """

        Returns:

        """
        if not self.structure and self._generic_input['structure'] == 'atoms':
            raise ValueError('This job does not contain a valid structure: {}'.format(self.job_name))

    def db_entry(self):
        """
        Generate the initial database entry

        Returns:
            (dict): db_dict
        """
        db_dict = super(AtomisticGenericJob, self).db_entry()
        if self.structure:
            parent_structure = self.structure.get_parent_basis()
            db_dict["ChemicalFormula"] = parent_structure.get_chemical_formula()
        return db_dict

    def restart(self, snapshot=-1, job_name=None, job_type=None):
        """
        Restart a new job created from an existing calculation.
        Args:
            project (pyiron.project.Project instance): Project instance at which the new job should be created
            snapshot (int): Snapshot of the calculations which would be the initial structure of the new job
            job_name (str): Job name
            job_type (str): Job type

        Returns:
            new_ham: New job
        """
        new_ham = super(AtomisticGenericJob, self).restart(snapshot=snapshot, job_name=job_name, job_type=job_type)
        if isinstance(new_ham, GenericMaster) and not isinstance(self, GenericMaster):
            new_child = self.restart(snapshot=snapshot, job_name=None, job_type=None)
            new_ham.append(new_child)
        if self.status.finished:
            new_ham.structure = self.get_structure(iteration_step=snapshot)
            new_ham._generic_input['structure'] = 'atoms'
        else:
            new_ham._generic_input['structure'] = 'continue_final'
        return new_ham

    # Required functions
    def continue_with_restart_files(self, job_type=None, job_name=None):
        """

        Args:
            job_type:
            job_name:

        Returns:

        """
        if job_name is None:
            job_name = "{}_continue".format(self.job_name)
        new_ham = self.restart(job_type=job_type, job_name=job_name)
        if self.status.initialized:
            self._job_id = self.save()
        new_ham.parent_id = self.job_id
        new_ham._generic_input['structure'] = 'continue_final'
        return new_ham

    def continue_with_final_structure(self, job_type=None, job_name=None):
        """

        Args:
            job_type:
            job_name:

        Returns:

        """
        if job_name is None:
            job_name = "{}_continue".format(self.job_name)
        if job_type is None:
            job_type = self.__name__
        new_ham = self.create_job(job_type, job_name)
        if self.status.initialized:
            self._job_id = self.save()
        new_ham.parent_id = self.job_id
        if self.status.finished:
            new_ham.structure = self.get_final_structure()
            new_ham._generic_input['structure'] = 'atoms'
        else:
            new_ham._generic_input['structure'] = 'continue_final'
        return new_ham

    def trajectory(self, stride=1, center_of_mass=False):
        """

        Args:
            stride:
            center_of_mass:

        Returns:

        """
        return Trajectory(self['output/generic/positions'][::stride], self.structure.get_parent_basis(),
                          center_of_mass=center_of_mass, cells=self['output/generic/cells'][::stride])

    def write_traj(self, filename, format=None, parallel=True, append=False, stride=1, center_of_mass=False, **kwargs):
        """
        Writes the trajectory in a given file format based on the `ase.io.write`_ function.

        Args:
            filename (str): Filename of the output
            format (str): The specific format of the output
            parallel (bool):
            append (bool):
            stride (int): Writes trajectory every `stride` steps
            center_of_mass (bool): True if the positions are centered on the COM
            **kwargs: Additional ase arguments

        .. _ase.io.write: https://wiki.fysik.dtu.dk/ase/_modules/ase/io/formats.html#write
        """
        traj = self.trajectory(stride=stride, center_of_mass=center_of_mass)
        # Using thr ASE output writer
        ase_write(filename=filename, images=traj, format=format,  parallel=parallel, append=append, **kwargs)

    def _run_if_lib_save(self, job_name=None, structure=None, db_entry=True):
        """

        Args:
            job_name:
            structure:
            db_entry:

        Returns:

        """
        if job_name:
            with self.project_hdf5.open(job_name + '/input') as hdf5_input:
                if structure:
                    structure.to_hdf(hdf5_input)
                else:
                    self.structure.to_hdf(hdf5_input)
        else:
            self.to_hdf()
        return super(AtomisticGenericJob, self)._run_if_lib_save(job_name=job_name, db_entry=db_entry)

    # Compatibility functions
    def get_final_structure(self):
        """

        Returns:

        """
        return self.get_structure()

    def set_kpoints(self, mesh=None, scheme='MP', center_shift=None, symmetry_reduction=True, manual_kpoints=None,
                    weights=None, reciprocal=True):
        raise NotImplementedError("The set_kpoints function is not implemented for this code.")

    def set_encut(self, encut):
        raise NotImplementedError("The set_encut function is not implemented for this code.")

    def get_encut(self):
        raise NotImplementedError("The set_encut function is not implemented for this code.")

    def get_structure(self, iteration_step=-1):
        """
        Gets the structure from a given iteration step of the simulation (MD/ionic relaxation). For static calculations
        there is only one ionic iteration step
        Args:
            iteration_step (int): Step for which the structure is requested

        Returns:
            atomistics.structure.atoms.Atoms object


        """
        if not (self.structure is not None):
            raise AssertionError()
        snapshot = self.structure.copy()
        snapshot.cell = self.get("output/generic/cells")[iteration_step]
        snapshot.positions = self.get("output/generic/positions")[iteration_step]
        return snapshot

    def gui(self):
        """

        Returns:

        """
        ProjectGUI(self)

    def _structure_to_hdf(self):
        with self.project_hdf5.open("input") as hdf5_input:
            if self._generic_input['structure'] == 'atoms':
                self.structure.to_hdf(hdf5_input)

    def _structure_from_hdf(self):
        with self.project_hdf5.open("input") as hdf5_input:
            if self._generic_input['structure'] == 'atoms':
                self.structure = Atoms().from_hdf(hdf5_input)

    def _write_chemical_formular_to_database(self):
        if self.structure:
            parent_structure = self.structure.get_parent_basis()
            self.project.db.item_update({"ChemicalFormula": parent_structure.get_chemical_formula()}, self._job_id)

    def _before_successor_calc(self, ham):
        if ham._generic_input['structure'] == 'continue_final':
            ham.structure = self.get_final_structure()
            ham.to_hdf()


class Trajectory(object):
    """
    A trajectory instance compatible with the ase.io class

    Args:
        positions (numpy.ndarray): The array of the trajectory in cartesian coordinates
        structure (pyiron.atomistics.structure.atoms.Atoms): The initial structure instance from which the species info
                                                             is derived
        center_of_mass (bool): False (default) if the specified positions are w.r.t. the origin
        cells (numpy.ndarray): Optional argument of the cell shape at every time step (Nx3x3 array) when the volume
                                varies
    """

    def __init__(self, positions, structure, center_of_mass=False, cells=None):
        if center_of_mass:
            pos = np.copy(positions)
            pos[:, :, 0] = (pos[:, :, 0].T - np.mean(pos[:, :, 0], axis=1)).T
            pos[:, :, 1] = (pos[:, :, 1].T - np.mean(pos[:, :, 1], axis=1)).T
            pos[:, :, 2] = (pos[:, :, 2].T - np.mean(pos[:, :, 2], axis=1)).T
            self._positions = pos
        else:
            self._positions = positions
        self._structure = structure
        self._cells = cells

    def __getitem__(self, item):
        new_structure = self._structure.copy()
        if self._cells is not None:
            new_structure.cell = self._cells[item]
        new_structure.positions = self._positions[item]
        # This step is necessary for using ase.io.write for trajectories
        new_structure.arrays['positions'] = new_structure.positions
        new_structure.arrays['cells'] = new_structure.cell
        return new_structure

    def __len__(self):
        return len(self._positions)


class GenericInput(GenericParameters):
    def __init__(self, input_file_name=None, table_name="generic"):
        super(GenericInput, self).__init__(input_file_name=input_file_name, table_name=table_name, comment_char="#",
                                           separator_char="=")

    def load_default(self):
        """
        Loads the default file content
        """
        file_content = '''\
calc_mode=static # static, minimize, md
structure=atoms # atoms, continue_final
'''
        self.load_string(file_content)


class GenericOutput(object):
    def __init__(self, job):
        self._job = job

    @property
    def cells(self):
        return self._job['output/generic/cells']

    @property
    def energy_pot(self):
        return self._job['output/generic/energy_pot']

    @property
    def energy_tot(self):
        return self._job['output/generic/energy_tot']

    @property
    def forces(self):
        return self._job['output/generic/forces']

    @property
    def positions(self):
        return self._job['output/generic/positions']

    @property
    def pressures(self):
        return self._job['output/generic/pressures']

    @property
    def steps(self):
        return self._job['output/generic/steps']

    @property
    def temperature(self):
        return self._job['output/generic/temperature']

    @property
    def computation_time(self):
        return self._job['output/generic/computation_time']

    @property
    def unwrapped_positions(self):
        return self._job['output/generic/unwrapped_positions']

    @property
    def volume(self):
        return self._job['output/generic/volume']

    def __dir__(self):
        hdf5_path = self._job['output/generic']
        if hdf5_path is not None:
            return hdf5_path.list_nodes()
        else:
            return []
