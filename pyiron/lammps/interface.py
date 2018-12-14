import ast
import h5py
import numpy as np
import os
from pyiron.lammps.structure import LammpsStructure, UnfoldingPrism
from pyiron.base.job.interface import JobInterface
from pyiron.base.pyio.parser import Logstatus, extract_data_from_file
from pyiron.atomistics.md_analysis.trajectory_analysis import unwrap_coordinates


class LammpsInterface(JobInterface):
    def write_input(self, job):
        """
        Call routines that generate the code specific input files

        Returns:

        """
        lmp_structure = self._get_lammps_structure(job=job, structure=job.structure, cutoff_radius=job.cutoff_radius)
        lmp_structure.write_file(file_name="structure.inp", cwd=job.working_directory)
        if job.executable.version and 'dump_modify' in job.input.control._dataset['Parameter'] and \
                (int(job.executable.version.split('.')[0]) < 2016 or
                 (int(job.executable.version.split('.')[0]) == 2016 and
                  int(job.executable.version.split('.')[1]) < 11)):
            job.input.control['dump_modify'] = job.input.control['dump_modify'].replace(' line ', ' ')
        if not all(job.structure.pbc):
            job.input.control['boundary'] = ' '.join(['p' if coord else 'f' for coord in job.structure.pbc])
        job._set_selective_dynamics()
        job.input.control.write_file(file_name="control.inp", cwd=job.working_directory)
        job.input.potential.write_file(file_name="potential.inp", cwd=job.working_directory)
        job.input.potential.copy_pot_files(job.working_directory)

    def collect_output(self, job):
        """

        Returns:

        """
        job.input.from_hdf(job._hdf5)
        if os.path.isfile(job.job_file_name(file_name="dump.h5", cwd=job.working_directory)):
            self._collect_h5md_file(job=job, file_name="dump.h5", cwd=job.working_directory)
        else:
            self._collect_dump_file(job=job, file_name="dump.out", cwd=job.working_directory)
        self._collect_output_log(job=job, file_name="log.lammps", cwd=job.working_directory)
        final_structure = job.get_final_structure()
        with job.project_hdf5.open("output") as hdf_output:
            final_structure.to_hdf(hdf_output)

    @staticmethod
    def _get_lammps_structure(job, structure=None, cutoff_radius=None):
        lmp_structure = LammpsStructure()
        lmp_structure.potential = job.input.potential
        lmp_structure.atom_type = job.input.control["atom_style"]
        if cutoff_radius is not None:
            lmp_structure.cutoff_radius = cutoff_radius
        else:
            lmp_structure.cutoff_radius = job.cutoff_radius
        lmp_structure.el_eam_lst = job.input.potential.get_element_lst()
        if structure is not None:
            lmp_structure.structure = structure
        else:
            lmp_structure.structure = job.structure
        if not set(lmp_structure.structure.get_species_symbols()).issubset(set(lmp_structure.el_eam_lst)):
            raise ValueError('The selected potentials do not support the given combination of elements.')
        return lmp_structure

    @staticmethod
    def _collect_dump_file(job, file_name="dump.out", cwd=None):
        """
        general purpose routine to extract static from a lammps dump file

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = job.job_file_name(file_name=file_name, cwd=cwd)
        tag_dict = {"ITEM: TIMESTEP": {"arg": "0",
                                       "rows": 1,
                                       "h5": "time"},
                    # "ITEM: NUMBER OF ATOMS": {"arg": "0",
                    #                          "rows": 1,
                    #                          "h5": "number_of_atoms"},
                    "ITEM: BOX BOUNDS": {"arg": "0",
                                         "rows": 3,
                                         "h5": "cells",
                                         "func": to_amat},
                    "ITEM: ATOMS": {"arg": ":,:",
                                    "rows": len(job.structure),
                                    "splitArg": True}
                    }

        h5_dict = {"id": "id",
                   "type": "type",
                   "xsu": "coord_xs",
                   "ysu": "coord_ys",
                   "zsu": "coord_zs",
                   "f_ave[1]": "coord_xs",
                   "f_ave[2]": "coord_ys",
                   "f_ave[3]": "coord_zs",
                   "fx": "force_x",
                   "fy": "force_y",
                   "fz": "force_z",
                   }

        lammps_dict = None

        lf = Logstatus()
        lf.extract_file(file_name=file_name,
                        tag_dict=tag_dict,
                        h5_dict=h5_dict,
                        key_dict=lammps_dict)
        lf.combine_xyz('force_x', 'force_y', 'force_z', 'forces')
        lf.combine_xyz('coord_xs', 'coord_ys', 'coord_zs', 'positions')

        prism = UnfoldingPrism(job.structure.cell, digits=15)

        rel_positions = list()

        for ind, (pos, forc, cel) in enumerate(
                zip(lf.status_dict["positions"], lf.status_dict["forces"], lf.status_dict["cells"])):
            cell = cel[1]
            positions = pos[1]
            forces = forc[1]

            # rotation matrix from lammps(unfolded) cell to original cell
            rotation_lammps2orig = np.linalg.inv(prism.R)

            # convert from scaled positions to absolute in lammps cell
            positions = np.array([np.dot(cell.T, r) for r in positions])
            # rotate positions from lammps to original
            positions_atoms = np.array([np.dot(np.array(r), rotation_lammps2orig) for r in positions])

            # rotate forces from lammps to original cell
            forces_atoms = np.array([np.dot(np.array(f), rotation_lammps2orig) for f in forces])

            # unfold cell
            cell = prism.unfold_cell(cell)
            # rotate cell from unfolded lammps to original
            cell_atoms = np.array([np.dot(np.array(f), rotation_lammps2orig) for f in cell])

            lf.status_dict["positions"][ind][1] = positions_atoms

            rel_positions.append(np.dot(positions_atoms, np.linalg.inv(cell_atoms)))

            lf.status_dict["forces"][ind][1] = forces_atoms
            lf.status_dict["cells"][ind][1] = cell_atoms

        del lf.status_dict['id']
        del lf.status_dict['type']
        unwrapped_rel_pos = unwrap_coordinates(positions=np.array(rel_positions), is_relative=True)
        unwrapped_pos = list()
        # print(np.shape(unwrapped_rel_pos))
        for i, cell in enumerate(lf.status_dict["cells"]):
            unwrapped_pos.append(np.dot(np.array(unwrapped_rel_pos[i]), cell[1]))
        lf.status_dict["unwrapped_positions"] = list()
        for pos in unwrapped_pos:
            lf.status_dict["unwrapped_positions"].append([[0], pos])
        with job.project_hdf5.open("output/generic") as hdf_output:
            lf.to_hdf(hdf_output)
        return lf

    # TODO: make rotation of all vectors back to the original as in self.collect_dump_file
    @staticmethod
    def _collect_h5md_file(job, file_name="dump.h5", cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = job.job_file_name(file_name=file_name, cwd=cwd)
        with h5py.File(file_name) as h5md:
            positions = [pos_i.tolist() for pos_i in h5md['/particles/all/position/value']]
            time = [time_i.tolist() for time_i in h5md['/particles/all/position/step']]
            forces = [for_i.tolist() for for_i in h5md['/particles/all/force/value']]
            # following the explanation at: http://nongnu.org/h5md/h5md.html
            cell = [np.eye(3) * np.array(cell_i.tolist()) for cell_i in h5md['/particles/all/box/edges/value']]
        with job.project_hdf5.open("output/generic") as h5_file:
            h5_file['forces'] = np.array(forces)
            h5_file['positions'] = np.array(positions)
            h5_file['time'] = np.array(time)
            h5_file['cells'] = cell

    @staticmethod
    def _collect_errors(job, file_name, cwd=None):
        """

        Args:
            file_name:
            cwd:

        Returns:

        """
        file_name = job.job_file_name(file_name=file_name, cwd=cwd)
        error = extract_data_from_file(file_name, tag="ERROR", num_args=1000)
        if len(error) > 0:
            error = " ".join(error[0])
            raise RuntimeError("Run time error occurred: " + str(error))
        else:
            return True

    def _collect_output_log(self, job, file_name="log.lammps", cwd=None):
        """
        general purpose routine to extract static from a lammps log file

        Args:
            file_name:
            cwd:

        Returns:

        """
        log_file = job.job_file_name(file_name=file_name, cwd=cwd)
        self._collect_errors(job, file_name)
        # log_fileName = "logfile.out"
        # self.collect_errors(log_fileName)
        attr = job.input.control.dataset["Parameter"]
        tag_dict = {"Loop time of": {"arg": "0",
                                     "type": "float",
                                     "h5": "time_loop"},

                    "Memory usage per processor": {"arg": "1",
                                                   "h5": "memory"}
                    }

        if 'minimize' in attr:
            tag_dict["Step Temp PotEng TotEng Pxx Pxy Pxz Pyy Pyz Pzz Volume"] = {"arg": ":,:",
                                                                                  "rows": "Loop",
                                                                                  "splitTag": True}

        elif 'run' in attr:
            num_iterations = ast.literal_eval(extract_data_from_file(log_file, tag="run")[0])
            num_thermo = ast.literal_eval(extract_data_from_file(log_file, tag="thermo")[0])
            num_iterations = num_iterations / num_thermo + 1

            tag_dict["thermo_style custom"] = {"arg": ":",
                                               "type": "str",
                                               "h5": "thermo_style"}
            tag_dict["$thermo_style custom"] = {"arg": ":,:",
                                                "rows": num_iterations,
                                                "splitTag": True}
            # print "tag_dict: ", tag_dict

        h5_dict = {"Step": "steps",
                   "Temp": "temperatures",
                   "PotEng": "energy_pot",
                   "TotEng": "energy_tot",
                   "Pxx": "pressure_x",
                   "Pxy": "pressure_xy",
                   "Pxz": "pressure_xz",
                   "Pyy": "pressure_y",
                   "Pyz": "pressure_yz",
                   "Pzz": "pressure_z",
                   "Volume": "volume",
                   "E_pair": "E_pair",
                   "E_mol": "E_mol"
                   }

        lammps_dict = {"step": "Step",
                       "temp": "Temp",
                       "pe": "PotEng",
                       "etotal": "TotEng",
                       "pxx": "Pxx",
                       "pxy": "Pxy",
                       "pxz": "Pxz",
                       "pyy": "Pyy",
                       "pyz": "Pyz",
                       "pzz": "Pzz",
                       "vol": "Volume"
                       }

        lf = Logstatus()
        lf.extract_file(file_name=log_file,
                        tag_dict=tag_dict,
                        h5_dict=h5_dict,
                        key_dict=lammps_dict)

        lf.store_as_vector = ['energy_tot', 'temperatures', 'steps', 'volume', 'energy_pot']
        # print ("lf_keys: ", lf.status_dict['energy_tot'])

        lf.combine_mat('pressure_x', 'pressure_xy', 'pressure_xz',
                       'pressure_y', 'pressure_yz', 'pressure_z', 'pressures')
        lf.convert_unit('pressures', 0.0001)  # bar -> GPa

        if 'minimize' not in attr:
            if 'thermo_style' in lf.status_dict.keys():
                del lf.status_dict['thermo_style']
        if 'time_loop' in lf.status_dict.keys():
            del lf.status_dict['time_loop']
        if 'memory' in lf.status_dict.keys():
            del lf.status_dict['memory']
        with job.project_hdf5.open("output/generic") as hdf_output:
            lf.to_hdf(hdf_output)


def to_amat(l_list):
    """

    Args:
        l_list:

    Returns:

    """
    lst = np.reshape(l_list, -1)
    if len(lst) == 9:
        xlo_bound, xhi_bound, xy, ylo_bound, yhi_bound, xz, zlo_bound, zhi_bound, yz = lst

    elif len(lst) == 6:
        xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound = lst
        xy, xz, yz = 0., 0., 0.
    else:
        raise ValueError("This format for amat not yet implemented: " + str(len(lst)))

    # > xhi_bound - xlo_bound = xhi -xlo  + MAX(0.0, xy, xz, xy + xz) - MIN(0.0, xy, xz, xy + xz)
    # > xhili = xhi -xlo   = xhi_bound - xlo_bound - MAX(0.0, xy, xz, xy + xz) + MIN(0.0, xy, xz, xy + xz)
    xhilo = (xhi_bound - xlo_bound) - max([0.0, xy, xz, xy + xz]) + min([0.0, xy, xz, xy + xz])

    # > yhilo = yhi -ylo = yhi_bound -ylo_bound - MAX(0.0, yz) + MIN(0.0, yz)
    yhilo = (yhi_bound - ylo_bound) - max([0.0, yz]) + min([0.0, yz])

    # > zhi - zlo = zhi_bound- zlo_bound
    zhilo = (zhi_bound - zlo_bound)

    cell = [[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]]
    return cell
