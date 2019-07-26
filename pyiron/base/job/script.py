# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
import os
import shutil
from pathlib2 import Path
from pyiron.base.job.generic import GenericJob
from pyiron.base.generic.hdfio import FileHDFio
from pyiron.base.generic.parameters import GenericParameters


"""
Jobclass to execute python scripts and jupyter notebooks
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class ScriptJob(GenericJob):
    """
    The ScriptJob class allows to submit Python scripts and Jupyter notebooks to the pyiron job management system.

    Args:
        project (ProjectHDFio): ProjectHDFio instance which points to the HDF5 file the job is stored in
        job_name (str): name of the job, which has to be unique within the project

    Attributes:

        attribute: job_name

            name of the job, which has to be unique within the project

        .. attribute:: status

            execution status of the job, can be one of the following [initialized, appended, created, submitted, running,
                                                                      aborted, collect, suspended, refresh, busy, finished]

        .. attribute:: job_id

            unique id to identify the job in the pyiron database

        .. attribute:: parent_id

            job id of the predecessor job - the job which was executed before the current one in the current job series

        .. attribute:: master_id

            job id of the master job - a meta job which groups a series of jobs, which are executed either in parallel or in
            serial.

        .. attribute:: child_ids

            list of child job ids - only meta jobs have child jobs - jobs which list the meta job as their master

        .. attribute:: project

            Project instance the jobs is located in

        .. attribute:: project_hdf5

            ProjectHDFio instance which points to the HDF5 file the job is stored in

        .. attribute:: job_info_str

            short string to describe the job by it is job_name and job ID - mainly used for logging

        .. attribute:: working_directory

            working directory of the job is executed in - outside the HDF5 file

        .. attribute:: path

            path to the job as a combination of absolute file system path and path within the HDF5 file.

        .. attribute:: version

            Version of the hamiltonian, which is also the version of the executable unless a custom executable is used.

        .. attribute:: executable

            Executable used to run the job - usually the path to an external executable.

        .. attribute:: library_activated

            For job types which offer a Python library pyiron can use the python library instead of an external executable.

        .. attribute:: server

            Server object to handle the execution environment for the job.

        .. attribute:: queue_id

            the ID returned from the queuing system - it is most likely not the same as the job ID.

        .. attribute:: logger

            logger object to monitor the external execution and internal pyiron warnings.

        .. attribute:: restart_file_list

            list of files which are used to restart the calculation from these files.

        .. attribute:: job_type

            Job type object with all the available job types: ['ExampleJob', 'SerialMaster', 'ParallelMaster', 'ScriptJob',
                                                               'ListMaster']

        .. attribute:: script_path

            the absolute path to the python script
    """
    def __init__(self, project, job_name):
        super(ScriptJob, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "Script"
        self._script_path = None
        self.input = GenericParameters(table_name='custom_dict')

    @property
    def script_path(self):
        """
        Python script path

        Returns:
            str: absolute path to the python script
        """
        return self._script_path

    @script_path.setter
    def script_path(self, path):
        """
        Python script path

        Args:
            path (str): relative or absolute path to the python script or a corresponding notebook
        """
        if isinstance(path, str):
            self._script_path = self._get_abs_path(path)
            self.executable = self._executable_command(working_directory=self.working_directory,
                                                       script_path=self._script_path)
        else:
            raise TypeError('path should be a string, but ', path, ' is a ', type(path), ' instead.')

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        self.input.read_only = True

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ScriptJob in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(ScriptJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input['path'] = self._script_path
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ScriptJob from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(ScriptJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            try:
                self.script_path = hdf5_input['path']
                self.input.from_hdf(hdf5_input)
            except TypeError:
                pass

    def write_input(self):
        """
        Copy the script to the working directory - only python scripts and jupyter notebooks are supported
        """
        file_name = os.path.basename(self._script_path)
        shutil.copyfile(src=self._script_path, dst=os.path.join(self.working_directory, file_name))

    def collect_output(self):
        """
        Collect output function updates the master ID entries for all the child jobs created by this script job, if the
        child job is already assigned to an master job nothing happens - master IDs are not overwritten.
        """
        for job in self.project.iter_jobs(recursive=False, convert_to_object=False):
            pr_job = self.project.open(os.path.relpath(job.working_directory, self.project.path))
            for subjob_id in pr_job.get_job_ids(recursive=False):
                if pr_job.db.get_item_by_id(subjob_id)['masterid'] is None:
                    pr_job.db.item_update({'masterid': str(job.job_id)}, subjob_id)

    def run_if_lib(self):
        """
        Compatibility function - but library run mode is not available
        """
        raise NotImplementedError("Library run mode is not implemented for script jobs.")

    def collect_logfiles(self):
        """
        Compatibility function - but no log files are being collected
        """
        pass

    @staticmethod
    def _executable_command(working_directory, script_path):
        """
        internal function to generate the executable command to either use jupyter or python

        Args:
            working_directory (str): working directory of the current job
            script_path (str): path to the script which should be executed in the working directory

        Returns:
            str: executable command
        """
        file_name = os.path.basename(script_path)
        path = os.path.join(working_directory, file_name)
        if file_name[-6:] == '.ipynb':
            return 'jupyter nbconvert --ExecutePreprocessor.timeout=9999999 --to notebook --execute ' + path
        elif file_name[-3:] == '.py':
            return 'python ' + path
        else:
            raise ValueError('Filename not recognized: ', path)

    def _executable_activate_mpi(self):
        """
        Internal helper function to switch the executable to MPI mode
        """
        pass

    @staticmethod
    def _get_abs_path(path):
        """
        internal function to convert absolute or relative paths to absolute paths, using os.path.normpath,
        os.path.abspath and os.path.curdir

        Args:
           path (str): relative or absolute path

        Returns:
            str: absolute path
        """
        return os.path.normpath(os.path.join(os.path.abspath(os.path.curdir), path))


class Notebook(object):
    """
    class for pyiron notebook objects
    """
    @staticmethod
    def get_custom_dict():
        folder = Path('.').cwd().parts[-1]
        hdf_file = Path('.').cwd().parents[1]/folder
        hdf_file = str(hdf_file)+'.h5'
        if Path(hdf_file).exists():
            hdf = FileHDFio(hdf_file)
            custom_dict = GenericParameters()
            for k, v in zip(hdf[folder+'/input/custom_dict/data_dict']['Parameter'],
                            hdf[folder+'/input/custom_dict/data_dict']['Value']):
                custom_dict[k] = v
            return custom_dict
        else:
            print(hdf_file, 'not found')
            return None

    @staticmethod
    def store_custom_output_dict(output_dict):
        folder = Path('.').cwd().parts[-1]
        hdf_file = Path('.').cwd().parents[1] / folder
        hdf_file = str(hdf_file) + '.h5'
        hdf = FileHDFio(hdf_file)
        hdf[folder].create_group('output')
        for k, v in output_dict.items():
            hdf[folder + '/output'][k] = v
