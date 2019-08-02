# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
# import copy
import signal
from datetime import datetime
import os
# import sys
import posixpath
import psutil
import multiprocessing
from pyiron.base.job.wrapper import JobWrapper
from pyiron.base.settings.generic import Settings
from pyiron.base.job.executable import Executable
from pyiron.base.job.jobstatus import JobStatus
from pyiron.base.job.core import JobCore
from pyiron.base.generic.util import static_isinstance
from pyiron.base.server.generic import Server
import subprocess
import shutil
import warnings

"""
Generic Job class extends the JobCore class with all the functionality to run the job object.
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()

intercepted_signals = [signal.SIGINT, signal.SIGTERM, signal.SIGABRT]  # , signal.SIGQUIT]


class GenericJob(JobCore):
    """
    Generic Job class extends the JobCore class with all the functionality to run the job object. From this class
    all specific Hamiltonians are derived. Therefore it should contain the properties/routines common to all jobs.
    The functions in this module should be as generic as possible.
     
    Args:
        project (ProjectHDFio): ProjectHDFio instance which points to the HDF5 file the job is stored in
        job_name (str): name of the job, which has to be unique within the project

    Attributes:

        .. attribute:: job_name

            name of the job, which has to be unique within the project

        .. attribute:: status

            execution status of the job, can be one of the following [initialized, appended, created, submitted,
                                                                      running, aborted, collect, suspended, refresh,
                                                                      busy, finished]

        .. attribute:: job_id

            unique id to identify the job in the pyiron database

        .. attribute:: parent_id

            job id of the predecessor job - the job which was executed before the current one in the current job series

        .. attribute:: master_id

            job id of the master job - a meta job which groups a series of jobs, which are executed either in parallel
            or in serial.

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

            For job types which offer a Python library pyiron can use the python library instead of an external
            executable.

        .. attribute:: server

            Server object to handle the execution environment for the job.

        .. attribute:: queue_id

            the ID returned from the queuing system - it is most likely not the same as the job ID.

        .. attribute:: logger

            logger object to monitor the external execution and internal pyiron warnings.

        .. attribute:: restart_file_list

            list of files which are used to restart the calculation from these files.

        .. attribute:: exclude_nodes_hdf

            list of nodes which are excluded from storing in the hdf5 file.

        .. attribute:: exclude_groups_hdf

            list of groups which are excluded from storing in the hdf5 file.

        .. attribute:: job_type

            Job type object with all the available job types: ['ExampleJob', 'SerialMaster', 'ParallelMaster',
                                                               'ScriptJob', 'ListMaster']
    """
    def __init__(self, project, job_name):
        super(GenericJob, self).__init__(project, job_name)
        self.__name__ = "GenericJob"
        self.__version__ = "0.4"
        self._server = Server()
        self._logger = s.logger
        self._executable = None
        self._import_directory = None
        self._status = JobStatus(db=project.db, job_id=self.job_id)
        self.refresh_job_status()
        self._restart_file_list = list()
        self._restart_file_dict = dict()
        self._exclude_nodes_hdf = list()
        self._exclude_groups_hdf = list()
        self._process = None
        self._compress_by_default = False
        self.interactive_cache = None

        for sig in intercepted_signals:
            signal.signal(sig,  self.signal_intercept)

    def signal_intercept(self,sig,frame):
        try:
            self._logger.info('Job {} intercept signal {}, job is shutting down'.format(self._job_id, sig))
            self.drop_status_to_aborted()
        except:
            raise
        # finally:
        #     if sig in intercepted_signals:
        #         sys.exit(0)

    def drop_status_to_aborted(self):
        self.refresh_job_status()
        if not (self.status.finished or self.status.suspended):
            self.status.aborted = True

    @property
    def version(self):
        """
        Get the version of the hamiltonian, which is also the version of the executable unless a custom executable is
        used.
        
        Returns:
            str: version number
        """
        if self.__version__:
            return self.__version__
        else:
            self._executable_activate()
            if self._executable is not None:
                return self._executable.version
            else:
                return None

    @version.setter
    def version(self, new_version):
        """
        Set the version of the hamiltonian, which is also the version of the executable unless a custom executable is
        used.
        
        Args:
            new_version (str): version
        """
        self._executable_activate()
        self._executable.version = new_version

    @property
    def executable(self):
        """
        Get the executable used to run the job - usually the path to an external executable.
        
        Returns:
            (str/pyiron.base.job.executable.Executable): exectuable path
        """
        self._executable_activate()
        return self._executable

    @executable.setter
    def executable(self, exe):
        """
        Set the executable used to run the job - usually the path to an external executable.
        
        Args:
            exe (str): executable path, if no valid path is provided an executable is chosen based on version.
        """
        self._executable_activate()
        self._executable.executable_path = exe

    @property
    def server(self):
        """
        Get the server object to handle the execution environment for the job.

        Returns:
            Server: server object
        """
        return self._server

    @server.setter
    def server(self, server):
        """
        Set the server object to handle the execution environment for the job.
        Args:
            server (Server): server object
        """
        self._server = server

    @property
    def queue_id(self):
        """
        Get the queue ID, the ID returned from the queuing system - it is most likely not the same as the job ID.

        Returns:
            int: queue ID
        """
        return self.server.queue_id

    @queue_id.setter
    def queue_id(self, qid):
        """
        Set the queue ID, the ID returned from the queuing system - it is most likely not the same as the job ID.

        Args:
            qid (int): queue ID
        """
        self.server.queue_id = qid

    @property
    def logger(self):
        """
        Get the logger object to monitor the external execution and internal pyiron warnings.

        Returns:
            logging.getLogger(): logger object
        """
        return self._logger

    @property
    def restart_file_list(self):
        """
        Get the list of files which are used to restart the calculation from these files.

        Returns:
            list: list of files
        """
        return self._restart_file_list

    @restart_file_list.setter
    def restart_file_list(self, filenames):
        """
        Append new files to the restart file list - the list of files which are used to restart the calculation from.
        
        Args:
            filenames (list):
        """
        for f in filenames:
            if not (os.path.isfile(f)):
                raise IOError("File: {} does not exist".format(f))
            self.restart_file_list.append(f)

    @property
    def restart_file_dict(self):
        """
        A dictionary of the new name of the copied restart files
        """
        for actual_name in [os.path.basename(f) for f in self._restart_file_list]:
            if actual_name not in self._restart_file_dict.keys():
                self._restart_file_dict[actual_name] = actual_name
        return self._restart_file_dict

    @restart_file_dict.setter
    def restart_file_dict(self, val):
        if not isinstance(val, dict):
            raise ValueError("restart_file_dict should be a dictionary!")
        else:
            self._restart_file_dict = val

    @property
    def exclude_nodes_hdf(self):
        """
        Get the list of nodes which are excluded from storing in the hdf5 file

        Returns:
            nodes(list)
        """
        return self._exclude_nodes_hdf

    @exclude_nodes_hdf.setter
    def exclude_nodes_hdf(self, val):
        if isinstance(val, str):
            val = [val]
        elif not hasattr(val, '__len__'):
            raise ValueError('Wrong type of variable.')
        self._exclude_nodes_hdf = val

    @property
    def exclude_groups_hdf(self):
        """
        Get the list of groups which are excluded from storing in the hdf5 file

        Returns:
            groups(list)
        """
        return self._exclude_groups_hdf

    @exclude_groups_hdf.setter
    def exclude_groups_hdf(self, val):
        if isinstance(val, str):
            val = [val]
        elif not hasattr(val, '__len__'):
            raise ValueError('Wrong type of variable.')
        self._exclude_groups_hdf = val

    @property
    def job_type(self):
        """
        Job type object with all the available job types: ['ExampleJob', 'SerialMaster', 'ParallelMaster', 'ScriptJob',
                                                           'ListMaster']
        Returns:
            JobTypeChoice: Job type object
        """
        return self.project.job_type

    @property
    def working_directory(self):
        """
        Get the working directory of the job is executed in - outside the HDF5 file. The working directory equals the
        path but it is represented by the filesystem:
            /absolute/path/to/the/file.h5/path/inside/the/hdf5/file
        becomes:
            /absolute/path/to/the/file_hdf5/path/inside/the/hdf5/file

        Returns:
            str: absolute path to the working directory
        """
        if self._import_directory:
            return self._import_directory
        elif not self.project_hdf5.working_directory:
            self._create_working_directory()
        return self.project_hdf5.working_directory

    def collect_logfiles(self):
        """
        Collect the log files of the external executable and store the information in the HDF5 file. This method has
        to be implemented in the individual hamiltonians.
        """
        pass

    def write_input(self):
        """
        Write the input files for the external executable. This method has to be implemented in the individual
        hamiltonians.
        """
        raise NotImplementedError("write procedure must be defined for derived Hamilton!")

    def collect_output(self):
        """
        Collect the output files of the external executable and store the information in the HDF5 file. This method has
        to be implemented in the individual hamiltonians.
        """
        raise NotImplementedError("read procedure must be defined for derived Hamilton!")

    def append(self, job):
        """
        Metajobs like GenericMaster, ParallelMaster, SerialMaser or ListMaster allow other jobs to be appended. In the
        GenericJob definition this is only a template function.
        """
        raise NotImplementedError("append procedure must be defined for derived Hamilton!")

    def suspend(self):
        """
        Suspend the job by storing the object and its state persistently in HDF5 file and exit it.
        """
        self.to_hdf()
        self.status.suspended = True
        self._logger.info('{}, status: {}, job has been suspended'.format(self.job_info_str, self.status))
        self.clear_job()

    def refresh_job_status(self):
        """
        Refresh job status by updating the job status with the status from the database if a job ID is available.
        """
        if self.job_id:
            self._status = JobStatus(initial_status=self.project.db.get_item_by_id(self.job_id)["status"],
                                     db=self.project.db, job_id=self.job_id)

    def clear_job(self):
        """
        Convenience function to clear job info after suspend. Mimics deletion of all the job info after suspend in a
        local test environment.
        """
        del self.__name__
        del self.__version__
        del self._executable
        del self._name
        del self._server
        del self._logger
        del self._parent_id
        del self._master_id
        del self._import_directory
        del self._status
        del self._restart_file_list
        del self._restart_file_dict
        # del self._process
        # del self._hdf5
        del self._job_id
        del self._status

    def copy(self):
        """
        Copy the GenericJob object which links to the job and its HDF5 file

        Returns:
            GenericJob: New GenericJob object pointing to the same job
        """
        if not self.project_hdf5.file_exists:
            delete_file_after_copy = True
        else:
            delete_file_after_copy = False
        self.to_hdf()
        self_class = self.__class__
        copied_self = self_class(job_name=self.job_name, project=self.project_hdf5.open('..'))
        copied_self.from_hdf()
        if delete_file_after_copy:
            self.project_hdf5.remove_file()
        copied_self._job_id = None
        return copied_self

    def copy_to(self, project=None, new_job_name=None, input_only=False, new_database_entry=True):
        """
        Copy the content of the job including the HDF5 file to a new location

        Args:
            project (ProjectHDFio): project to copy the job to
            new_job_name (str): to duplicate the job within the same porject it is necessary to modify the job name
                                - optional
            input_only (bool): [True/False] to copy only the input - default False
            new_database_entry (bool): [True/False] to create a new database entry - default True

        Returns:
            GenericJob: GenericJob object pointing to the new location.
        """
        if project is None and new_job_name is None:
            raise ValueError('copy_to requires either a new project or a new_job_name.')
        if not self.project_hdf5.file_exists:
            self.to_hdf()
            delete_file_after_copy = True
        else:
            delete_file_after_copy = False
        if project is None and new_job_name is not None:
            new_generic_job = self.copy()
            new_generic_job.reset_job_id()
            if len(self.project_hdf5.h5_path.split('/')) > 2:
                new_location = self.project_hdf5.open('../' + new_job_name)
            else:
                new_location = self.project_hdf5.__class__(self.project, new_job_name, h5_path='/' + new_job_name)
            new_generic_job._name = new_job_name
            new_generic_job.project_hdf5.copy_to(new_location, maintain_name=False)
            new_generic_job.project_hdf5 = new_location
            if new_database_entry:
                new_generic_job.save()
        else:
            new_generic_job = super(GenericJob, self).copy_to(project, new_database_entry=new_database_entry)
            new_generic_job.reset_job_id(job_id=new_generic_job.job_id)
            new_generic_job.from_hdf()
        if input_only:
            if 'output' in new_generic_job.project_hdf5.list_groups():
                del new_generic_job.project_hdf5[posixpath.join(new_generic_job.project_hdf5.h5_path, 'output')]
        if delete_file_after_copy:
            self.project_hdf5.remove_file()
        if project is not None and new_job_name:
            new_generic_job.job_name = new_job_name
        return new_generic_job

    def copy_file_to_working_directory(self, file):
        """
        Copy a specific file to the working directory before the job is executed.

        Args:
            file (str): path of the file to be copied.
        """
        if os.path.isabs(file):
            self.restart_file_list.append(file)
        else: 
            self.restart_file_list.append(file)
    
    def copy_template(self, project, new_job_name=None):
        """
        Copy the content of the job including the HDF5 file but without the output data to a new location

        Args:
            project (ProjectHDFio): project to copy the job to
            new_job_name (str): to duplicate the job within the same porject it is necessary to modify the job name
                                - optional

        Returns:
            GenericJob: GenericJob object pointing to the new location.
        """
        return self.copy_to(project=project, new_job_name=new_job_name, input_only=True, new_database_entry=False)

    def _kill_child(self):
        """
        Internal helper function to kill a child process.
        """
        if not self.server.run_mode.queue and (self.status.running or self.status.submitted):
            for proc in psutil.process_iter():
                try:
                    pinfo = proc.as_dict(attrs=['pid', 'cwd'])
                except psutil.NoSuchProcess:
                    pass
                else:
                    if pinfo['cwd'] is not None and pinfo['cwd'].startswith(self.working_directory):
                        job_process = psutil.Process(pinfo['pid'])
                        job_process.kill()
    
    def remove_child(self):
        """
        internal function to remove command that removes also child jobs.
        Do never use this command, since it will destroy the integrity of your project.
        """
        self._kill_child()
        super(GenericJob, self).remove_child()
    
    def kill(self):
        if self.status.running or self.status.submitted:
            master_id, parent_id = self.master_id, self.parent_id
            self.remove()
            self.reset_job_id()
            self.master_id, self.parent_id = master_id, parent_id
        else:
            raise ValueError('The kill() function is only available during the execution of the job.')
    
    def validate_ready_to_run(self):
        """
        Validate that the calculation is ready to be executed. By default no generic checks are performed, but one could
        check that the input information is complete or validate the consistency of the input at this point.
        """
        pass

    def check_setup(self):
        """
        Checks whether certain parameters (such as plane wave cutoff radius in DFT) are changed from the pyiron standard
        values to allow for a physically meaningful results. This function is called manually or only when the job is
        submitted to the queueing system.
        """
        pass

    def reset_job_id(self, job_id=None):
        """
        Reset the job id sets the job_id to None in the GenericJob as well as all connected modules like JobStatus.
        """
        if job_id is not None:
            job_id = int(job_id)
        self._job_id = job_id
        self._status = JobStatus(db=self.project.db, job_id=self._job_id)

    def run(self, run_again=False, repair=False, debug=False, run_mode=None):
        """
        This is the main run function, depending on the job status ['initialized', 'created', 'submitted', 'running',
        'collect','finished', 'refresh', 'suspended'] the corresponding run mode is chosen.
        
        Args:
            run_again (bool): Delete the existing job and run the simulation again.
            repair (bool): Set the job status to created and run the simulation again.
            debug (bool): Debug Mode - defines the log level of the subprocess the job is executed in.
            run_mode (str): ['modal', 'non_modal', 'queue', 'manual'] overwrites self.server.run_mode
        """
        try:
            self._logger.info('run {}, status: {}'.format(self.job_info_str, self.status))
            status = self.status.string
            if run_mode is not None:
                self.server.run_mode = run_mode
            if run_again and self.job_id:
                self._logger.info("run repair "+str(self.job_id))
                status = 'initialized'
                master_id, parent_id = self.master_id, self.parent_id
                self.remove()
                self.reset_job_id()
                self.master_id, self.parent_id = master_id, parent_id
            if repair and self.job_id and not self.status.finished:
                status = 'created'
            if status == 'initialized':
                self._run_if_new(debug=debug)
            elif status == 'created':
                que_id = self._run_if_created()
                if que_id:
                    self._logger.info('{}, status: {}, submitted: queue id {}'.format(self.job_info_str,
                                                                                      self.status,
                                                                                      que_id))
                    # print('job was submitted, queue id: ', que_id)
            elif status == 'submitted':
                self._run_if_submitted()
            elif status == 'running':
                self._run_if_running()
            elif status == 'collect':
                self._run_if_collect()
            elif status == 'suspend':
                self._run_if_suspended()
            elif status == 'refresh':
                self.run_if_refresh()
            elif status == 'busy':
                self._run_if_busy()
            elif status == 'finished':
                self._run_if_finished(run_again=run_again)
        except Exception:
            self.drop_status_to_aborted()
            raise
        except KeyboardInterrupt:
            self.drop_status_to_aborted()
            raise
        except SystemExit:
            self.drop_status_to_aborted()
            raise

    def run_if_modal(self):
        """
        The run if modal function is called by run to execute the simulation, while waiting for the output. For this we
        use subprocess.check_output()
        """
        self.run_static()
        
    def run_static(self):
        """
        The run static function is called by run to execute the simulation.
        """
        self._logger.info('{}, status: {}, run job (modal)'.format(self.job_info_str, self.status))
        if self.executable.executable_path == '':
            self.status.aborted = True
            raise ValueError('No executable set!')
        self.status.running = True
        self.project.db.item_update({"timestart": datetime.now()}, self.job_id)
        job_crashed, out = False, None
        try:
            if self.server.cores == 1 or not self.executable.mpi:
                out = subprocess.check_output(str(self.executable), cwd=self.project_hdf5.working_directory, shell=True,
                                              stderr=subprocess.STDOUT, universal_newlines=True)
            else:
                out = subprocess.check_output([self.executable.executable_path,
                                               str(self.server.cores),
                                               str(self.server.threads)],
                                              cwd=self.project_hdf5.working_directory, shell=False,
                                              stderr=subprocess.STDOUT, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            if not self.server.accept_crash:
                self._logger.warn("Job aborted")
                self._logger.warn(e.output)
                self.status.aborted = True
                error_file = posixpath.join(self.project_hdf5.working_directory, "error.msg")
                with open(error_file, "w") as f:
                    f.write(e.output)
                if self.server.run_mode.non_modal:
                    s.close_connection()
                raise RuntimeError("Job aborted")
            else:
                job_crashed = True

        self.set_input_to_read_only()
        self.status.collect = True
        self._logger.info('{}, status: {}, output: {}'.format(self.job_info_str, self.status, out))
        self.run()
        if job_crashed:
            self.status.aborted = True

    def run_if_interactive(self):
        """
        For jobs which executables are available as Python library, those can also be executed with a library call
        instead of calling an external executable. This is usually faster than a single core python job.
        """
        raise NotImplementedError("This function needs to be implemented in the specific class.")

    def interactive_close(self):
        """
        For jobs which executables are available as Python library, those can also be executed with a library call
        instead of calling an external executable. This is usually faster than a single core python job. After the
        interactive execution, the job can be closed using the interactive_close function.
        """
        raise NotImplementedError("This function needs to be implemented in the specific class.")

    def interactive_fetch(self):
        """
        For jobs which executables are available as Python library, those can also be executed with a library call
        instead of calling an external executable. This is usually faster than a single core python job. To access the
        output data during the execution the interactive_fetch function is used.
        """
        raise NotImplementedError("This function needs to be implemented in the specific class.")

    def interactive_flush(self, path="generic", include_last_step=True):
        """
        For jobs which executables are available as Python library, those can also be executed with a library call
        instead of calling an external executable. This is usually faster than a single core python job. To write the
        interactive cache to the HDF5 file the interactive flush function is used.
        """
        raise NotImplementedError("This function needs to be implemented in the specific class.")

    def run_if_interactive_non_modal(self):
        """
        For jobs which executables are available as Python library, those can also be executed with a library call
        instead of calling an external executable. This is usually faster than a single core python job.
        """
        raise NotImplementedError("This function needs to be implemented in the specific class.")

    # def run_if_non_modal(self):
    #     """
    #     The run if non modal function is called by run to execute the simulation in the background. For this we use
    #     subprocess.Popen()
    #     """
    #     shell = (os.name == 'nt')
    #     try:
    #         file_name = posixpath.join(self.project_hdf5.working_directory, "run_job.py")
    #         self._logger.info("{}, status: {}, script: {}".format(self.job_info_str, self.status, file_name))
    #         with open(posixpath.join(self.project_hdf5.working_directory, 'out.txt'), mode='w') as f_out:
    #             with open(posixpath.join(self.project_hdf5.working_directory, 'error.txt'), mode='w') as f_err:
    #                 self._process = subprocess.Popen(['python', '-m', 'pyiron.base.job.wrappercmd', '-p',
    #                                                   self.working_directory, '-j', str(self.job_id)],
    #                                                  cwd=self.project_hdf5.working_directory, shell=shell, stdout=f_out,
    #                                                  stderr=f_err, universal_newlines=True)
    #         self._logger.info("{}, status: {}, job submitted".format(self.job_info_str, self.status))
    #     except subprocess.CalledProcessError as e:
    #         self._logger.warn("Job aborted")
    #         self._logger.warn(e.output)
    #         self.status.aborted = True
    #         raise ValueError("run_job.py crashed")
    #     s.logger.info('submitted run %s', self.job_name)
    #     self._logger.info('job status: %s', self.status)

    def run_if_non_modal(self):
        """
        The run if non modal function is called by run to execute the simulation in the background. For this we use
        multiprocessing.Process()
        """
        p = multiprocessing.Process(target=multiprocess_wrapper, args=(self.job_id,
                                                                       self.project_hdf5.working_directory,
                                                                       False))
        if self.master_id:
            del self
        p.start()

    def run_if_manually(self, _manually_print=True):
        """
        The run if manually function is called by run if the user decides to execute the simulation manually - this
        might be helpful to debug a new job type or test updated executables.
        
        Args:
            _manually_print (bool): Print explanation how to run the simulation manually - default=True.
        """
        if _manually_print:
            abs_working = posixpath.abspath(self.project_hdf5.working_directory)
            print('You have selected to start the job manually. ' +
                  'To run it, go into the working directory {} and '.format(abs_working) +
                  'call \'python -m pyiron.base.job.wrappercmd -p {}'.format(abs_working) +
                  ' -j {} \' '.format(self.job_id))

    def run_if_scheduler(self):
        """
        The run if queue function is called by run if the user decides to submit the job to and queing system. The job
        is submitted to the queuing system using subprocess.Popen()

        Returns:
            int: Returns the queue ID for the job.
        """
        if s.queue_adapter is None:
            raise TypeError('No queue adapter defined.')
        try:
            que_id = s.queue_adapter.submit_job(queue=self.server.queue,
                                                job_name='pi_' + str(self.job_id),
                                                working_directory=self.project_hdf5.working_directory,
                                                cores=self.server.cores,
                                                run_time_max=self.server.run_time,
                                                memory_max=self.server.memory_limit,
                                                command='python -m pyiron.base.job.wrappercmd -p ' +
                                                        self.working_directory + ' -j ' + str(self.job_id))
            self.server.queue_id = que_id
            self._server.to_hdf(self._hdf5)
            print('Queue system id: ', que_id)
        except subprocess.CalledProcessError as e:
            self._logger.warn("Job aborted")
            self._logger.warn(e.output)
            self.status.aborted = True
            raise ValueError("run_queue.sh crashed")
        s.logger.debug('submitted %s', self.job_name)
        self._logger.debug('job status: %s', self.status)
        return que_id

    def send_to_database(self):
        """
        if the jobs should be store in the external/public database this could be implemented here, but currently it is
        just a placeholder.
        """
        if self.server.send_to_db:
            pass

    def create_job(self, job_type, job_name):
        """
        Create one of the following jobs:
        - 'StructureContainer’:
        - ‘StructurePipeline’:
        - ‘AtomisticExampleJob’: example job just generating random number
        - ‘ExampleJob’: example job just generating random number
        - ‘Lammps’:
        - ‘KMC’:
        - ‘Sphinx’:
        - ‘Vasp’:
        - ‘GenericMaster’:
        - ‘SerialMaster’: series of jobs run in serial
        - ‘AtomisticSerialMaster’:
        - ‘ParallelMaster’: series of jobs run in parallel
        - ‘KmcMaster’:
        - ‘ThermoLambdaMaster’:
        - ‘RandomSeedMaster’:
        - ‘MeamFit’:
        - ‘Murnaghan’:
        - ‘MinimizeMurnaghan’:
        - ‘ElasticMatrix’:
        - ‘ConvergenceVolume’:
        - ‘ConvergenceEncutParallel’:
        - ‘ConvergenceKpointParallel’:
        - ’PhonopyMaster’:
        - ‘DefectFormationEnergy’:
        - ‘LammpsASE’:
        - ‘PipelineMaster’:
        - ’TransformationPath’:
        - ‘ThermoIntEamQh’:
        - ‘ThermoIntDftEam’:
        - ‘ScriptJob’: Python script or jupyter notebook job container
        - ‘ListMaster': list of jobs

        Args:
            job_type (str): job type can be ['StructureContainer’, ‘StructurePipeline’, ‘AtomisticExampleJob’,
                                             ‘ExampleJob’, ‘Lammps’, ‘KMC’, ‘Sphinx’, ‘Vasp’, ‘GenericMaster’,
                                             ‘SerialMaster’, ‘AtomisticSerialMaster’, ‘ParallelMaster’, ‘KmcMaster’,
                                             ‘ThermoLambdaMaster’, ‘RandomSeedMaster’, ‘MeamFit’, ‘Murnaghan’,
                                             ‘MinimizeMurnaghan’, ‘ElasticMatrix’, ‘ConvergenceVolume’,
                                             ‘ConvergenceEncutParallel’, ‘ConvergenceKpointParallel’, ’PhonopyMaster’,
                                             ‘DefectFormationEnergy’, ‘LammpsASE’, ‘PipelineMaster’,
                                             ’TransformationPath’, ‘ThermoIntEamQh’, ‘ThermoIntDftEam’, ‘ScriptJob’,
                                             ‘ListMaster']
            job_name (str): name of the job

        Returns:
            GenericJob: job object depending on the job_type selected
        """
        job = self.project.create_job(job_type=job_type, job_name=job_name)
        if static_isinstance(obj=job.__class__,
                             obj_type=['pyiron.base.master.parallel.ParallelMaster',
                                       'pyiron.base.master.serial.SerialMasterBase',
                                       'pyiron.atomistic.job.interactivewrapper.InteractiveWrapper']):
            job.ref_job = self
            if self.server.run_mode.non_modal:
                job.server.run_mode.non_modal = True
            elif self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal:
                job.server.run_mode.interactive = True
        return job

    def update_master(self):
        """
        After a job is finished it checks whether it is linked to any metajob - meaning the master ID is pointing to
        this jobs job ID. If this is the case and the master job is in status suspended - the child wakes up the master
        job, sets the status to refresh and execute run on the master job. During the execution the master job is set to
        status refresh. If another child calls update_master, while the master is in refresh the status of the master is
        set to busy and if the master is in status busy at the end of the update_master process another update is
        triggered.
        """
        master_id = self.master_id
        project = self.project
        self._logger.info("update master: {} {}".format(master_id, self.get_job_id()))
        if master_id is not None and not self.server.run_mode.modal and not self.server.run_mode.interactive:
            queue_flag = self.server.run_mode.queue
            master_db_entry = project.db.get_item_by_id(master_id)
            if master_db_entry['status'] == 'suspended':
                project.db.item_update({'status': 'refresh'}, master_id)
                self._logger.info("run_if_refresh() called")
                # p = multiprocessing.Process(target=multiprocess_master, args=(master_id,
                #                                                               self.project.path,
                #                                                               self.server.run_mode.thread,
                #                                                               False))
                # del self
                # p.start()
                del self
                master_inspect = project.inspect(master_id)
                if master_inspect["server"]["run_mode"] == "non_modal" or \
                        (master_inspect["server"]["run_mode"] == "modal" and queue_flag):
                    master = project.load(master_id)
                    # master = master_inspect.load_object()
                    master.run_if_refresh()
                # if master.server.run_mode.non_modal or master.server.run_mode.queue:
                #     master._run_if_refresh()
                #     if master.server.run_mode.queue and master._process:
                #         master._process.communicate()
            elif master_db_entry['status'] == 'refresh':
                project.db.item_update({'status': 'busy'}, master_id)
                self._logger.info("busy master: {} {}".format(master_id, self.get_job_id()))
                del self

    def job_file_name(self, file_name, cwd=None):
        """
        combine the file name file_name with the path of the current working directory

        Args:
            file_name (str): name of the file
            cwd (str): current working directory - this overwrites self.project_hdf5.working_directory - optional

        Returns:
            str: absolute path to the file in the current working directory
        """
        if not cwd:
            cwd = self.project_hdf5.working_directory
        return posixpath.join(cwd, file_name)

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the GenericJob in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        if hdf is not None:
            self._hdf5 = hdf
        if group_name is not None:
            self._hdf5.open(group_name)
        self._executable_activate_mpi()
        self._type_to_hdf()
        self._server.to_hdf(self._hdf5)
        with self._hdf5.open('input') as hdf_input:
            generic_dict = {"restart_file_list": self._restart_file_list,
                            "restart_file_dict": self._restart_file_dict,
                            "exclude_nodes_hdf": self._exclude_nodes_hdf,
                            "exclude_groups_hdf": self._exclude_groups_hdf}
            hdf_input["generic_dict"] = generic_dict

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the GenericJob from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        if hdf is not None:
            self._hdf5 = hdf
        if group_name is not None:
            self._hdf5.open(group_name)
        self._type_from_hdf()
        self._server.from_hdf(self._hdf5)
        with self._hdf5.open('input') as hdf_input:
            if "generic_dict" in hdf_input.list_nodes():
                generic_dict = hdf_input["generic_dict"]
                self._restart_file_list = generic_dict["restart_file_list"]
                self._restart_file_dict = generic_dict["restart_file_dict"]
                self._exclude_nodes_hdf = generic_dict["exclude_nodes_hdf"]
                self._exclude_groups_hdf = generic_dict["exclude_groups_hdf"]
            # Backwards compatbility
            if "restart_file_list" in hdf_input.list_nodes():
                self._restart_file_list = hdf_input["restart_file_list"]
            if "restart_file_dict" in hdf_input.list_nodes():
                self._restart_file_dict = hdf_input["restart_file_dict"]
            if "exclude_nodes_hdf" in hdf_input.list_nodes():
                self._exclude_nodes_hdf = hdf_input["exclude_nodes_hdf"]
            if "exclude_groups_hdf" in hdf_input.list_nodes():
                self._exclude_groups_hdf = hdf_input["exclude_groups_hdf"]

    def save(self):
        """
        Save the object, by writing the content to the HDF5 file and storing an entry in the database.

        Returns:
            (int): Job ID stored in the database
        """
        self.to_hdf()
        job_id = self.project.db.add_item_dict(self.db_entry())
        self._job_id = job_id
        self.refresh_job_status()
        return job_id

    def convergence_check(self):
        """
        Validate the convergence of the calculation.

        Returns:
             (bool): If the calculation is converged
        """
        return True

    def db_entry(self):
        """
        Generate the initial database entry for the current GenericJob

        Returns:
            (dict): database dictionary {"username", "projectpath", "project", "job", "subjob", "hamversion",
                                         "hamilton", "status", "computer", "timestart", "masterid", "parentid"}
        """
        db_dict = {"username": s.login_user,
                   "projectpath": self.project_hdf5.root_path,
                   "project": self.project_hdf5.project_path,
                   "job": self.job_name,
                   "subjob": self.project_hdf5.h5_path,
                   "hamversion": self.version,
                   "hamilton": self.__name__,
                   "status": self.status.string,
                   "computer": self._db_server_entry(),
                   "timestart": datetime.now(),
                   "masterid": self.master_id,
                   "parentid": self.parent_id}
        return db_dict

    def restart(self, snapshot=-1, job_name=None, job_type=None):
        """
        Create an restart calculation from the current calculation - in the GenericJob this is the same as create_job().
        A restart is only possible after the current job has finished. If you want to run the same job again with
        different input parameters use job.run(run_again=True) instead.

        Args:
            snapshot (int): time step from which to restart the calculation - default=-1 - the last time step
            job_name (str): job name of the new calculation - default=<job_name>_restart
            job_type (str): job type of the new calculation - default is the same type as the exeisting calculation

        Returns:

        """
        if not self.job_id:
            self._create_job_structure(debug=False)
        if job_name is None:
            job_name = "{}_restart".format(self.job_name)
        if job_type is None:
            job_type = self.__name__
        if job_type == self.__name__:
            new_ham = self.copy_to(new_job_name=job_name, new_database_entry=False)
        else:
            new_ham = self.create_job(job_type, job_name)
        new_ham.parent_id = self.job_id
        # ensuring that the new job does not inherit the restart_file_list from the old job
        new_ham._restart_file_list = list()
        new_ham._restart_file_dict = dict()
        return new_ham

    def _copy_restart_files(self):
        """
        Internal helper function to copy the files required for the restart job.
        """
        if not (os.path.isdir(self.working_directory)):
            raise ValueError("The working directory is not yet available to copy restart files")
        for i, actual_name in enumerate([os.path.basename(f) for f in self._restart_file_list]):
            if actual_name in self.restart_file_dict.keys():
                new_name = self.restart_file_dict[actual_name]
                shutil.copy(self.restart_file_list[i], posixpath.join(self.working_directory, new_name))
            else:
                shutil.copy(self.restart_file_list[i], self.working_directory)

    def _run_manually(self, _manually_print=True):
        """
        Internal helper function to run a job manually.

        Args:
            _manually_print (bool): [True/False] print command for execution - default=True
        """
        if _manually_print:
            print('You have selected to start the job manually. ' +
                  'To run it, go into the working directory {} and ' +
                  'call \'python run_job.py\' '.format(posixpath.abspath(self.project_hdf5.working_directory)))

    def _run_if_new(self, debug=False):
        """
        Internal helper function the run if new function is called when the job status is 'initialized'. It prepares
        the hdf5 file and the corresponding directory structure.

        Args:
            debug (bool): Debug Mode
        """
        self.validate_ready_to_run()
        if self.server.run_mode.queue:
            self.check_setup()
        if self.check_if_job_exists():
            print('job exists already and therefore was not created!')
        else:
            self._create_job_structure(debug=debug)
            self.run()

    def _run_if_created(self):
        """
        Internal helper function the run if created function is called when the job status is 'created'. It executes
        the simulation, either in modal mode, meaning waiting for the simulation to finish, manually, or submits the
        simulation to the que.

        Returns:
            int: Queue ID - if the job was send to the queue
        """
        self.status.submitted = True

        # Different run modes
        if self.server.run_mode.manual:
            self.run_if_manually()
        elif self.server.run_mode.modal:
            self.run_static()
        elif self.server.run_mode.non_modal or self.server.run_mode.thread:
            self.run_if_non_modal()
        elif self.server.run_mode.queue:
            return self.run_if_scheduler()
        elif self.server.run_mode.interactive:
            self.run_if_interactive()
        elif self.server.run_mode.interactive_non_modal:
            self.run_if_interactive_non_modal()
        return None

    def _run_if_submitted(self):  # Submitted jobs are handled by the job wrapper!
        """
        Internal helper function the run if submitted function is called when the job status is 'submitted'. It means
        the job is waiting in the queue. ToDo: Display a list of the users jobs in the queue.
        """
        if self.server.run_mode.queue and not self.project.queue_check_job_is_waiting_or_running(self.job_id):
            self.run(run_again=True)
        else:
            print('Job ' + str(self.job_id) + ' is waiting in the que!')

    def _run_if_running(self):
        """
        Internal helper function the run if running function is called when the job status is 'running'. It allows the
        user to interact with the simulation while it is running.
        """
        if self.server.run_mode.queue and not self.project.queue_check_job_is_waiting_or_running(self.job_id):
            self.run(run_again=True)
        elif self.server.run_mode.interactive:
            self.run_if_interactive()
        elif self.server.run_mode.interactive_non_modal:
            self.run_if_interactive_non_modal()
        else:
            print('Job ' + str(self.job_id) + ' is running!')

    def run_if_refresh(self):
        """
        Internal helper function the run if refresh function is called when the job status is 'refresh'. If the job was
        suspended previously, the job is going to be started again, to be continued.
        """
        raise NotImplementedError('Refresh is not supported for this job type for job  ' + str(self.job_id))

    def set_input_to_read_only(self):
        """
        This function enforces read-only mode for the input classes, but it has to be implement in the individual
        classes.
        """
        pass

    def _run_if_busy(self):
        """
        Internal helper function the run if busy function is called when the job status is 'busy'.
        """
        raise NotImplementedError('Refresh is not supported for this job type for job  ' + str(self.job_id))

    def _run_if_collect(self):
        """
        Internal helper function the run if collect function is called when the job status is 'collect'. It collects
        the simulation output using the standardized functions collect_output() and collect_logfiles(). Afterwards the
        status is set to 'finished'
        """
        self.collect_output()
        self.collect_logfiles()
        self.project.db.item_update(self._runtime(), self.job_id)
        if self.status.collect:
            if not self.convergence_check():
                self.status.not_converged = True
            else:
                if self._compress_by_default:
                    self.compress()
                self.status.finished = True
        self._calculate_successor()
        self.send_to_database()
        self.update_master()

    def _run_if_suspended(self):
        """
        Internal helper function the run if suspended function is called when the job status is 'suspended'. It
        restarts the job by calling the run if refresh function after setting the status to 'refresh'.
        """
        self.status.refresh = True
        self.run()

    def _run_if_finished(self, run_again=False):
        """
        Internal helper function the run if finished function is called when the job status is 'finished'. It loads
        the existing job.

        Args:
            run_again (bool): Delete the existing job and run the simulation again.
        """
        if run_again:
            parent_id = self.parent_id
            self.parent_id = None
            self.remove()
            self._job_id = None
            self.status.initialized = True
            self.parent_id = parent_id
            self.run()
        else:
            self.from_hdf()

    def _executable_activate(self, enforce=False):
        """
        Internal helper function to koad the executable object, if it was not loaded already.
        """
        if not self._executable or enforce:
            if len(self.__module__.split('.')) > 1:
                self._executable = Executable(codename=self.__name__,
                                              module=self.__module__.split('.')[1],
                                              path_binary_codes=s.resource_paths)
            else:
                self._executable = Executable(codename=self.__name__, path_binary_codes=s.resource_paths)

    def _type_to_hdf(self):
        """
        Internal helper function to save type and version in HDF5 file root
        """
        self._hdf5["NAME"] = self.__name__
        self._hdf5["TYPE"] = str(type(self))
        if self._executable:
            self._hdf5["VERSION"] = self.executable.version
        else:
            self._hdf5["VERSION"] = self.__version__

    def _type_from_hdf(self):
        """
        Internal helper function to load type and version from HDF5 file root
        """
        self.__obj_type__ = self._hdf5["TYPE"]
        if self._executable:
            try:
                self.executable.version = self._hdf5["VERSION"]
            except ValueError:
                self.executable.executable_path = self._hdf5["VERSION"]
        else:
            self.__obj_version__ = self._hdf5["VERSION"]

    def _runtime(self):
        """
        Internal helper function to calculate runtime by substracting the starttime, from the stoptime.

        Returns:
            (dict): Database dictionary db_dict
        """
        start_time = self.project.db.get_item_by_id(self.job_id)["timestart"]
        stop_time = datetime.now()
        return {"timestop": stop_time, "totalcputime": int((stop_time - start_time).total_seconds())}

    def _db_server_entry(self):
        """
        Internal helper function to connect all the info regarding the server into a single word that can be used
        e.g. as entry in a database

        Returns:
            (str): server info as single word

        """
        return self._server.db_entry()

    def _executable_activate_mpi(self):
        """
        Internal helper function to switch the executable to MPI mode
        """
        try:
            if self.server.cores > 1:
                self.executable.mpi = True
        except ValueError:
            self.server.cores = 1
            warnings.warn('No multi core executable found falling back to the single core executable.', RuntimeWarning)

    def _calculate_predecessor(self):
        """
        Internal helper function to calculate the predecessor of the current job if it was not calculated before. This
        function is used to execute a series of jobs based on their parent relationship - marked by the parent ID.
        Mainly used by the ListMaster job type.
        """
        parent_id = self.parent_id
        if parent_id is not None:
            if self._hdf5.db.get_item_by_id(parent_id)['status'] in ['initialized', 'created']:
                self.status.suspended = True
                parent_job = self._hdf5.load(parent_id)
                parent_job.run()

    def _calculate_successor(self):
        """
        Internal helper function to calculate the successor of the current job. This function is used to execute a
        series of jobs based on their parent relationship - marked by the parent ID. Mainly used by the ListMaster job
        type.
        """
        for child_id in sorted([job["id"] for job in self.project.db.get_items_dict({'parentid': str(self.job_id)},
                                                                                    return_all_columns=False)]):
            if self._hdf5.db.get_item_by_id(child_id)['status'] in ['suspended']:
                child = self._hdf5.load(child_id)
                child.status.created = True
                self._before_successor_calc(child)
                child.run()

    def _create_job_structure(self, debug=False):
        """
        Internal helper function to create the input directories, save the job in the database and write the wrapper.

        Args:
            debug (bool): Debug Mode
        """
        self._job_id = self.save()
        print('The job ' + self.job_name + ' was saved and received the ID: ' + str(self._job_id))
        if self._check_if_input_should_be_written():
            self.project_hdf5.create_working_directory()
            self.write_input()
            self._copy_restart_files()
        self.status.created = True
        self._calculate_predecessor()

    def _check_if_input_should_be_written(self):
        return not (self.server.run_mode.interactive or self.server.run_mode.interactive_non_modal)

    def _before_successor_calc(self, ham):
        """
        Internal helper function which is executed based on the hamiltonian of the successor job, before it is executed.
        This function is used to execute a series of jobs based on their parent relationship - marked by the parent ID.
        Mainly used by the ListMaster job type.
        """
        pass


def multiprocess_wrapper(job_id, working_dir, debug=False):
    job_wrap = JobWrapper(working_directory=str(working_dir), job_id=int(job_id), debug=debug)
    job_wrap.job.run_static()


# def multiprocess_master(job_id, working_dir, is_thread_mode=False, debug=False):
#     job_wrap = JobWrapper(working_directory=str(working_dir), job_id=int(job_id), debug=debug)
#     job_wrap.job._run_if_refresh()
#     if is_thread_mode and job_wrap.job._process:
#         job_wrap.job._process.communicate()
