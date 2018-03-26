# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from collections import OrderedDict
from datetime import datetime
import pandas
from pyironbase.objects.job.generic import GenericJob
from pyironbase.objects.job.master import GenericMaster
from pyironbase.objects.job.submissionstatus import SubmissionStatus
from pyironbase.objects.generic.parameters import GenericParameters
from pyironbase.objects.job.jobstatus import JobStatus

"""
The parallel master class is a metajob consisting of a list of jobs which are executed in parallel. 
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class ParallelMaster(GenericMaster):
    """
    MasterJob that handles the creation and analysis of several parallel jobs (including master and
    continuation jobs), Examples are Murnaghan or Phonon calculations

    Args:
        project (ProjectHDFio): ProjectHDFio instance which points to the HDF5 file the job is stored in
        job_name (str): name of the job, which has to be unique within the project

    Attributes:

        .. attribute:: job_name

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

        .. attribute:: child_names

            Dictionary matching the child ID to the child job name.

        .. attribute:: ref_job

            Reference job template from which all jobs within the ParallelMaster are generated.

        .. attribute:: number_jobs_total

            Total number of jobs
    """
    def __init__(self, project, job_name):
        self.input = GenericParameters("parameters")
        super(ParallelMaster, self).__init__(project, job_name=job_name)
        self.__name__ = "ParallelMaster"
        self.__version__ = "0.3"
        self._ref_job = None
        self._output = GenericOutput()
        self.submission_status = SubmissionStatus(db=project.db, job_id=self.job_id)
        self.refresh_submission_status()

    @property
    def ref_job(self):
        """
        Get the reference job template from which all jobs within the ParallelMaster are generated.

        Returns:
            GenericJob: reference job
        """
        if self._ref_job:
            return self._ref_job
        try:
            if isinstance(self[0], GenericJob):
                self._ref_job = self[0]
                self._ref_job._job_id = None
                self._ref_job._status = JobStatus(db=self.project.db)
                return self._ref_job
            else:
                return None
        except IndexError:
            return None

    @ref_job.setter
    def ref_job(self, ref_job):
        """
        Set the reference job template from which all jobs within the ParallelMaster are generated.

        Args:
            ref_job (GenericJob): reference job
        """
        self.append(ref_job)

    @property
    def number_jobs_total(self):
        """
        Get number of total jobs

        Returns:
            int: number of total jobs
        """
        return self.submission_status.total_jobs

    @number_jobs_total.setter
    def number_jobs_total(self, num_jobs):
        """
        Set number of total jobs (optional: default = None)

        Args:
            num_jobs (int): number of submitted jobs
        """
        self.input['num_points'] = num_jobs
        self.submission_status.total_jobs = num_jobs

    def reset_job_id(self, job_id=None):
        """
        Reset the job id sets the job_id to None as well as all connected modules like JobStatus and SubmissionStatus.
        """
        super(ParallelMaster, self).reset_job_id(job_id=job_id)
        self.submission_status = SubmissionStatus(db=self.project.db, job_id=self.job_id)

    def to_hdf(self, hdf=None, group_name=None):
        """
        Store the ParallelMaster in an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(ParallelMaster, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open('input') as hdf5_input:
            self.input.to_hdf(hdf5_input)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore the ParallelMaster from an HDF5 file

        Args:
            hdf (ProjectHDFio): HDF5 group object - optional
            group_name (str): HDF5 subgroup name - optional
        """
        super(ParallelMaster, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open('input') as hdf5_input:
            self.input.from_hdf(hdf5_input)

    def write_input(self):
        """
        Write the input files - this contains the GenericInput of the ParallelMaster as well as reseting the submission
        status.
        """
        self.submission_status.total_jobs = self.input['num_points']
        self.submission_status.submitted_jobs = 0
        self.input.write_file(file_name="input.inp", cwd=self.working_directory)

    def collect_output(self):
        """
        Collect the output files of the external executable and store the information in the HDF5 file. This method has
        to be implemented in the individual meta jobs derived from the ParallelMaster.
        """
        raise ValueError("Implement in derived class")

    def collect_logfiles(self):
        """
        Collect the log files of the external executable and store the information in the HDF5 file. This method is
        currently not implemented for the ParallelMaster.
        """
        pass

    def output_to_pandas(self, sort_by=None, h5_path="output"):
        """
        Convert output of all child jobs to a pandas Dataframe object.

        Args:
            sort_by (str): sort the output using pandas.DataFrame.sort_values(by=sort_by)
            h5_path (str): select child output to include - default='output'

        Returns:
            pandas.Dataframe: output as dataframe
        """
        # TODO: The output to pandas function should no longer be required
        with self.project_hdf5.open(h5_path) as hdf:
            for key in hdf.list_nodes():
                self._output[key] = hdf[key]
        df = pandas.DataFrame(self._output)
        if sort_by is not None:
            df = df.sort_values(by=sort_by)
        return df

    # TODO: make it more general and move it then into genericJob
    def show_hdf(self):
        """
        Display the output of the child jobs in a human readable print out
        """
        try:
            from IPython import display
            for nn in self.project_hdf5.list_groups():
                with self.project_hdf5.open(nn) as hdf_dir:
                    display.display(nn)
                    if nn.strip() == 'output':
                        display.display(self.output_to_pandas(h5_path=nn))
                        continue
                    for n in hdf_dir.list_groups():
                        display.display("-->" + n)
                        try:
                            display.display(hdf_dir.get_pandas(n))
                        except Exception as e:
                            print(e)
                            print("Not a pandas object")
        except ImportError:
            print('show_hdf() requires IPython to be installed.')

    def save(self):
        """
        Save the object, by writing the content to the HDF5 file and storing an entry in the database.

        Returns:
            (int): Job ID stored in the database
        """
        job_id = super(ParallelMaster, self).save()
        self.refresh_submission_status()
        return job_id

    def refresh_submission_status(self):
        """
        Refresh the submission status - if a job ID job_id is set then the submission status is loaded from the database.
        """
        if self.job_id:
            self.submission_status = SubmissionStatus(db=self.project.db, job_id=self.job_id)
            self.submission_status.refresh()

    def copy(self):
        """
        Copy the GenericJob object which links to the job and its HDF5 file

        Returns:
            GenericJob: New GenericJob object pointing to the same job
        """
        new_job = super(ParallelMaster, self).copy()
        new_job.ref_job = self.ref_job
        return new_job

    def copy_to(self, project, new_job_name=None, input_only=False, new_database_entry=True):
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
        new_generic_job = super(ParallelMaster, self).copy_to(project, new_job_name=new_job_name,
                                                              input_only=input_only,
                                                              new_database_entry=new_database_entry)
        new_generic_job.submission_status = SubmissionStatus(db=new_generic_job._hdf5.project.db,
                                                             job_id=new_generic_job.job_id)
        return new_generic_job

    def create_jobs(self):
        """
        The create_jobs method needs to be implemented in the derived classes to define the rules how to create the job
        series which should then be executed in parallel.
        """
        raise NotImplementedError('define in derived class to create jobs')

    def is_finished(self):
        """
        Check if the ParallelMaster job is finished - by checking the job status and the submission status.

        Returns:
            bool: [True/False]
        """
        if self.status.finished:
            return True
        # self.status.busy = True
        self.submission_status.refresh()
        if not self.submission_status.finished:
            return False
        else:
            return set([self.project.db.get_item_by_id(child_id)['status'] for child_id in self.child_ids])\
                   <{'finished', 'busy', 'refresh'}

    def run_if_modal(self):
        """
        The run if modal function is called by run to execute the simulation, while waiting for the output. For the
        ParallelMaster this means executing all the childs appened in parallel.
        """
        self._logger.info('{} run parallel master (modal)'.format(self.job_info_str))
        self.status.running = True
        self.create_jobs()
        if self.is_finished():
            self.status.collect = True
            self.run()

    def iter_jobs(self, convert_to_object=True):
        """
        Iterate over the jobs within the ListMaster

        Args:
            convert_to_object (bool): load the full GenericJob object (default) or just the HDF5 / JobCore object

        Returns:
            yield: Yield of GenericJob or JobCore
        """
        for job_id in self.child_ids:
            yield self.project.load(job_id, convert_to_object=convert_to_object)

    def __getitem__(self, item):
        """
        Get/ read data from the HDF5 file

        Args:
            item (str, slice): path to the data or key of the data object

        Returns:
            dict, list, float, int: data or data object
        """
        child_id_lst = self.child_ids
        child_name_lst = [self.project.db.get_item_by_id(child_id)["job"] for child_id in self.child_ids]
        if isinstance(item, str):
            name_lst = item.split("/")
            item_obj = name_lst[0]
            if item_obj in child_name_lst:
                child_id = child_id_lst[child_name_lst.index(item_obj)]
                child = self.project.load(child_id, convert_to_object=True)
                if len(name_lst) == 1:
                    return child
                else:
                    return child['/'.join(name_lst[1:])]
            if item_obj in self._job_list:
                child = getattr(self, item_obj)
                if len(name_lst) == 1:
                    return child
                else:
                    return child['/'.join(name_lst[1:])]
            return super(GenericMaster, self).__getitem__(item)
        elif isinstance(item, int):
            total_lst = self._job_list + child_name_lst
            job_name = total_lst[item]
            if job_name in self._job_list:
                job_name = self._job_list[item]
                return getattr(self, job_name)
            else:
                child_id = child_id_lst[child_name_lst.index(job_name)]
                return self.project.load(child_id, convert_to_object=True)

    def __len__(self):
        """
        Length of the ListMaster equal the number of childs appended.

        Returns:
            int: length of the ListMaster
        """
        return len(self.child_ids)

    def _run_if_refresh(self):
        """
        Internal helper function the run if refresh function is called when the job status is 'refresh'. If the job was
        suspended previously, the job is going to be started again, to be continued.
        """
        self._logger.info("{}, status: {}, finished: {} parallel master "
                          "refresh".format(self.job_info_str, self.status, self.is_finished()))
        if self.is_finished():
            self.status.collect = True
            self.run()  # self.run_if_collect()
        elif (self.server.run_mode.non_modal or self.server.run_mode.queue) and not self.submission_status.finished:
            self.run_if_modal()
        else:
            self.refresh_job_status()
            if self.status.refresh:
                self.status.suspended = True

    def _run_if_collect(self):
        """
        Internal helper function the run if collect function is called when the job status is 'collect'. It collects
        the simulation output using the standardized functions collect_output() and collect_logfiles(). Afterwards the
        status is set to 'finished'.
        """
        self._logger.info("{}, status: {}, finished".format(self.job_info_str, self.status))
        self.collect_output()

        job_id = self.get_job_id()
        db_dict = {}
        start_time = self.project.db.get_item_by_id(job_id)["timestart"]
        db_dict["timestop"] = datetime.now()
        db_dict["totalcputime"] = (db_dict["timestop"] - start_time).seconds
        self.project.db.item_update(db_dict, job_id)
        self.status.finished = True
        self.update_master()
        self.send_to_database()
        self._logger.info("{}, status: {}, parallel master".format(self.job_info_str, self.status))

    def _create_child_job(self, job_name):
        """
        Internal helper function to create the next child job from the reference job template - usually this is called
        as part of the create_jobs() function.

        Args:
            job_name (str): name of the next job

        Returns:
            GenericJob: next job
        """
        project = self.project.open(self.job_name + '_hdf5')
        job_id = project.get_job_id(job_name)
        if job_id is not None:
            ham = project.load(job_id)
            print("job ", job_name, " found, status:", ham.status)
            if ham.server.run_mode.queue:
                self.project.refresh_job_status_based_on_job_id(job_id, que_mode=True)
            else:
                self.project.refresh_job_status_based_on_job_id(job_id, que_mode=False)
            if ham.status.aborted:
                ham.status.created = True

            print("job - status:", ham.status)
            return ham

        job = self.ref_job.copy()

        if self.server.new_hdf:
            job._hdf5 = self.project_hdf5.create_hdf(path=self._hdf5._project.open(self.job_name + '_hdf5').path,
                                                     job_name=job_name)
        else:
            job._hdf5 = self.project_hdf5.open(job_name)
        try:
            self.ref_job.project_hdf5.copy_to(job._hdf5, maintain_name=False)
        except ValueError:
            pass
        self._logger.debug(
            "create_job:: {} {} {} {}".format(self.project_hdf5.path, self._name, self.project_hdf5.h5_path,
                                              str(self.get_job_id())))
        job._name = job_name
        job.master_id = self.get_job_id()
        job.status.initialized = True
        if self.server.run_mode.non_modal and job.server.run_mode.modal:
            job.server.run_mode.non_modal = True
        elif self.server.run_mode.queue:
            job.server.run_mode.thread = True
        return job

    def _db_server_entry(self):
        """
        connect all the info regarding the server into a single word that can be used e.g. as entry in a database

        Returns:
            (str): server info as single word

        """
        db_entry = super(ParallelMaster, self)._db_server_entry()
        if self.submission_status.total_jobs:
            return db_entry + '#' + str(self.submission_status.submitted_jobs) + '/' + \
                   str(self.submission_status.total_jobs)
        else:
            return db_entry + '#' + str(self.submission_status.submitted_jobs)


class GenericOutput(OrderedDict):
    """
    Generic Output just a place holder to store the output of the last child directly in the ParallelMaster.
    """
    def __init__(self):
        super(GenericOutput, self).__init__()
