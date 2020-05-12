# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.base.database.generic import DatabaseAccess
from pyiron.base.database.filetable import FileTable

"""
The SubmissionStatus class belongs to the GenericJob object. It is presently used only for the parallel master class.
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class SubmissionStatus(object):
    """
    The SubmissionStatus object handles the different submission states a job could have. The available states are:
        initialized: No jobs have been submitted.
        sub_m_n: m out of n jobs have been submitted
        finished: The job and all connected sub jobs are finished.

    Args:
        initial_status (str): If no initial status is provided the status is set to 'initialized'
        db (DatabaseAccess): The database which is responsible for this job.
        job_id: job ID

    Attributes:

        .. attribute:: database

            the database which is responsible for this job.

        .. attribute:: job_id

            Job ID

        .. attribute:: string

            job status as string

        .. attribute:: total_jobs

            number of jobs which have to be submitted in total

        .. attribute:: submitted_jobs

            number of jobs which have been submitted
    """

    STATUS = ["initialized", "finished"]

    def __init__(self, initial_status="initialized", db=None, job_id=None):
        self._submitted_jobs = 0
        self._total_jobs = None
        self._string = initial_status
        self.database = db
        self.job_id = job_id

    @property
    def database(self):
        """
        Get the database which is responsible for this job. If no database is linked it returns None.

        Returns:
            DatabaseAccess: The database which is responsible for this job.
        """
        return self._db

    @database.setter
    def database(self, db):
        """
        Set the database which is responsible for this job.

        Args:
            db (DatabaseAccess): The database which should be responsible for this job.
        """
        if db and not isinstance(db, (DatabaseAccess, FileTable)):
            raise TypeError("The database has to be an DatabaseAccess object.")
        self._db = db

    @property
    def initialized(self):
        """
        Check if the status is 'initialized', meaning the object for the corresponding job was just created.

        Returns:
            bool: [True/False]
        """
        return self.string == "initialized"

    @property
    def submitted(self):
        """
        Check if the status is 'submitted', meaning the job has not yet submitted all jobs.

        Returns:
            bool: [True/False]
        """
        return "submitted_" in self.string

    @property
    def submitted_jobs(self):
        """
        Get the number of jobs which have been submitted.

        Returns:
            int: number of submitted jobs
        """
        self.refresh()
        return self._submitted_jobs

    @submitted_jobs.setter
    def submitted_jobs(self, num_submitted):
        """
        Set the status to 'submitted', meaning some but not all jobs have been submitted

        Args:
            num_submitted (int): number of submitted jobs
        """
        self._submitted_jobs = num_submitted
        self._update_db()
        if self._total_jobs:
            if self._submitted_jobs > self._total_jobs:
                raise ValueError("Number of submitted jobs exceed number of total jobs")

    @property
    def total_jobs(self):
        """
        Get number of total jobs

        Returns:
            int: number of total jobs
        """
        self.refresh()
        return self._total_jobs

    @total_jobs.setter
    def total_jobs(self, num_total):
        """
        Set number of total jobs (optional: default = None)

        Args:
            num_total (int): number of submitted jobs
        """
        self._total_jobs = num_total
        self._update_db()

    @property
    def finished(self):
        """
        Check if the status is 'finished', meaning the job and all connected sub jobs are finished.

        Returns:
            bool: [True/False]
        """
        return self.string == "finished"

    @property
    def string(self):
        """
        Get the current status as string, it can be:
            initialized: No jobs have been submitted.
            sub_m_n: m out of n jobs have been submitted
            finished: The job and all connected sub jobs are finished.

        Returns:
            str: status [initialized, finished]
        """
        total_jobs = self.total_jobs
        submitted_jobs = self.submitted_jobs
        if submitted_jobs == 0:
            return "initialized"
        elif total_jobs and submitted_jobs == total_jobs:
            return "finished"
        elif total_jobs:
            return "submitted_{}_{}".format(submitted_jobs, total_jobs)
        else:
            return "submitted_{}".format(submitted_jobs)

    def submit_next(self):
        """
        Increasing the number of submitted jobs by one: self.submitted_jobs += 1
        """
        self.submitted_jobs += 1

    def refresh(self):
        """
        Refresh the submission status, if a job_id is present load the current submission status from the database.
        """
        if self.job_id:
            computer = self.database.get_item_by_id(self.job_id)["computer"]
            if computer is not None:
                submission_status_lst = (
                    computer
                    .split("#")[-1]
                    .split("/")
                )
            else:
                submission_status_lst = []
            if len(submission_status_lst) == 2:
                self._submitted_jobs = int(submission_status_lst[0])
                self._total_jobs = int(submission_status_lst[1])
            elif len(submission_status_lst) == 0:
                self._submitted_jobs = 0
                self._total_jobs = None
            else:
                self._submitted_jobs = int(submission_status_lst[0])
                self._total_jobs = None

    def __repr__(self):
        """
        Human readable representation of the submission status

        Returns:
            str: human readable string
        """
        return repr(self.string)

    def __str__(self):
        """
        Machine readable representation of the submission status

        Returns:
            str: machine readable string
        """
        return str(self.string)

    def _update_db(self):
        """
        Internal function to update the database, with the current number of submitted jobs.
        """
        if self.job_id:
            db_entry = self.database.get_item_by_id(self.job_id)
            if db_entry["computer"] is not None:
                split_str = db_entry["computer"].split("#")
                if len(split_str) > 2:
                    computer = split_str[:-1]
                else:
                    computer = split_str
                if self._total_jobs:
                    status = (
                        computer[0]
                        + "#"
                        + computer[1]
                        + "#"
                        + str(self._submitted_jobs)
                        + "/"
                        + str(self._total_jobs)
                    )
                else:
                    status = (
                        computer[0] + "#" + computer[1] + "#" + str(self._submitted_jobs)
                    )
                self.database.item_update({"computer": status}, self.job_id)
