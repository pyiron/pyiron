# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import six
from pyiron.base.database.generic import DatabaseAccess
from pyiron.base.database.filetable import FileTable

"""
The JobStatus class belongs to the GenericJob object.
"""

__author__ = "Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


job_status_lst = [
    "initialized",
    "appended",
    "created",
    "submitted",
    "running",
    "aborted",
    "collect",
    "suspended",
    "refresh",
    "busy",
    "finished",
    "not_converged",
    "warning"
]


def format_docstring_with_statuses(n_tabs=1):
    """
    Replaces a '{}' in the decorated object's docstring with the documentation for all possible job status.

    The intent here is to have a *single place* where all status-value documentation can be maintained.

    Args:
        n_tabs (int): The number of tabs preceding each status type. (Default is one)
    """

    status_docs = """
    {0}initialized: The object for the corresponding job was just created.
    {0}appended: The job was appended to an master job.
    {0}created: The files required for the simulation were written to the harddisk.
    {0}submitted: The job was submitted to the jobscheduler and is waiting to be executed.
    {0}running: The job is currently executed.
    {0}aborted: The job failed to execute.
    {0}collect: The job finished successfully and the written files are being collected.
    {0}suspended: The job was set to sleep, waiting until other related jobs are finished, before it continous.
    {0}refresh: The job was suspended before and it is currently checking if there are new tasks it can execute.
    {0}busy: The job is refreshing, but during the refresh more related jobs finished so another refresh is necessary.
    {0}finished: The job and all connected sub jobs are finished.
    {0}warning: The job finished but encountered warnings in its output.
    """

    def decorator(obj):
        obj.__doc__ = obj.__doc__.format(status_docs.format('\t' * n_tabs))
        return obj
    return decorator


@format_docstring_with_statuses()
class JobStatus(object):
    """
    The JobStatus object handles the different states a job could have. The available states are: {}

    Args:
        initial_status (str): If no initial status is provided the status is set to 'initialized'
        db (DatabaseAccess): The database which is responsible for this job.
        job_id (int): job ID

    Attributes:

        .. attribute:: database

            the database which is responsible for this job.

        .. attribute:: job_id

            Job ID

        .. attribute:: string

            job status as string
    """

    def __init__(self, initial_status="initialized", db=None, job_id=None):
        super(JobStatus, self).__setattr__("_status_dict", {})
        self._db = None
        self._job_id = None
        self.string = initial_status
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
        if db is not None and not isinstance(db, (DatabaseAccess, FileTable)):
            raise TypeError("The database has to be an DatabaseAccess object.")
        self._db = db

    @property
    def job_id(self):
        """
        Get the job id of the job this jobstatus is associated to.
        Returns:
            int: job id
        """
        return self._job_id

    @job_id.setter
    def job_id(self, unique_id):
        """
        Get the job id of the job this jobstatus is associated to.
        Args:
            unique_id (int): job id
        """
        if unique_id and not isinstance(unique_id, int):
            raise TypeError("The Job_ID should be an integer.")
        self._job_id = unique_id
        self.refresh_status()

    @format_docstring_with_statuses(n_tabs=2)
    @property
    def string(self):
        """
        Get the current status as string, it can be: {}

        Returns:
            (str): status [initialized, appended, created, submitted, running, aborted, collect, suspended, refresh,
                   busy, finished, warning]
        """
        self.refresh_status()
        return self._get_status_from_dict()

    @format_docstring_with_statuses(n_tabs=2)
    @string.setter
    def string(self, status):
        """
        Set the current status, to one of the following: {}

        Args:
            status (str): status [initialized, appended, created, submitted, running, aborted, collect, suspended,
                          refresh, busy, finished, warning]
        """
        self._reset()
        if isinstance(status, six.string_types) and status in self._status_dict.keys():
            self._status_dict[status] = True
            self._status_write()
        else:
            raise (
                "No valid job status: ",
                status,
                " Instead use [initialized, appended, created, submitted, running,"
                "aborted, collect, suspended, refresh, busy, finished, not_converged].",
            )

    def refresh_status(self):
        """
        Refresh the job status - check if the database and job_id are set and if this is the case load the job status
        from the database.
        """
        if self.database and self.job_id and not any([self._status_dict[i] for i in ["finished", "aborted", "warning", "not_converged"]]):
            try:
                status = self.database.get_job_status(job_id=self.job_id)
            except IndexError:
                raise (
                    "The job with the job ID "
                    + str(self.job_id)
                    + " is not listed in the database anymore."
                )
            self._reset()
            self._status_dict[status] = True

    def _status_write(self):
        """
        Private function: Write the job status to the internal variable _key and store it in the database.
        """
        if self.database and self.job_id:
            current_status = str(self._get_status_from_dict())
            if self.database.get_item_by_id(self.job_id)["status"] != current_status:
                self.database.set_job_status(job_id=self.job_id, status=current_status)

    def _reset(self):
        """
        internal function to reset the run mode - sets all run modes to false.
        """
        self._status_dict = {status: False for status in job_status_lst}

    @staticmethod
    def _bool_check(boolean):
        """
        Private function: Raise TypeError if boolean is not type bool and raise a ValueError if it is not True.
        Args:
            boolean (bool): True
        """
        if not isinstance(boolean, bool):
            raise TypeError(
                "The JobStatus can only be set to a boolean, more specifically to True only."
            )
        if boolean is False:
            raise ValueError("The JobStatus can only be set to True.")

    def _get_status_from_dict(self):
        return [key for key, val in self._status_dict.items() if val][0]

    def __repr__(self):
        """
        Human readable representation of the job status
        Returns:
            str: human readable string
        """
        return repr(self.string)

    def __str__(self):
        """
        Machine readable representation of the job status
        Returns:
            str: machine readable string
        """
        return str(self.string)

    def __getattr__(self, name):
        if name in self._status_dict.keys():
            self.refresh_status()
            return self._status_dict[name]
        raise AttributeError("'{}' object has no attribute '{}'".format(
                                self.__class__.__name__, name))

    def __setattr__(self, name, value):
        if name in self._status_dict.keys():
            if not isinstance(value, bool):
                raise TypeError("A run mode can only be activated using [True].")
            if value:
                self.string = name
            else:
                raise ValueError("A run mode can only be activated using [True].")
        else:
            super(JobStatus, self).__setattr__(name, value)

    def __dir__(self):
        return list(self._status_dict.keys())

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return other._status_dict == self._status_dict
        elif isinstance(other, str):
            return other == self.string
        else:
            return super(JobStatus, self).__eq__(other)

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return other._status_dict != self._status_dict
        elif isinstance(other, str):
            return other != self.string
        else:
            return super(JobStatus, self).__ne__(other)
