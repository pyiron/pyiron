# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron_base.core.settings.database import DatabaseAccess

"""
The JobStatus class belongs to the GenericJob object.
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class JobStatus(object):
    """
    The JobStatus object handles the different states a job could have. The available states are:
        initialized: The object for the corresponding job was just created.
        appended: The job was appended to an master job.
        created: The files required for the simulation were written to the harddisk.
        submitted: The job was submitted to the jobscheduler and is waiting to be executed.
        running: The job is currently executed.
        aborted: The job failed to execute.
        collect: The job finished successfully and the written files are being collected.
        suspended: The job was set to sleep, waiting until other related jobs are finished, before it continous.
        refresh: The job was suspended before and it is currently checking if there are new tasks it can execute.
        busy: The job is refreshing, but during the refresh more related jobs finished so another refresh is necessary.
        finished: The job and all connected sub jobs are finished.

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

    def __init__(self, initial_status='initialized', db=None, job_id=None):
        self._key = None
        self._db = None
        self._job_id = None
        self.string = initial_status
        self.database = db
        self.job_id = job_id

        self.INITIALIZED = 'initialized'
        self.APPEND = 'append'
        self.CREATED = 'created'
        self.SUBMITTED = 'submitted'
        self.RUNNING = 'running'
        self.ABORTED = 'aborted'
        self.COLLECT = 'collect'
        self.SUSPEND = 'suspend'
        self.REFRESH = 'refresh'
        self.BUSY = 'busy'
        self.FINISHED = 'finished'

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
        if db and not isinstance(db, DatabaseAccess):
            raise TypeError('The database has to be an DatabaseAccess object.')
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
            raise TypeError('The Job_ID should be an integer.')
        self._job_id = unique_id
        self.refresh_status()

    @property
    def initialized(self):
        """
        Check if the status is 'initialized', meaning the object for the corresponding job was just created.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='initialized')

    @initialized.setter
    def initialized(self, boolean):
        """
        Set the status to 'initialized', meaning the object for the corresponding job was just created.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='initialized')

    @property
    def appended(self):
        """
        Check if the status is 'appended', meaning the job was appended to an master job.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='appended')

    @appended.setter
    def appended(self, boolean):
        """
        Set the status to 'appended', meaning the job was appended to an master job.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='appended')

    @property
    def created(self):
        """
        Check if the status is 'created', meaning the files required for the simulation were written to the harddisk.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='created')

    @created.setter
    def created(self, boolean):
        """
        Set the status to 'created', meaning the files required for the simulation were written to the harddisk.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='created')

    @property
    def submitted(self):
        """
        Check if the status is 'submitted', meaning the job was submitted to the jobscheduler and is waiting to be
        executed.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='submitted')

    @submitted.setter
    def submitted(self, boolean):
        """
        Set the status to 'submitted', meaning the job was submitted to the jobscheduler and is waiting to be executed.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='submitted')

    @property
    def running(self):
        """
        Check if the status is 'running', meaning the job is currently executed.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='running')

    @running.setter
    def running(self, boolean):
        """
        Set the status to 'running', meaning the job is currently executed.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='running')

    @property
    def aborted(self):
        """
        Check if the status is 'aborted', meaning the job failed to execute.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='aborted')

    @aborted.setter
    def aborted(self, boolean):
        """
        Set the status to 'aborted', meaning the job failed to execute.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='aborted')

    @property
    def collect(self):
        """
        Check if the status is 'collect', meaning the job finished successfully and the written files are being
        collected.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='collect')

    @collect.setter
    def collect(self, boolean):
        """
        Set the status to 'collect', meaning the job finished successfully and the written files are being collected.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='collect')

    @property
    def suspended(self):
        """
        Check if the status is 'suspended', meaning the job was set to sleep, waiting until other related jobs are
        finished, before it continous.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='suspended')

    @suspended.setter
    def suspended(self, boolean):
        """
        Set the status to 'suspended', meaning the job was set to sleep, waiting until other related jobs are finished,
        before it continous.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='suspended')

    @property
    def refresh(self):
        """
        Check if the status is 'refresh', meaning the job was suspended before and it is currently checking if there are
        new tasks it can execute.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='refresh')

    @refresh.setter
    def refresh(self, boolean):
        """
        Set the status to 'refresh', meaning the job was suspended before and it is currently checking if there are new
        tasks it can execute.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='refresh')

    @property
    def busy(self):
        """
        Check if the status is 'busy', meaning the job is refreshing, but during the refresh more related jobs finished
        so another refresh is necessary.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='busy')

    @busy.setter
    def busy(self, boolean):
        """
        Set the status to 'busy', meaning The job is refreshing, but during the refresh more related jobs finished so
        another refresh is necessary.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='busy')

    @property
    def finished(self):
        """
        Check if the status is 'finished', meaning the job and all connected sub jobs are finished.

        Returns:
            (bool): [True/False]
        """
        return self._status_validate(status='finished')

    @finished.setter
    def finished(self, boolean):
        """
        Set the status to 'finished', meaning the job and all connected sub jobs are finished.

        Args:
            boolean (bool): [True]
        """
        self._bool_check(boolean)
        self._status_write(status='finished')

    @property
    def string(self):
        """
        Get the current status as string, it can be:
            initialized: The object for the corresponding job was just created.
            appended: The job was appended to an master job.
            created: The files required for the simulation were written to the harddisk.
            submitted: The job was submitted to the jobscheduler and is waiting to be executed.
            running: The job is currently executed.
            aborted: The job failed to execute.
            collect: The job finished successfully and the written files are being collected.
            suspended: The job was set to sleep, waiting until other related jobs are finished, before it continous.
            refresh: The job was suspended before and it is currently checking if there are new tasks it can execute.
            busy: The job is refreshing, but during the refresh more related jobs finished so another refresh is 
                  necessary.
            finished: The job and all connected sub jobs are finished.

        Returns:
            (str): status [initialized, appended, created, submitted, running, aborted, collect, suspended, refresh,
                   busy, finished]
        """
        if not self._key:
            return 'initialized'
        return self._key

    @string.setter
    def string(self, status):
        """
        Set the current status, to one of the following:
            initialized: The object for the corresponding job was just created.
            appended: The job was appended to an master job.
            created: The files required for the simulation were written to the harddisk.
            submitted: The job was submitted to the jobscheduler and is waiting to be executed.
            running: The job is currently executed.
            aborted: The job failed to execute.
            collect: The job finished successfully and the written files are being collected.
            suspended: The job was set to sleep, waiting until other related jobs are finished, before it continous.
            refresh: The job was suspended before and it is currently checking if there are new tasks it can execute.
            busy: The job is refreshing, but during the refresh more related jobs finished so another refresh is 
                  necessary.
            finished: The job and all connected sub jobs are finished.

        Args:
            status (str): status [initialized, appended, created, submitted, running, aborted, collect, suspended,
                          refresh, busy, finished]
        """
        if status not in ['initialized', 'appended', 'created', 'submitted', 'running', 'aborted', 'collect',
                          'suspended', 'refresh', 'busy', 'finished']:
            raise ('No valid job status: ', status, ' Instead use [initialized, appended, created, submitted, running,'
                   'aborted, collect, suspended, refresh, busy, finished].')
        if status == 'initialized':
            self.initialized = True
        elif status == 'appended':
            self.appended = True
        elif status == 'created':
            self.created = True
        elif status == 'submitted':
            self.submitted = True
        elif status == 'running':
            self.running = True
        elif status == 'aborted':
            self.aborted = True
        elif status == 'collect':
            self.collect = True
        elif status == 'suspended':
            self.suspended = True
        elif status == 'refresh':
            self.refresh = True
        elif status == 'busy':
            self.busy = True
        elif status == 'finished':
            self.finished = True

    def refresh_status(self):
        """
        Refresh the job status - check if the database and job_id are set and if this is the case load the job status
        from the database.
        """
        if self.database and self.job_id:
            try:
                self.string = self.database.get_item_by_id(self.job_id)["status"]
            except IndexError:
                raise('The job with the job ID ' + str(self.job_id) + ' is not listed in the database anymore.')

    def _status_validate(self, status):
        """
        Private function: Get the current job status from the database and check if it is the same like the one provided
        in the status variable.

        Args:
            status (str): [initialized, appended, created, submitted, running, aborted, collect, suspended, refresh,
                          busy, finished]

        Returns:
            bool: [True/False]
        """
        self.refresh_status()
        if self._key == status:
            return True
        return False

    def _status_write(self, status):
        """
        Private function: Write the job status to the internal variable _key and store it in the database.

        Args:
            status (str): [initialized, appended, created, submitted, running, aborted, collect, suspended, refresh,
                          busy, finished]
        """
        self._key = status
        if self.database and self.job_id:
            if self.database.get_item_by_id(self.job_id)["status"] != str(self.string):
                self.database.item_update({'status': str(self.string)}, self.job_id)

    @staticmethod
    def _bool_check(boolean):
        """
        Private function: Raise TypeError if boolean is not type bool and raise a ValueError if it is not True.

        Args:
            boolean (bool): True
        """
        if not isinstance(boolean, bool):
            raise TypeError('The JobStatus can only be set to a boolean, more specifically to True only.')
        if boolean is False:
            raise ValueError('The JobStatus can only be set to True.')

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
