# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import pandas
import time
from pyiron.base.settings.generic import Settings
from pyiron.base.job.jobtype import static_isinstance

"""
Set of functions to interact with the queuing system directly from within pyiron - optimized for the Sun grid engine.
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

QUEUE_SCRIPT_PREFIX = "pi_"

s = Settings()


def queue_table(job_ids=[], project_only=True, full_table=False):
    """
    Display the queuing system table as pandas.Dataframe

    Args:
        job_ids (list): check for a specific list of job IDs - empty list by default
        project_only (bool): Query only for jobs within the current project - True by default

    Returns:
        pandas.DataFrame: Output from the queuing system - optimized for the Sun grid engine
    """
    if project_only and not job_ids:
        return []
    if s.queue_adapter is not None:
        if full_table:
            pandas.set_option('display.max_rows', None)
            pandas.set_option('display.max_columns', None)
        else:
            pandas.reset_option('display.max_rows')
            pandas.reset_option('display.max_columns')
        df = s.queue_adapter.get_status_of_my_jobs()
        if not project_only:
            return df[
                [
                    True if QUEUE_SCRIPT_PREFIX in job_name else False
                    for job_name in list(df.jobname)
                ]
            ]
        else:
            job_name_lst = [QUEUE_SCRIPT_PREFIX + str(job_id) for job_id in job_ids]
            return df[
                [
                    True if job_name in job_name_lst else False
                    for job_name in list(df.jobname)
                ]
            ]
    else:
        return None


def queue_check_job_is_waiting_or_running(item):
    """
    Check if a job is still listed in the queue system as either waiting or running.

    Args:
        item (int, GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        bool: [True/False]
    """
    que_id = _validate_que_request(item)
    if s.queue_adapter is not None:
        return s.queue_adapter.get_status_of_job(process_id=que_id) in [
            "pending",
            "running",
        ]
    else:
        return None


def queue_info_by_job_id(job_id):
    """
    Display the queuing system info of job by qstat | grep  shell command
    as dictionary

    Args:
        requested_id (int): query for a specific job_id

    Returns:
        dict: Dictionary with the output from the queuing system - optimized for the Sun grid engine
    """
    if s.queue_adapter is not None:
        return s.queue_adapter.get_status_of_job(process_id=job_id)
    else:
        return None


def queue_is_empty():
    """
    Check if the queue table is currently empty - no more jobs to wait for.

    Returns:
        bool: True if the table is empty, else False - optimized for the Sun grid engine
    """
    if s.queue_adapter is not None:
        return len(s.queue_adapter.get_status_of_my_jobs()) == 0
    else:
        return True


def queue_delete_job(item):
    """
    Delete a job from the queuing system

    Args:
        item (int, pyiron.base.job.generic.GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        str: Output from the queuing system as string - optimized for the Sun grid engine
    """
    que_id = _validate_que_request(item)
    if s.queue_adapter is not None:
        return s.queue_adapter.delete_job(process_id=que_id)
    else:
        return None


def queue_enable_reservation(item):
    """
    Enable a reservation for a particular job within the queuing system

    Args:
        item (int, pyiron.base.job.generic.GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        str: Output from the queuing system as string - optimized for the Sun grid engine
    """
    que_id = _validate_que_request(item)
    if s.queue_adapter is not None:
        if isinstance(que_id, list):
            return [s.queue_adapter.enable_reservation(process_id=q) for q in que_id]
        else:
            return s.queue_adapter.enable_reservation(process_id=que_id)
    else:
        return None


def wait_for_job(job, interval_in_s=5, max_iterations=100):
    """
    Sleep until the job is finished but maximum interval_in_s * max_iterations seconds.

    Args:
        job (pyiron.base.job.generic.GenericJob): Job to wait for
        interval_in_s (int): interval when the job status is queried from the database - default 5 sec.
        max_iterations (int): maximum number of iterations - default 100
    """
    if s.queue_adapter is not None and s.queue_adapter.remote_flag and job.server.queue is not None:
        finished = False
        for _ in range(max_iterations):
            if not job.project.queue_check_job_is_waiting_or_running(job):
                job.transfer_from_remote()
                finished = True
                break
            time.sleep(interval_in_s)
        if not finished:
            raise ValueError("Maximum iterations reached, but the job was not finished.")
    else:
        finished = False
        for _ in range(max_iterations):
            job.refresh_job_status()
            if job.status.finished or job.status.aborted or job.status.not_converged:
                finished = True
                break
            time.sleep(interval_in_s)
        if not finished:
            raise ValueError("Maximum iterations reached, but the job was not finished.")


def _validate_que_request(item):
    """
    Internal function to convert the job_ID or hamiltonian to the queuing system ID.

    Args:
        item (int, pyiron.base.job.generic.GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        int: queuing system ID
    """

    if isinstance(item, int):
        que_id = item
    elif static_isinstance(item.__class__, "pyiron.base.master.generic.GenericMaster"):
        if item.server.queue_id:
            que_id = item.server.queue_id
        else:
            queue_id_lst = [item.project.load(child_id).server.queue_id for child_id in item.child_ids]
            que_id = [queue_id for queue_id in queue_id_lst if queue_id is not None]
            if len(que_id) == 0:
                raise ValueError("This job does not have a queue ID.")
    elif static_isinstance(item.__class__, "pyiron.base.job.generic.GenericJob"):
        if item.server.queue_id:
            que_id = item.server.queue_id
        else:
            raise ValueError("This job does not have a queue ID.")
    elif static_isinstance(item.__class__, "pyiron.base.job.core.JobCore"):
        if "server" in item.project_hdf5.list_nodes():
            server_hdf_dict = item.project_hdf5["server"]
            if "qid" in server_hdf_dict.keys():
                que_id = server_hdf_dict["qid"]
            else:
                raise ValueError("This job does not have a queue ID.")
        else:
            raise ValueError("This job does not have a queue ID.")
    else:
        raise TypeError(
            "The queue can either query for IDs or for pyiron GenericJobObjects."
        )
    return que_id
