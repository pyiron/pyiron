# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import subprocess
import pandas
import time
from pyiron_base.core.settings.generic import Settings
from pyiron_base.objects.job.generic import GenericJob

"""
Set of functions to interact with the queuing system directly from within pyiron - optimized for the Sun grid engine.
"""

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


def _queue_status_command(user):
    return ['qstat', '-u', str(user)]


def _queue_del_command(que_id):
    return ['qdel', str(que_id)]


def _queue_enable_reservation(que_id):
    return ['qalter', '-R', 'y', str(que_id)]


def _queue_job_status_command(que_id):
    return ['qstat', '-j', str(que_id)]


def _queue_job_details_command(user, que_id):
    return ['qacct', '-o', str(user), '-j', str(que_id)]


def _queue_function(funct, user=None, que_id=None):
    if user is not None and que_id is not None:
        return subprocess.check_output(funct(user=user, que_id=que_id), stderr=subprocess.STDOUT, 
                                       universal_newlines=True).split('\n')
    elif user is not None:
        return subprocess.check_output(funct(user=user), stderr=subprocess.STDOUT, 
                                       universal_newlines=True).split('\n')
    elif que_id is not None:
        return subprocess.check_output(funct(que_id=que_id), stderr=subprocess.STDOUT, 
                                       universal_newlines=True).split('\n')
    else: 
        return subprocess.check_output(funct(), stderr=subprocess.STDOUT, 
                                       universal_newlines=True).split('\n')

    
def queue_table(job_ids=[], project_only=True):
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
    output = _queue_function(funct=_queue_status_command, user=s.login_user)
    if not output == ['']:
        job_dict_lst = []
        job_title = output[0].split()
        for line in output[2:]:
            job = line.split()
            while not len(job_title) == len(job):
                job.append('')
            try:
                in_list = int(job[2].split('_')[1].split('.')[0]) in job_ids
            except (IndexError, ValueError):
                in_list = False
            if (project_only and in_list) or not project_only:
                job_dict_lst.append(dict(zip(job_title, job)))
        if project_only:
            return pandas.DataFrame(job_dict_lst)
        else:
            return pandas.DataFrame(job_dict_lst[:-1])
    else:
        return None


def queue_id_table(requested_id=None):
    """
    Display the queuing system table as dictionary

    Args:
        requested_id (int): query for a specific job_id - optional

    Returns:
        dict: Dictionary with the output from the queuing system - optimized for the Sun grid engine
    """
    output = _queue_function(funct=_queue_status_command, user=s.login_user)
    if not output == ['']:
        jobs_dict = {}
        for line in output[2:]:
            job = line.split()
            if len(job) > 4 and job[2].startswith('pi_'):
                job_id = job[2].split('.sh')[0].split('py_')[1]
                job_status = job[4]
                jobs_dict[job_id] = (int(job[0]), job_status)
                if requested_id is not None and job_id == requested_id:
                    return {job_id: (int(job[0]), job_status)}
        return jobs_dict
    else:
        return {}


def queue_is_empty():
    """
    Check if the queue table is currently empty - no more jobs to wait for.

    Returns:
        bool: True if the table is empty, else False - optimized for the Sun grid engine
    """
    output = _queue_function(funct=_queue_status_command, user=s.login_user)
    if output == ['']:
        return True
    else:
        return False


def queue_delete_job(item):
    """
    Delete a job from the queuing system

    Args:
        item (int, GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        str: Output from the queuing system as string - optimized for the Sun grid engine
    """
    que_id = _validate_que_request(item)
    return _queue_function(funct=_queue_del_command, que_id=que_id)[0]


def queue_enable_reservation(item):
    """
    Enable a reservation for a particular job within the queuing system

    Args:
        item (int, GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        str: Output from the queuing system as string - optimized for the Sun grid engine
    """
    que_id = _validate_que_request(item)
    return _queue_function(funct=_queue_enable_reservation, que_id=que_id)[0]


def queue_report(item):
    """
    Detailed reporting for a particular job - using the qacct command of the sun grid engine.

    Args:
        item (int, GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        pandas.DataFrame: Detailed report returned from the queuing system - optimized for the Sun grid engine
    """
    que_id = _validate_que_request(item)
    try:
        output = _queue_function(funct=_queue_job_details_command, user=s.login_user, que_id=que_id)
        return pandas.DataFrame(
            [[line.split()[0], line.split()[1:]] for line in output[1:]])
    except subprocess.CalledProcessError:
        pass


def queue_job_info(item):
    """
    Short reporting for a particular job - using the qstat command of the sun grid engine.

    Args:
        item (int, GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        pandas.DataFrame: Short report returned from the queuing system - optimized for the Sun grid engine
    """
    que_id = _validate_que_request(item)
    try:
        output = _queue_function(funct=_queue_job_status_command, que_id=que_id)
        return pandas.DataFrame(
            [[line_arg.replace('  ', '') for line_arg in line.split(': ')] for line in output[1:]][0:24])
    except subprocess.CalledProcessError:
        pass


def wait_for_job(job, interval_in_s=5, max_iterations=100):
    """
    Sleep until the job is finished but maximum interval_in_s * max_iterations seconds.

    Args:
        job (GenericJob): Job to wait for
        interval_in_s (int): interval when the job status is queried from the database - default 5 sec.
        max_iterations (int): maximum number of iterations - default 100
    """
    for _ in range(max_iterations):
        job.refresh_job_status()
        if job.status.finished or job.status.aborted:
            break
        time.sleep(interval_in_s)


def _validate_que_request(item):
    """
    Internal function to convert the job_ID or hamiltonian to the queuing system ID.

    Args:
        item (int, GenericJob): Provide either the job_ID or the full hamiltonian

    Returns:
        int: queuing system ID
    """
    if isinstance(item, int):
        que_id = item
    elif isinstance(item, GenericJob):
        if item.server.queue_id:
            que_id = item.server.queue_id
        else:
            raise ValueError('This job does not have a queue ID.')
    else:
        raise TypeError('The queue can either query for IDs or for pyiron GenericJobObjects.')
    return que_id
