# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import pandas
import sys
import os
import numpy as np
from pyiron.base.settings.generic import Settings

"""
The Jobtable module provides a set of top level functions to interact with the database.
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

s = Settings()


def _job_dict(
    database,
    sql_query,
    user,
    project_path,
    recursive,
    job=None,
    sub_job_name="%",
    element_lst=None,
):
    """
    Internal function to access the database from the project directly.

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        recursive (bool): search subprojects [True/False]
        job (str): job_name - by default None
        sub_job_name (str): path inside the HDF5 file - "%" by default to accept any path
        element_lst (list): list of elements required in the chemical formular - by default None

    Returns:
        list: the function returns a list of dicts like get_items_sql, but it does not format datetime:
             [{'chemicalformula': u'Ni108',
              'computer': u'mapc157',
              'hamilton': u'LAMMPS',
              'hamversion': u'1.1',
              'id': 24,
              'job': u'DOF_1_0',
              'parentid': 21L,
              'project': u'lammps.phonons.Ni_fcc',
              'projectpath': u'D:/PyIron/PyIron_data/projects',
              'status': u'finished',
              'timestart': datetime.datetime(2016, 6, 24, 10, 17, 3, 140000),
              'timestop': datetime.datetime(2016, 6, 24, 10, 17, 3, 173000),
              'totalcputime': 0.033,
              'username': u'test'},
             {'chemicalformula': u'Ni108',
              'computer': u'mapc157',
              'hamilton': u'LAMMPS',
              'hamversion': u'1.1',
              'id': 21,
              'job': u'ref',
              'parentid': 20L,
              'project': u'lammps.phonons.Ni_fcc',
              'projectpath': u'D:/PyIron/PyIron_data/projects',
              'status': u'finished',
              'timestart': datetime.datetime(2016, 6, 24, 10, 17, 2, 429000),
              'timestop': datetime.datetime(2016, 6, 24, 10, 17, 2, 463000),
              'totalcputime': 0.034,
              'username': u'test'},.......]

    """
    dict_clause = {}
    # FOR GET_ITEMS_SQL: clause = []
    if user is not None:
        dict_clause["username"] = str(user)
        # FOR GET_ITEMS_SQL: clause.append("username = '" + self.user + "'")
    if sql_query is not None:
        # FOR GET_ITEMS_SQL: clause.append(self.sql_query)
        if "AND" in sql_query:
            cl_split = sql_query.split(" AND ")
        elif "and" in sql_query:
            cl_split = sql_query.split(" and ")
        else:
            cl_split = [sql_query]
        dict_clause.update(
            {str(element.split()[0]): element.split()[2] for element in cl_split}
        )
    if job is not None:
        dict_clause["job"] = str(job)

    if recursive:
        dict_clause["project"] = str(project_path) + "%"
    else:
        dict_clause["project"] = str(project_path)
    if sub_job_name is None:
        dict_clause["subjob"] = None
    elif sub_job_name != "%":
        dict_clause["subjob"] = str(sub_job_name)
    if element_lst is not None:
        dict_clause["element_lst"] = element_lst

    s.logger.debug("sql_query: %s", str(dict_clause))
    return database.get_items_dict(dict_clause)


def get_db_columns(database):
    """
    Get column names

    Args:
        database (DatabaseAccess): Database object

    Returns:
        list: list of column names like:
             ['id',
             'parentid',
             'masterid',
             'projectpath',
             'project',
             'job',
             'subjob',
             'chemicalformula',
             'status',
             'hamilton',
             'hamversion',
             'username',
             'computer',
             'timestart',
             'timestop',
             'totalcputime']
    """
    return database.get_table_headings()


def job_table(
    database,
    sql_query,
    user,
    project_path,
    recursive=True,
    columns=None,
    all_columns=False,
    sort_by="id",
    max_colwidth=200,
    full_table=False,
    element_lst=None,
    job_name_contains='',
):
    """
    Access the job_table

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        recursive (bool): search subprojects [True/False]
        columns (list): by default only the columns ['job', 'project', 'chemicalformula'] are selected, but the
                        user can select a subset of ['id', 'status', 'chemicalformula', 'job', 'subjob', 'project',
                        'projectpath', 'timestart', 'timestop', 'totalcputime', 'computer', 'hamilton', 'hamversion',
                        'parentid', 'masterid']
        all_columns (bool): Select all columns - this overwrites the columns option.
        sort_by (str): Sort by a specific column
        max_colwidth (int): set the column width
        full_table (bool): Whether to show the entire pandas table
        element_lst (list): list of elements required in the chemical formular - by default None
        job_name_contains (str): a string which should be contained in every job_name

    Returns:
        pandas.Dataframe: Return the result as a pandas.Dataframe object
    """
    if columns is None:
        columns = ["job", "project", "chemicalformula"]
    all_db = [
        "id",
        "status",
        "chemicalformula",
        "job",
        "subjob",
        "projectpath",
        "project",
        "timestart",
        "timestop",
        "totalcputime",
        "computer",
        "hamilton",
        "hamversion",
        "parentid",
        "masterid",
    ]

    if all_columns:
        columns = all_db
    job_dict = _job_dict(
        database=database,
        sql_query=sql_query,
        user=user,
        project_path=project_path,
        recursive=recursive,
        element_lst=element_lst,
    )
    if full_table:
        pandas.set_option('display.max_rows', None)
        pandas.set_option('display.max_columns', None)
    pandas.set_option("display.max_colwidth", max_colwidth)
    df = pandas.DataFrame(job_dict)
    if len(job_dict) == 0:
        return df
    if job_name_contains != '':
        df = df[df.job.str.contains(job_name_contains)]
    if sort_by in columns:
        return df[columns].sort_values(by=sort_by)
    return df[columns]


def get_jobs(database, sql_query, user, project_path, recursive=True, columns=None):
    """
    Internal function to return the jobs as dictionary rather than a pandas.Dataframe

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        recursive (bool): search subprojects [True/False]
        columns (list): by default only the columns ['id', 'project'] are selected, but the user can select a subset
                        of ['id', 'status', 'chemicalformula', 'job', 'subjob', 'project', 'projectpath', 'timestart',
                        'timestop', 'totalcputime', 'computer', 'hamilton', 'hamversion', 'parentid', 'masterid']

    Returns:
        dict: columns are used as keys and point to a list of the corresponding values
    """
    if columns is None:
        columns = ["id", "project"]
    df = job_table(database, sql_query, user, project_path, recursive, columns=columns)
    if len(df) == 0:
        dictionary = {}
        for key in columns:
            dictionary[key] = list()
        return dictionary
        # return {key: list() for key in columns}
    dictionary = {}
    for key in df.keys():
        dictionary[key] = df[
            key
        ].tolist()  # ToDo: Check difference of tolist and to_list
    return dictionary
    # return {key: df[key].tolist() for key in df.keys()}


def get_job_ids(database, sql_query, user, project_path, recursive=True):
    """
    Return the job IDs matching a specific query

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        recursive (bool): search subprojects [True/False]

    Returns:
        list: a list of job IDs
    """
    return get_jobs(database, sql_query, user, project_path, recursive=recursive)["id"]


def get_child_ids(database, sql_query, user, project_path, job_specifier, status=None):
    """
    Get the childs for a specific job

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        job_specifier (str): name of the master job or the master jobs job ID
        status (str): filter childs which match a specific status - None by default

    Returns:
        list: list of child IDs
    """
    id_master = get_job_id(database, sql_query, user, project_path, job_specifier)
    if id_master is None:
        return []
    else:
        search_dict = {"masterid": str(id_master)}
        if status is not None:
            search_dict["status"] = status
        return sorted(
            [
                job["id"]
                for job in database.get_items_dict(
                    search_dict, return_all_columns=False
                )
            ]
        )


def get_job_id(database, sql_query, user, project_path, job_specifier):
    """
    get the job_id for job named job_name in the local project path from database

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        job_specifier (str): name of the job or job ID

    Returns:
        int: job ID of the job
    """
    if isinstance(job_specifier, (int, np.integer)):
        return job_specifier  # is id

    job_specifier.replace(".", "_")
    # if job_specifier[0] is not '/':
    #     sub_job_name = '/' + job_specifier
    # else:
    #     sub_job_name = job_specifier
    # job_dict = _job_dict(database, sql_query, user, project_path, recursive=False,  # job=job_specifier,
    #                      sub_job_name=sub_job_name)
    # if len(job_dict) == 0:
    #     job_dict = _job_dict(database, sql_query, user, project_path, recursive=True,  # job=job_specifier,
    #                          sub_job_name=sub_job_name)
    job_dict = _job_dict(
        database, sql_query, user, project_path, recursive=False, job=job_specifier
    )
    if len(job_dict) == 0:
        job_dict = _job_dict(
            database, sql_query, user, project_path, recursive=True, job=job_specifier
        )
    if len(job_dict) == 0:
        return None
    elif len(job_dict) == 1:
        return job_dict[0]["id"]
    else:
        raise ValueError(
            "job name '{0}' in this project '{1}' is not unique '{2}".format(job_specifier, project_path, job_dict)
        )


def set_job_status(database, sql_query, user, project_path, job_specifier, status):
    """
    Set the status of a particular job

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        job_specifier (str): name of the job or job ID
        status (str): job status can be one of the following ['initialized', 'appended', 'created', 'submitted',
                     'running', 'aborted', 'collect', 'suspended', 'refresh', 'busy', 'finished']

    """
    database.item_update(
        {"status": str(status)},
        get_job_id(database, sql_query, user, project_path, job_specifier),
    )


def get_job_status(database, sql_query, user, project_path, job_specifier):
    """
    Get the status of a particular job

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        job_specifier (str): name of the job or job ID

    Returns:
        str: job status can be one of the following ['initialized', 'appended', 'created', 'submitted', 'running',
             'aborted', 'collect', 'suspended', 'refresh', 'busy', 'finished']
    """
    try:
        return database.get_item_by_id(
            get_job_id(database, sql_query, user, project_path, job_specifier)
        )["status"]
    except KeyError:
        return None


def get_job_working_directory(database, sql_query, user, project_path, job_specifier):
    """
    Get the working directory of a particular job

    Args:
        database (DatabaseAccess): Database object
        sql_query (str): SQL query to enter a more specific request
        user (str): username of the user whoes user space should be searched
        project_path (str): root_path - this is in contrast to the project_path in GenericPath
        job_specifier (str): name of the job or job ID

    Returns:
        str: working directory as absolute path
    """
    try:
        db_entry = database.get_item_by_id(
            get_job_id(database, sql_query, user, project_path, job_specifier)
        )
        if db_entry:
            job_name = db_entry["subjob"][1:]
            return os.path.join(
                db_entry["projectpath"],
                db_entry["project"],
                job_name + "_hdf5",
                job_name,
            )
        else:
            return None
    except KeyError:
        return None
