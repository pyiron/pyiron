import numpy as np
import os
import pandas
import datetime
import h5io
import sys
from pyfileindex import PyFileIndex


def filter_function(file_name):
    return '.h5' in file_name


class FileTable(object):
    def __init__(self, project):
        self._fileindex = PyFileIndex(path=project, filter_function=filter_function)
        df = pandas.DataFrame(self.init_table(fileindex=self._fileindex.dataframe))
        self._project = os.path.abspath(project)
        self._columns = ['id', 'status', 'chemicalformula', 'job', 'subjob', 'projectpath', 'project', 'timestart',
                         'timestop', 'totalcputime', 'computer', 'hamilton', 'hamversion', 'parentid', 'masterid']
        self._job_table = df[self._columns]

    def init_table(self, fileindex, working_dir_lst=None):
        if working_dir_lst is None:
            working_dir_lst = []
        fileindex = fileindex[~fileindex.is_directory]
        fileindex = fileindex.iloc[fileindex.path.values.argsort()]
        job_lst = []
        for path, mtime in zip(fileindex.path, fileindex.mtime):
            job_dict = self.get_extract(path, mtime)
            job_dict['id'] = len(working_dir_lst) + 1
            working_dir_lst.append(job_dict['project'] + job_dict['subjob'] + '_hdf5')
            if job_dict['project'] in working_dir_lst:
                job_dict['masterid'] = working_dir_lst.index(job_dict['project']) + 1
            else:
                job_dict['masterid'] = None
            job_lst.append(job_dict)
        return job_lst

    @staticmethod
    def get_extract(path, mtime):
        basename = os.path.basename(path)
        job = os.path.splitext(basename)[0]
        time = datetime.datetime.fromtimestamp(mtime)
        return {'status': h5io.read_hdf5(path, job + '/status'),
                'chemicalformula': None,
                'job': job,
                'subjob': '/' + job,
                'projectpath': None,
                'project': os.path.dirname(path),
                'timestart': time,
                'timestop': time,
                'totalcputime': 0.0,
                'computer': None,
                'parentid': None,
                'hamilton': h5io.read_hdf5(path, job + '/TYPE').split(".")[-1].split("'")[0],
                'hamversion': h5io.read_hdf5(path, job + '/VERSION')}

    def update(self):
        self._fileindex.update()
        files_lst, working_dir_lst = zip(*[[project + subjob + '.h5', project + subjob + '_hdf5']
                                           for project, subjob in zip(self._job_table.project.values,
                                                                      self._job_table.subjob.values)])
        df_new = self._fileindex.dataframe[
            ~self._fileindex.dataframe.is_directory & ~self._fileindex.dataframe.path.isin(files_lst)]
        if len(df_new) > 0:
            job_lst = self.init_table(fileindex=df_new, working_dir_lst=list(working_dir_lst))
            df = pandas.DataFrame(job_lst)[self._columns]
            self._job_table = pandas.concat([self._job_table, df])

    def get_db_columns(self):
        return self._job_table.columns.values

    def job_table(self, project=None, recursive=True, columns=None, all_columns=False, sort_by="id", max_colwidth=200,
                  job_name_contains=''):
        if project is None:
            project = self._project
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
        if recursive:
            df = self._job_table[self._job_table.project.str.contains(project)]
        else:
            df = self._job_table[self._job_table.project == project]
        pandas.set_option("display.max_colwidth", max_colwidth)
        if len(df) == 0:
            return df
        if job_name_contains != '':
            df = df[df.job.str.contains(job_name_contains)]
        if sort_by in columns:
            return df[columns].sort_values(by=sort_by)
        return df[columns]

    def get_jobs(self, project=None, recursive=True, columns=None):
        if project is None:
            project = self._project
        if columns is None:
            columns = ["id", "project"]
        df = self.job_table(project=project, recursive=recursive, columns=columns)
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

    def get_job_ids(self, project=None, recursive=True):
        if project is None:
            project = self._project
        return self.get_jobs(project=project, recursive=recursive, columns=['id'])["id"]

    def get_job_id(self, job_specifier, project=None):
        if project is None:
            project = self._project
        if sys.version_info.major == 2:
            if isinstance(job_specifier, (int, long, np.integer)):
                return int(job_specifier)  # is id
        else:
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
        job_id_lst = self._job_table[
            (self._job_table.project == project) & (self._job_table.job == job_specifier)].id.values
        if len(job_id_lst) == 0:
            job_id_lst = self._job_table[
                self._job_table.project.str.contains(project) & (self._job_table.job == job_specifier)].id.values
        if len(job_id_lst) == 0:
            return None
        elif len(job_id_lst) == 1:
            return job_id_lst[0]
        else:
            raise ValueError(
                "job name '{0}' in this project is not unique".format(job_specifier)
            )

    def get_child_ids(self, job_specifier, project=None, status=None):
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
        if project is None:
            project = self._project
        id_master = self.get_job_id(project=project, job_specifier=job_specifier)
        if id_master is None:
            return []
        else:
            if status is not None:
                id_lst = self._job_table[
                    (self._job_table.masterid == id_master) & (self._job_table.status == status)].id.values
            else:
                id_lst = self._job_table[(self._job_table.masterid == id_master)].id.values
            return sorted(id_lst)

    def set_job_status(self, job_specifier, status, project=None):
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
        if project is None:
            project = self._project
        job_id = self.get_job_id(project=project, job_specifier=job_specifier)
        self._job_table.loc[self._job_table.id == job_id, 'status'] = status
        db_entry = self._job_table[self._job_table.id == job_id].to_dict()
        h5io.write_hdf5(db_entry["project"][0] + db_entry["subjob"][0] + '.h5',
                        status,
                        title=db_entry["subjob"][0][1:] + '/status',
                        overwrite="update")

    def get_job_status(self, job_specifier, project=None):
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
        if project is None:
            project = self._project
        try:
            return self._job_table[
                self._job_table.id == self.get_job_id(project=project, job_specifier=job_specifier)].status.values[0]
        except KeyError:
            return None

    def get_job_working_directory(self, job_specifier, project=None):
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
        if project is None:
            project = self._project
        try:
            db_entry = self._job_table[
                self._job_table.id == self.get_job_id(project=project, job_specifier=job_specifier)].to_dict()
            if db_entry:
                job_name = db_entry["subjob"][0][1:]
                return os.path.join(
                    db_entry["project"][0],
                    job_name + "_hdf5",
                    job_name,
                )
            else:
                return None
        except KeyError:
            return None