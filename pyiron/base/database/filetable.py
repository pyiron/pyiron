import numpy as np
import os
import pandas
import datetime
import h5io
from six import with_metaclass
from pyfileindex import PyFileIndex


class Singleton(type):
    """
    Implemented with suggestions from

    http://stackoverflow.com/questions/6760685/creating-a-singleton-in-python

    """

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


def filter_function(file_name):
    return '.h5' in file_name


class FileTable(with_metaclass(Singleton)):
    def __init__(self, project):
        self._fileindex = None
        self._job_table = None
        self._project = os.path.abspath(project)
        self._columns = ['id', 'status', 'chemicalformula', 'job', 'subjob', 'projectpath', 'project', 'timestart',
                         'timestop', 'totalcputime', 'computer', 'hamilton', 'hamversion', 'parentid', 'masterid',
                         'username']
        self.force_reset()

    def force_reset(self):
        self._fileindex = PyFileIndex(path=self._project, filter_function=filter_function)
        df = pandas.DataFrame(self.init_table(fileindex=self._fileindex.dataframe))
        if len(df) != 0:
            self._job_table = df[np.array(self._columns)]
        else:
            self._job_table = pandas.DataFrame({k: [] for k in self._columns})

    def init_table(self, fileindex, working_dir_lst=None):
        if working_dir_lst is None:
            working_dir_lst = []
        fileindex = fileindex[~fileindex.is_directory]
        fileindex = fileindex.iloc[fileindex.path.values.argsort()]
        job_lst = []
        for path, mtime in zip(fileindex.path, fileindex.mtime):
            job_dict = self.get_extract(path, mtime)
            job_dict['id'] = len(working_dir_lst) + 1
            working_dir_lst.append(job_dict['project'][:-1] + job_dict['subjob'] + '_hdf5/')
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
        return {'status': get_job_status_from_file(hdf5_file=path, job_name=job),
                'chemicalformula': None,
                'job': job,
                'subjob': '/' + job,
                'projectpath': None,
                'project': os.path.dirname(path) + '/',
                'timestart': time,
                'timestop': time,
                'totalcputime': 0.0,
                'computer': None,
                'username': None,
                'parentid': None,
                'hamilton': get_hamilton_from_file(hdf5_file=path, job_name=job),
                'hamversion': get_hamilton_version_from_file(hdf5_file=path, job_name=job)}

    def add_item_dict(self, par_dict):
        par_dict = dict((key.lower(), value) for key, value in par_dict.items())
        if len(self._job_table) != 0:
            job_id = np.max(self._job_table.id.values) + 1
        else:
            job_id = 1
        default_values = {'id': job_id,
                          'status': 'initialized',
                          'chemicalformula': None,
                          'timestart': datetime.datetime.now(),
                          'computer': None,
                          'parentid': None,
                          'username': None,
                          'timestop': None,
                          'totalcputime': None,
                          'masterid': None}
        for k, v in default_values.items():
            if k not in par_dict.keys():
                par_dict[k] = v
        self._job_table = pandas.concat([self._job_table,
                                         pandas.DataFrame([par_dict])[self._columns]]).reset_index(drop=True)
        return int(par_dict['id'])

    def item_update(self, par_dict, item_id):
        if isinstance(item_id, list):
            item_id = item_id[0]
        if isinstance(item_id, str):
            item_id = float(item_id)
        for k, v in par_dict.items():
            self._job_table.loc[self._job_table.id == int(item_id), k] = v

    def delete_item(self, item_id):
        item_id = int(item_id)
        if item_id in [int(v) for v in self._job_table.id.values]:
            self._job_table = self._job_table[self._job_table.id != item_id].reset_index(drop=True)
        else:
            raise ValueError

    def get_item_by_id(self, item_id):
        item_id = int(item_id)
        return {k: list(v.values())[0] for k, v in self._job_table[self._job_table.id == item_id].to_dict().items()}

    def get_items_dict(self, item_dict, return_all_columns=True):
        df = self._job_table
        if not isinstance(item_dict, dict):
            raise TypeError
        for k, v in item_dict.items():
            if k in ['id', 'parentid', 'masterid']:
                df = df[df[k] == int(v)]
            elif "%" not in str(v):
                df = df[df[k] == v]
            else:
                df = df[df[k].str.contains(v.replace('%', ''))]
        df_dict = df.to_dict()
        if return_all_columns:
            return [{k: v[i] for k, v in df_dict.items()} for i in df_dict['id'].keys()]
        else:
            return [{'id': i} for i in df_dict['id'].values()]

    def update(self):
        self._job_table.status = [
            self._get_job_status_from_hdf5(job_id)
            for job_id in self._job_table.id.values
        ]
        self._fileindex.update()
        if len(self._job_table) != 0:
            files_lst, working_dir_lst = zip(*[[project + subjob[1:] + '.h5', project + subjob[1:] + '_hdf5']
                                               for project, subjob in zip(self._job_table.project.values,
                                                                          self._job_table.subjob.values)])
            df_new = self._fileindex.dataframe[
                ~self._fileindex.dataframe.is_directory & ~self._fileindex.dataframe.path.isin(files_lst)]
        else:
            files_lst, working_dir_lst = [], []
            df_new = self._fileindex.dataframe[~self._fileindex.dataframe.is_directory]
        if len(df_new) > 0:
            job_lst = self.init_table(fileindex=df_new, working_dir_lst=list(working_dir_lst))
            df = pandas.DataFrame(job_lst)[self._columns]
            if len(files_lst) != 0 and len(working_dir_lst) != 0:
                self._job_table = pandas.concat([self._job_table, df]).reset_index(drop=True)
            else:
                self._job_table = df

    def get_db_columns(self):
        return self.get_table_headings()

    def get_table_headings(self):
        return self._job_table.columns.values

    def job_table(
        self, 
        project=None, 
        recursive=True, 
        columns=None, 
        all_columns=False, 
        sort_by="id", 
        max_colwidth=200, 
        full_table=False, 
        job_name_contains=''
    ):
        if project is None:
            project = self._project
        if columns is None:
            columns = ["job", "project", "chemicalformula"]
        if all_columns:
            columns = self._columns
        if len(self._job_table) != 0:
            if recursive:
                df = self._job_table[self._job_table.project.str.contains(project)]
            else:
                df = self._job_table[self._job_table.project == project]
        else:
            df = self._job_table
        if full_table:
            pandas.set_option('display.max_rows', None)
            pandas.set_option('display.max_columns', None)
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
        return self.get_jobs(project=project, recursive=recursive, columns=['id'])["id"]

    def get_job_id(self, job_specifier, project=None):
        if project is None:
            project = self._project
        if isinstance(job_specifier, (int, np.integer)):
            return job_specifier  # is id

        job_specifier.replace(".", "_")
        job_id_lst = self._job_table[
            (self._job_table.project == project) & (self._job_table.job == job_specifier)].id.values
        if len(job_id_lst) == 0:
            job_id_lst = self._job_table[
                self._job_table.project.str.contains(project) & (self._job_table.job == job_specifier)].id.values
        if len(job_id_lst) == 0:
            return None
        elif len(job_id_lst) == 1:
            return int(job_id_lst[0])
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

    def get_job_working_directory(self, job_id):
        """
        Get the working directory of a particular job

        Args:
            job_id (int):

        Returns:
            str: working directory as absolute path
        """
        try:
            db_entry = self.get_item_by_id(job_id)
            if db_entry and len(db_entry) > 0:
                job_name = db_entry["subjob"][1:]
                return os.path.join(
                    db_entry["project"],
                    job_name + "_hdf5",
                    job_name,
                )
            else:
                return None
        except KeyError:
            return None

    def _get_job_status_from_hdf5(self, job_id):
        db_entry = self.get_item_by_id(job_id)
        job_name = db_entry["subjob"][1:]
        return get_job_status_from_file(
            hdf5_file=os.path.join(db_entry["project"], job_name + ".h5"),
            job_name=job_name
        )

    def get_job_status(self, job_id):
        return self._job_table[self._job_table.id == job_id].status.values[0]

    def set_job_status(self, job_id, status):
        db_entry = self.get_item_by_id(item_id=job_id)
        self._job_table.loc[self._job_table.id == job_id, 'status'] = status
        h5io.write_hdf5(db_entry["project"] + db_entry["subjob"] + '.h5',
                        status,
                        title=db_entry["subjob"][1:] + '/status',
                        overwrite="update")


def get_hamilton_from_file(hdf5_file, job_name):
    return h5io.read_hdf5(hdf5_file, job_name + '/TYPE').split(".")[-1].split("'")[0]


def get_hamilton_version_from_file(hdf5_file, job_name):
    return h5io.read_hdf5(hdf5_file, job_name + '/VERSION')


def get_job_status_from_file(hdf5_file, job_name):
    return h5io.read_hdf5(hdf5_file, job_name + '/status')
