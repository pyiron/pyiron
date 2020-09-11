# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import codecs
from collections import OrderedDict
from datetime import datetime
import dill as pickle
import inspect
import json
import numpy as np
import os
import pandas
from pandas.errors import EmptyDataError
from tqdm import tqdm
import types

from pyiron_base.job.generic import GenericJob
from pyiron_base.generic.hdfio import FileHDFio
from pyiron_base.master.generic import get_function_from_string
from pyiron.table.funct import (
    get_incar,
    get_sigma,
    get_total_number_of_atoms,
    get_elements,
    get_convergence_check,
    get_number_of_species,
    get_number_of_ionic_steps,
    get_ismear,
    get_encut,
    get_n_kpts,
    get_n_equ_kpts,
    get_number_of_final_electronic_steps,
    get_majority_species,
    get_job_name,
    get_job_id,
    get_energy_tot,
    get_energy_free,
    get_energy_int,
    get_energy_tot_per_atom,
    get_energy_free_per_atom,
    get_energy_int_per_atom,
    get_e_conv_level,
    get_f_states,
    get_e_band,
    get_majority_crystal_structure,
    get_equilibrium_parameters,
    get_structure,
    get_forces,
    get_magnetic_structure,
    get_average_waves,
    get_plane_waves,
    get_ekin_error,
    get_volume,
    get_volume_per_atom,
)


__author__ = "Uday Gajera, Jan Janssen, Joerg Neugebauer"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "0.0.1"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2018"


class FunctionContainer(object):
    """
    Class which is able to append, store and retreive a set of functions.

    """
    def __init__(self):
        self._user_function_dict = {}
        self._system_function_lst = [
            get_incar,
            get_sigma,
            get_total_number_of_atoms,
            get_elements,
            get_convergence_check,
            get_number_of_species,
            get_number_of_ionic_steps,
            get_ismear,
            get_encut,
            get_n_kpts,
            get_n_equ_kpts,
            get_number_of_final_electronic_steps,
            get_majority_species,
            get_job_name,
            get_job_id,
            get_energy_tot,
            get_energy_free,
            get_energy_int,
            get_energy_tot_per_atom,
            get_energy_free_per_atom,
            get_energy_int_per_atom,
            get_e_conv_level,
            get_f_states,
            get_e_band,
            get_majority_crystal_structure,
            get_equilibrium_parameters,
            get_structure,
            get_forces,
            get_magnetic_structure,
            get_average_waves,
            get_plane_waves,
            get_ekin_error,
            get_volume,
            get_volume_per_atom,
        ]
        self._system_function_dict = {
            func.__name__: False for func in self._system_function_lst
        }
        self._system_function_dict["get_job_id"] = True

    @property
    def _function_lst(self):
        return [
            funct
            for funct in self._system_function_lst
            if funct.__name__ in self._system_function_dict.keys()
            and self._system_function_dict[funct.__name__]
        ] + list(self._user_function_dict.values())

    def _to_hdf(self, hdf):
        self._to_pickle(
            hdf=hdf, key="user_function_dict", value=self._user_function_dict
        )
        self._to_pickle(
            hdf=hdf, key="system_function_dict", value=self._system_function_dict
        )

    def _from_hdf(self, hdf):
        self._user_function_dict = self._from_pickle(hdf=hdf, key="user_function_dict")
        self._system_function_dict = self._from_pickle(
            hdf=hdf, key="system_function_dict"
        )

    def __setitem__(self, key, item):
        if isinstance(item, str):
            self._user_function_dict[key] = eval(
                'lambda job: {"' + key + '":' + item + "}"
            )
        elif isinstance(item, types.FunctionType):
            self._user_function_dict[key] = lambda job: {key: item(job)}
        else:
            raise TypeError("unsupported function type!")

    def __getitem__(self, key):
        return self._user_function_dict[key]

    def __getattr__(self, name):
        if name in list(self._system_function_dict.keys()):
            self._system_function_dict[name] = True
            return self._system_function_dict[name]
        else:
            super(FunctionContainer, self).__getattr__(name)

    def __dir__(self):
        return list(self._system_function_dict.keys())

    @staticmethod
    def _to_pickle(hdf, key, value):
        hdf[key] = codecs.encode(pickle.dumps(value), "base64").decode()

    @staticmethod
    def _from_pickle(hdf, key):
        return pickle.loads(codecs.decode(hdf[key].encode(), "base64"))


class JobFilters(object):
    """
    Certain predefined job filters

    """
    @staticmethod
    def job_type(job_type):
        def filter_job_type(job):
            return job.__name__ == job_type
        return filter_job_type

    @staticmethod
    def job_name_contains(job_name_segment):
        def filter_job_name_segment(job):
            return job_name_segment in job.job_name
        return filter_job_name_segment


class PyironTable(object):
    """
    Class for easy, efficient, and pythonic analysis of data from pyiron projects

    Args:
        project (pyiron.project.Project/None): The project to analyze
        name (str): Name of the pyiron table
    """
    def __init__(self, project, name=None):
        self._project = project
        self._df = pandas.DataFrame({})
        self.convert_to_object = False
        self._name = name
        self._db_filter_function = always_true_pandas
        self._db_filter_function_str = inspect.getsource(always_true_pandas)
        self._filter_function = always_true
        self._filter_function_str = inspect.getsource(always_true)
        self._filter = JobFilters()
        self.add = FunctionContainer()
        self._csv_file = None
        if self._is_file():
            self.load()

        self.EMPTY_STR = "-"

    @property
    def filter(self):
        """
        Object containing pre-defined filter functions

        Returns:
            pyiron.table.datamining.JobFilters: The object containing the filters

        """
        return self._filter

    @property
    def _file_name_csv(self):
        if self._csv_file is None:
            return self._project.path + self.name + ".csv"
        else:
            return self._csv_file

    @property
    def _file_name_txt(self):
        return self._project.path + self.name + ".txt"

    @property
    def name(self):
        """
        Name of the table. Takes the project name if not specified

        Returns:
            str: Name of the table

        """
        if self._name is None:
            return self._project.name
        return self._name

    @property
    def db_filter_function(self):
        """
        Function to filter the a project database table before job specific functions are applied.

        The function must take a pyiron project table in the pandas.DataFrame format (project.job_table()) and return a
        boolean pandas.DataSeries with the same number of rows as the project table

        Example:

            def function(df):
                return (df["chemicalformula"=="H2"]) & (df["hamilton"=="Vasp"])

        """
        return self._db_filter_function

    @db_filter_function.setter
    def db_filter_function(self, funct):
        self._db_filter_function = funct
        try:
            self._db_filter_function_str = inspect.getsource(funct)
        except (OSError, IOError):
            pass

    @property
    def filter_function(self):
        """
        Function to filter each job before more expensive functions are applied
        """
        return self._filter_function

    @filter_function.setter
    def filter_function(self, funct):
        self._filter_function = funct
        try:
            self._filter_function_str = inspect.getsource(funct)
        except (OSError, IOError):
            pass

    def to_hdf(self):
        file = FileHDFio(file_name=self._project.path + self.name + ".h5", h5_path="/")
        self.add._to_hdf(file)

    def from_hdf(self):
        file = FileHDFio(file_name=self._project.path + self.name + ".h5", h5_path="/")
        self.add._from_hdf(file)

    def save(self, name=None):
        self._name = name
        self.to_hdf()
        self._save_csv()

    def load(self, name=None):
        self._name = name
        self.from_hdf()
        self._load_csv()

    def create_table(self, enforce_update=False, level=3, file=None, job_status_list=None):
        skip_table_update = False
        filter_funct = self.filter_function
        if job_status_list is None:
            job_status_list = ["finished"]
        if self._is_file():
            if file is None:
                file = FileHDFio(
                    file_name=self._project.path + self.name + ".h5", h5_path="/"
                )
            temp_user_function_dict, temp_system_function_dict = self._get_data_from_hdf5(
                hdf=file
            )
            job_update_lst = self._collect_job_update_lst(
                job_status_list=job_status_list,
                filter_funct=filter_funct,
                job_stored_ids=self._get_job_ids()
            )
            keys_update_user_lst = [
                key
                for key in self.add._user_function_dict.keys()
                if key not in temp_user_function_dict.keys()
            ]
            keys_update_system_lst = [
                k
                for k, v in self.add._system_function_dict.items()
                if v and not temp_system_function_dict[k]
            ]
            if (
                len(job_update_lst) == 0
                and len(keys_update_user_lst) == 0
                and keys_update_system_lst == 0
                and not enforce_update
            ):
                skip_table_update = True
        else:
            job_update_lst = self._collect_job_update_lst(
                job_status_list=job_status_list,
                filter_funct=filter_funct,
                job_stored_ids=None
            )
            keys_update_user_lst, keys_update_system_lst = [], []
        if not skip_table_update and len(job_update_lst) != 0:
            df_new_ids = self._iterate_over_job_lst(
                job_lst=job_update_lst, function_lst=self.add._function_lst, level=level
            )
        else:
            df_new_ids = pandas.DataFrame({})
        if not skip_table_update and (
            len(keys_update_user_lst) != 0 or len(keys_update_system_lst) != 0
        ):
            job_update_lst = self._collect_job_update_lst(
                job_status_list=job_status_list,
                filter_funct=filter_funct,
                job_stored_ids=None
            )
            function_lst = [
                v
                for k, v in self.add._user_function_dict.items()
                if k in keys_update_system_lst
            ] + [
                funct
                for funct in self.add._system_function_lst
                if funct.__name__ in keys_update_system_lst
            ]
            df_new_keys = self._iterate_over_job_lst(
                job_lst=job_update_lst, function_lst=function_lst, level=level
            )
        else:
            df_new_keys = pandas.DataFrame({})
        if len(self._df) > 0 and len(df_new_keys) > 0:
            self._df = pandas.concat(
                [self._df, df_new_keys], axis=1, sort=False
            ).reset_index(drop=True)
        if len(self._df) > 0 and len(df_new_ids) > 0:
            self._df = pandas.concat([self._df, df_new_ids], sort=False).reset_index(
                drop=True
            )
        elif len(df_new_ids) > 0:
            self._df = df_new_ids

    def convert_dict(self, input_dict):
        return {key: self.str_to_value(value) for key, value in input_dict.items()}

    def refill_dict(self, diff_dict_lst):
        total_key_lst = self.total_lst_of_keys(diff_dict_lst)
        for ind, sub_dict in enumerate(diff_dict_lst):
            for key in total_key_lst:
                if key not in sub_dict.keys():
                    sub_dict[key] = self.EMPTY_STR
                else:
                    sub_dict[key] = self.str_to_value(sub_dict[key])

    def col_to_value(self, col_name):
        val_lst, key_lst, ind_lst = [], [], []
        for ind, name in enumerate(self._df[col_name]):
            #             print ('name: ', ind, name)
            #             if name == self.EMPTY_STR:
            #                 continue
            name = name.split("_")
            ind_lst.append(ind)
            key_lst.append(name[0])
            val_lst.append(eval(".".join(name[1:])))
        if len(set(key_lst)) == 1:
            key = key_lst[0]
            self._df[key] = val_lst
        else:
            raise ValueError("key not unique: {}".format(set(key_lst)))

    def get_dataframe(self):
        return self._df

    def list_nodes(self):
        return list(self._df.columns)

    def list_groups(self):
        return list(set(self._df["col_0"]))

    @staticmethod
    def str_to_value(input_val):
        if not isinstance(input_val, str):
            return input_val
        else:
            try:
                return eval(input_val)
            except (TypeError, SyntaxError, NameError):
                return input_val

    @staticmethod
    def _apply_function_on_job(funct, job):
        try:
            return funct(job)
        except ValueError:
            return {}

    @staticmethod
    def total_lst_of_keys(diff_dict_lst):
        total_key_lst = []
        for sub_dict in diff_dict_lst:
            for key in sub_dict.keys():
                total_key_lst.append(key)
        return set(total_key_lst)

    def __getitem__(self, item, max_level=5):
        rename_dict = OrderedDict()
        if item in self.list_groups():
            for i in range(1, max_level):
                rename_dict["col_{}".format(i)] = "col_{}".format(i - 1)

            new_table = PyironTable(self._project[item])
            new_table._df = self._df.drop("col_0", axis=1)
            new_table._df.rename(index=str, columns=rename_dict, inplace=True)
            return new_table
        if item in self.list_nodes():
            return np.array(self._df[item])
        return None

    def __str__(self):
        return self._df.__str__()

    def __repr__(self):
        """
        Human readable string representation

        Returns:
            str: pandas Dataframe structure as string
        """
        return self._df.__repr__()

    def _is_file(self):
        return self._project is not None and os.path.isfile(self._file_name_csv)

    def _save_csv(self):
        self._df.to_csv(self._file_name_csv, index=False)

    def _load_csv(self):
        self._df = pandas.read_csv(self._file_name_csv)

    def _get_project_list(self, name, pr_len, level=3):
        lst = [self.EMPTY_STR for _ in range(level)]
        for i, p in enumerate(name.split("/")[pr_len - 1 : -1]):
            if len(lst) > i:
                lst[i] = p
        return lst

    def _get_data_from_hdf5(self, hdf):
        temp_user_function_dict = self.add._from_pickle(
            hdf=hdf, key="user_function_dict"
        )
        temp_system_function_dict = self.add._from_pickle(
            hdf=hdf, key="system_function_dict"
        )
        return temp_user_function_dict, temp_system_function_dict

    def _get_job_ids(self):
        if len(self._df) > 0:
            return self._df.job_id.values
        else:
            return np.array([])

    def _get_filtered_job_ids_from_project(self, recursive=True):
        project_table = self._project.job_table(recursive=recursive)
        filter_funct = self.db_filter_function
        return project_table[filter_funct(project_table)]["id"].tolist()

    def _apply_list_of_functions_on_job(self, job, fucntion_lst):
        diff_dict = {}
        for funct in fucntion_lst:
            funct_dict = self._apply_function_on_job(funct, job)
            for key, value in funct_dict.items():
                diff_dict[key] = value
        return diff_dict

    def _iterate_over_job_lst(self, job_lst, function_lst, level):
        pr_len = len(self._project.project_path.split("/"))
        diff_dict_lst = []
        for job_inspect in tqdm(job_lst):
            if self.convert_to_object:
                job = job_inspect.load_object()
            else:
                job = job_inspect
            diff_dict = self._apply_list_of_functions_on_job(
                job=job, fucntion_lst=function_lst
            )
            pr_lst = self._get_project_list(job.project.project_path, pr_len, level)
            for ic, col in enumerate(pr_lst):
                diff_dict["col_{}".format(ic)] = col
            diff_dict_lst.append(diff_dict)
        self.refill_dict(diff_dict_lst)
        return pandas.DataFrame(diff_dict_lst)

    def _collect_job_update_lst(self, job_status_list, filter_funct, job_stored_ids=None):
        """
        Collect jobs to update the pyiron table

        Args:
            job_status_list (list): List of job status to consider
            filter_funct (function): Filter function
            job_stored_ids (list/ None): List of already analysed job ids

        Returns:
            list: List of JobCore objects
        """
        if job_stored_ids is not None:
            job_id_lst = [
                job_id
                for job_id in self._get_filtered_job_ids_from_project()
                if job_id not in job_stored_ids
            ]
        else:
            job_id_lst = self._get_filtered_job_ids_from_project()

        job_update_lst = []
        for job_id in tqdm(job_id_lst):
            try:
                job = self._project.inspect(job_id)
            except IndexError:  # In case the job was deleted while the pyiron table is running
                job = None
            if job is not None and job.status in job_status_list and filter_funct(job):
                job_update_lst.append(job)
        return job_update_lst

    def _repr_html_(self):
        """
        Internal helper function to represent the GenericParameters object within the Jupyter Framework

        Returns:
            HTML: Jupyter HTML object
        """
        return self._df._repr_html_()


class TableJob(GenericJob):
    def __init__(self, project, job_name):
        super(TableJob, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "TableJob"
        self._analysis_project = None
        self._pyiron_table = PyironTable(project=None)
        self._enforce_update = False
        self._project_level = 0
        self.analysis_project = project.project

    @property
    def filter(self):
        return self._pyiron_table.filter

    @property
    def project_level(self):
        return self._project_level

    @project_level.setter
    def project_level(self, level):
        self._project_level = level

    @property
    def db_filter_function(self):
        return self._pyiron_table.db_filter_function

    @db_filter_function.setter
    def db_filter_function(self, funct):
        self._pyiron_table.db_filter_function = funct

    @property
    def filter_function(self):
        return self._pyiron_table.filter_function

    @filter_function.setter
    def filter_function(self, funct):
        self._pyiron_table.filter_function = funct

    @property
    def pyiron_table(self):
        return self._pyiron_table

    @property
    def ref_project(self):
        return self.analysis_project

    @ref_project.setter
    def ref_project(self, project):
        self.analysis_project = project

    @property
    def analysis_project(self):
        return self._analysis_project

    @analysis_project.setter
    def analysis_project(self, project):
        self._analysis_project = project
        self._pyiron_table = PyironTable(project=self._analysis_project)

    @property
    def add(self):
        return self._pyiron_table.add

    @property
    def convert_to_object(self):
        return self._pyiron_table.convert_to_object

    @convert_to_object.setter
    def convert_to_object(self, conv_to_obj):
        self._pyiron_table.convert_to_object = conv_to_obj

    @property
    def enforce_update(self):
        return self._enforce_update

    @enforce_update.setter
    def enforce_update(self, enforce):
        if isinstance(enforce, bool):
            if enforce:
                self._enforce_update = True
                if self.status.finished:
                    self.status.created = True
            else:
                self._enforce_update = False
        else:
            raise TypeError()

    @staticmethod
    def convert_numpy_to_list(table_dict):
        for k,v in table_dict.items():
            for k1,v1 in v.items():
                if isinstance(v1,np.ndarray):
                    v[k1] = v1.tolist()
        return table_dict


    def to_hdf(self, hdf=None, group_name=None):
        """
        Store pyiron table job in HDF5

        Args:
            hdf:
            group_name:

        """
        super(TableJob, self).to_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            hdf5_input["bool_dict"] = {
                "enforce_update": self._enforce_update,
                "convert_to_object": self._pyiron_table.convert_to_object,
            }
            self._pyiron_table.add._to_hdf(hdf5_input)
            if self._analysis_project is not None:
                hdf5_input["project"] = {
                    "path": self._analysis_project.path,
                    "user": self._analysis_project.user,
                    "sql_query": self._analysis_project.sql_query,
                    "filter": self._analysis_project._filter,
                    "inspect_mode": self._analysis_project._inspect_mode,
                }
            if self.pyiron_table._filter_function is not None:
                try:
                    hdf5_input["filter"] = inspect.getsource(
                        self.pyiron_table._filter_function
                    )
                except (OSError, IOError):
                    if self.pyiron_table._filter_function_str is not None:
                        hdf5_input["filter"] = self.pyiron_table._filter_function_str
            if self.pyiron_table._db_filter_function is not None:
                try:
                    hdf5_input["db_filter"] = inspect.getsource(
                        self.pyiron_table._db_filter_function
                    )
                except (OSError, IOError):
                    if self.pyiron_table._db_filter_function_str is not None:
                        hdf5_input["db_filter"] = self.pyiron_table._db_filter_function_str
        if len(self.pyiron_table._df) != 0:
            with self.project_hdf5.open("output") as hdf5_output:
                table_dict = self.convert_numpy_to_list(self.pyiron_table._df.to_dict())
                hdf5_output["table"] = json.dumps(table_dict)

    def from_hdf(self, hdf=None, group_name=None):
        """
        Restore pyiron table job from HDF5

        Args:
            hdf:
            group_name:
        """
        super(TableJob, self).from_hdf(hdf=hdf, group_name=group_name)
        with self.project_hdf5.open("input") as hdf5_input:
            if "project" in hdf5_input.list_nodes():
                project_dict = hdf5_input["project"]
                project = self.project.__class__(
                    path=project_dict["path"],
                    user=project_dict["user"],
                    sql_query=project_dict["sql_query"],
                )
                project._filter = project_dict["filter"]
                project._inspect_mode = project_dict["inspect_mode"]
                self.analysis_project = project
            if "filter" in hdf5_input.list_nodes():
                self.pyiron_table._filter_function_str = hdf5_input["filter"]
                self.pyiron_table.filter_function = get_function_from_string(
                    hdf5_input["filter"]
                )
            if "db_filter" in hdf5_input.list_nodes():
                self.pyiron_table._db_filter_function_str = hdf5_input["db_filter"]
                self.pyiron_table.db_filter_function = get_function_from_string(
                    hdf5_input["db_filter"]
                )
            bool_dict = hdf5_input["bool_dict"]
            self._enforce_update = bool_dict["enforce_update"]
            self._pyiron_table.convert_to_object = bool_dict["convert_to_object"]
            self._pyiron_table.add._from_hdf(hdf5_input)
        pyiron_table = os.path.join(self.working_directory, "pyirontable.csv")
        if os.path.exists(pyiron_table):
            try:
                self._pyiron_table._df = pandas.read_csv(pyiron_table)
                self._pyiron_table._csv_file = pyiron_table
            except EmptyDataError:
                pass
        else:
            with self.project_hdf5.open("output") as hdf5_output:
                if "table" in hdf5_output.list_nodes():
                    self._pyiron_table._df = pandas.DataFrame(
                        json.loads(hdf5_output["table"])
                    )

    def validate_ready_to_run(self):
        if self._analysis_project is None:
            raise ValueError("Analysis project not defined!")

    def run_static(self):
        self._create_working_directory()
        self.status.running = True
        self.update_table()
        self.status.finished = True
        self.run()

    def update_table(self, job_status_list=None):
        """
        Update the pyiron table object, add new columns if a new function was added or add new rows for new jobs

        Args:
            job_status_list (list/None): List of job status which are added to the table by default ["finished"]
        """
        if job_status_list is None:
            job_status_list = ["finished"]
        self.project.db.item_update({"timestart": datetime.now()}, self.job_id)
        with self.project_hdf5.open("input") as hdf5_input:
            self._pyiron_table.create_table(
                enforce_update=self._enforce_update,
                file=hdf5_input,
                level=self._project_level,
                job_status_list=job_status_list,
            )
        self.to_hdf()
        self._pyiron_table._df.to_csv(
            os.path.join(self.working_directory, "pyirontable.csv"), index=False
        )
        with self.project_hdf5.open("output") as hdf5_output:
            table_dict = self.convert_numpy_to_list(self.pyiron_table._df.to_dict())
            hdf5_output["table"] = json.dumps(table_dict)
        self.project.db.item_update(self._runtime(), self.job_id)

    def write_input(self):
        pass

    def get_dataframe(self):
        """

        Returns:
            pandas.Dataframe
        """
        return self.pyiron_table._df


def always_true_pandas(job_table):
    """
    A function which returns a pandas Series with all True values based on the size of the input pandas dataframe
    Args:
        job_table (pandas.DataFrame): Input dataframe

    Returns:
        pandas.Series: A series of True values

    """
    from pandas import Series
    return Series([True] * len(job_table), index=job_table.index)


def always_true(_):
    """
    A function that always returns True no matter what!

    Returns:
        bool: True

    """
    return True
