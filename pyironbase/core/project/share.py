# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import os
import subprocess

import pandas
from sqlalchemy import Column, create_engine, Integer, MetaData, String, Table
from sqlalchemy.sql import select, and_

import pyironbase.external.getent as getent
from pyironbase.core.settings.generic import Settings

try:  # If the graphical packages are not available, the GUI will not work.
    import ipywidgets as widgets
    from IPython.display import display
except ImportError:
    pass

__author__ = "Jan Janssen"
__copyright__ = "Copyright 2017, Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "development"
__date__ = "Sep 1, 2017"


s = Settings()

# Configuration:
PYIRON_USER_GROUP = 'pyiron'
PYIRON_SHARE_PATH = '/cmmc/u/'
PYIRON_SHARE_TABLE = 'share_cmmc'


def is_h5_file(name):
    return 'h5' == name.split('.')[-1]


def is_h5_dir(name):
    return 'hdf5' == name.split('_')[-1]


class MultiUserShare(object):
    def __init__(self):
        self._user_name = None
        self.user_name = s.login_user
        connection_str, user_path = self._select_environment()
        self._user_path = user_path
        engine = create_engine(connection_str)
        self._conn = engine.connect()
        self._meta_data = MetaData(bind=engine)
        self._meta_data.reflect(engine)
        self.table = Table(PYIRON_SHARE_TABLE, self._meta_data,
                           Column('id', Integer, primary_key=True, autoincrement=True),
                           Column('owner', String(20)),
                           Column('username', String(20)),
                           Column('project', String(255)),
                           Column('permission', String(1)),
                           extend_existing=True)
        self._meta_data.create_all()

    @property
    def user_name(self):
        return self._user_name

    @user_name.setter
    def user_name(self, new_user_name):
        self._validate_user_name(new_user_name)
        self._user_name = new_user_name

    @property
    def user_path(self):
        return self._user_path

    def _repr_html_(self):
        return self.list_shares_pandas(only_my_shares=False, user=None)._repr_html_()

    def _validate_user(self, user):
        if user == self.user_name:
            raise ValueError('Sharing with yourself is unnecessary!')
        if user not in self.get_users():
            raise ValueError('The user ' + user + ' does not exist within your pyiron installation.')

    def _validate_path(self, path):
        if not os.path.isabs(path):
            raise ValueError('The path should be absolute.')
        if self.user_path not in path:
            raise ValueError('You can only share files which are located in your home directory.')

    def _validate_share(self, path, user, permission):
        share_dict = self.list_shares_dict(only_my_shares=True, user=user)
        if share_dict:
            for shared_paths, shared_permission in zip(share_dict['project'], share_dict['permission']):
                if shared_paths in path:
                    if shared_permission == permission:
                        return True
                    else:
                        self.unshare_path(path, user)
        return False

    def check_share_permission(self, path, user):
        share_dict = self.list_shares_dict(only_my_shares=True, user=user)
        path_red = self._reduce_path(path)
        if share_dict:
            for shared_paths, shared_permission in zip(share_dict['project'], share_dict['permission']):
                if shared_paths in path_red and shared_permission in ['r', 'w']:
                    return shared_permission
        return False

    @staticmethod
    def _validate_permission(permission):
        if permission not in ['r']:  # ['r', 'w', 'rw']:
            raise ValueError('The only supported permissions are: [r, w]')
        if permission == 'rw':
            return 'w'
        else:
            return permission

    # def _validate_shares(self):  # Check that the shares given in the database actually exist on the filesystem

    def list_shares(self, only_my_shares=False, user=None):
        if user:
            query = select([self.table], and_(self.table.c['owner'] == self.user_name,
                           self.table.c['username'] == user))
        elif only_my_shares:
            query = select([self.table], self.table.c['owner'] == self.user_name)
        else:
            query = select([self.table])
        row = self._conn.execute(query).fetchall()
        return [dict(zip(col.keys(), col.values())) for col in row]

    def list_shares_dict(self, only_my_shares=False, user=None):
        share_lst = self.list_shares(only_my_shares=only_my_shares, user=user)
        share_dict = {}
        if share_lst:
            for dict_key in share_lst[-1].keys():
                share_dict[dict_key] = [share[dict_key] for share in share_lst]
        return share_dict

    def list_shares_pandas(self, only_my_shares=False, user=None):
        pandas.set_option('display.column_space', 50)
        return pandas.DataFrame(self.list_shares_dict(only_my_shares=only_my_shares, user=user))

    def _path_overlap(self, path, user):
        share_dict = self.list_shares_dict(only_my_shares=True, user=user)
        path_red = self._reduce_path(path)
        if share_dict:
            shared_path_lst = list(set([sub_path for shared_path in share_dict['project']
                                        for sub_path in self._get_subpaths(shared_path)]))
            shared_path_lst.sort(key=lambda shared_path: -len(shared_path))
            new_path_lst = self._get_subpaths(path_red)
            path_lst = [new_path_lst[:new_path_lst.index(path)] for path in new_path_lst if path in shared_path_lst][0]
            if path_lst:
                return path_lst
            else:
                return self._get_subpaths(path)[:-2]
        else:
            return self._get_subpaths(path)[:-2]

    def _share_subpath(self, path, user):
        path_lst = self._path_overlap(path, user)
        subprocess.call(['setfacl', '-m', 'user:' + user + ':rx'] + path_lst, stderr=subprocess.STDOUT,
                        universal_newlines=True)

    def _unshare_subpath(self, path, user):
        path_lst = self._path_overlap(path, user)
        subprocess.call(['setfacl', '-x', 'user:' + user] + path_lst[:-1], stderr=subprocess.STDOUT,
                        universal_newlines=True)

    @staticmethod
    def _get_subpaths(path):
        path_elements = path.split('/')
        return list(filter(None, ['/'.join(path_elements[:-i]) for i in range(len(path_elements))]))

    def _reduce_path(self, path):
        self._validate_path(path)
        try:
            return path.split(PYIRON_SHARE_PATH)[1]
        except ValueError:
            raise('path ', path, ' is not included in PYIRON_SHARE_PATH ', PYIRON_SHARE_PATH)

    def share_path(self, path, user, permission):
        self._validate_user(user)
        permission = self._validate_permission(permission)
        path_red = self._reduce_path(path)
        if not self._validate_share(path, user, permission):
            self._share_subpath(path, user)
            self._conn.execute(self.table.insert(
                {'owner': self.user_name, 'username': user, 'project': path_red, 'permission': permission}))
            if permission == 'w':
                permission = 'rw'
            subprocess.call(['setfacl', '-Rm', 'user:' + user + ':' + permission + 'x', path],
                            stderr=subprocess.STDOUT, universal_newlines=True)

    def unshare_path(self, path, user):
        self._validate_user(user)
        path_red = self._reduce_path(path)
        share_dict = self.list_shares_dict(only_my_shares=True, user=user)
        if not share_dict:
            return
        if path_red in share_dict['project']:
            subprocess.call(['setfacl', '-R', '-x', 'user:' + user, path], stderr=subprocess.STDOUT,
                            universal_newlines=True)
            share_id = share_dict['id'][share_dict['project'].index(path_red)]
            self._conn.execute(self.table.delete(self.table.c['id'] == int(share_id)))
            self._unshare_subpath(path, user)
        else:
            raise ValueError('No matching share found.')

    @staticmethod
    def _select_environment():
        if 'postgresql' in s.db_connection_string:
            top_lvl_lst = list(s.db_translate_dict[s.db_connection_name].keys())
            if PYIRON_SHARE_PATH in top_lvl_lst:
                return s.db_connection_string, PYIRON_SHARE_PATH
        raise TypeError('Multiuser environments require and Postgres database, but your pirmary database is: '
                        + str(s.db_connection_string))

    @staticmethod
    def _users_list():
        return sorted([dict(item)['name'] for item in getent.passwd() if dict(item)['gid'] == 12500])

    def get_users(self):
        return [user for user in self._users_list() if user != self._user_name]

    def _validate_user_name(self, new_user_name):
        if new_user_name not in self._users_list():
            raise ValueError('The user is not part of the pyiron user group.')


class ShareGUI(object):
    def __init__(self, project):
        self._ref_path = project.path
        self.project = project.copy()
        self.project._inspect_mode = True
        self.parent = None
        self.name = None
        self.fig, self.ax = None, None
        self.share = MultiUserShare()

        # Define widgets
        self.w_folder = widgets.Select(description='Folder:')
        self.w_user = widgets.Select(description='User:')
        self.w_access = widgets.Select(description='Access:')
        self.w_path = widgets.Text(name='Path: ')
        self.w_path.layout.width = '400px'

        # Arrange widgets
        self.w_tab_l = widgets.VBox([self.w_path, self.w_folder])
        self.w_tab_r = widgets.VBox([self.w_user, self.w_access])
        self.w_bottom = widgets.HTML()
        self.w_bottom.value = self.share.list_shares_pandas()._repr_html_()
        self.w_bottom.layout.width = '400px'
        # Refresh
        self.refresh_view()

        self.w_folder.observe(self.on_folder_change, names='value')
        self.w_user.observe(self.on_user_change, names='value')
        self.w_access.observe(self.on_access_change, names='value')
        menu_widgets = widgets.HBox([self.w_tab_l, self.w_tab_r])
        display(widgets.VBox([menu_widgets, self.w_bottom]))

    def on_folder_change(self, change):
        name = change['new']
        if name == '..':
            self.move_up()
            self.w_folder.value = '.'
        else:
            self.update_project(name)
        self.w_user.value = 'None'
        self.w_access.value = 'no access'

    def on_user_change(self, change):
        user = change['new']
        if user == 'None':
            self.w_access.value = 'no access'
        else:
            self.update_access_widget(user)

    def update_access_widget(self, user):
        permission = self.share.check_share_permission(self.w_path.value, user)
        if permission == 'r':
            self.w_access.value = 'read'
        elif permission == 'w':
            self.w_access.value = 'write'
        else:
            self.w_access.value = 'no access'

    def on_access_change(self, change):
        option = change['new']
        if self.w_user.value == 'None':
            return
        if option == 'read':
            self.share.share_path(self.w_path.value, self.w_user.value, 'r')
        elif option == 'write':
            self.share.share_path(self.w_path.value, self.w_user.value, 'w')
        else:
            self.share.unshare_path(self.w_path.value, self.w_user.value)
        self.w_bottom.value = self.share.list_shares_pandas()._repr_html_()

    def refresh_view(self):
        if isinstance(self.project, str):
            self._move_up()
            self.w_path.value = self.parent.path + '/' + self.name
        elif isinstance(self.project, dict):
            self._move_up()
            self.w_path.value = self.parent.path + '/' + self.name
        elif isinstance(self.project, list):
            self._move_up()
            self.w_path.value = self.parent.path + '/' + self.name
        elif self.project is None:
            raise ValueError('project is None: {}'.format(type(self.project)), self.parent)
        else:
            self.w_folder.options = ['.', '..']
            self.w_user.options = ['None'] + self.share.get_users()
            self.w_access.options = ['no access', 'read']  # ['no access', 'read', 'write']

            if hasattr(self.project, 'path'):
                self.w_path.value = self.project.path
            else:
                print('path missing: ', type(self.project))

            groups = tuple(sorted([el for el in self.project.list_groups() if not is_h5_dir(el)]))
            self.w_folder.options += groups

    def update_project(self, name):
        if name == '.':
            return
        self.name = name
        self.parent = self.project.copy()
        self.project = self.project[name]
        self.refresh_view()

    def _move_up(self):
        if hasattr(self.project, 'path'):
            self.project = self.project['..']
        else:
            self.project = self.parent
        self.w_path.value = '/'.join(self.w_path.value.split('/')[:-1])

    def move_up(self):
        self._move_up()
        self.refresh_view()
