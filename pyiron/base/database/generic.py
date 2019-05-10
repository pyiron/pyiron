# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import numpy as np
import re
import time
from datetime import datetime
from sqlalchemy import Column, create_engine, DateTime, Float, Integer, MetaData, String, Table, text, and_, or_
from sqlalchemy.pool import NullPool
from sqlalchemy.sql import select
from sqlalchemy.exc import OperationalError, DatabaseError

"""
DatabaseAccess class deals with accessing the database
"""

__author__ = "Murat Han Celik"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH" \
                " - Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


class AutorestoredConnection:
    def __init__(self, engine):
        self.engine = engine
        self._conn = None

    def execute(self, *args, **kwargs):
        try:
            if not self._conn or self._conn.closed:
                self._conn = self.engine.connect()
            result = self._conn.execute(*args, **kwargs)
        except OperationalError:
            time.sleep(5)
            result = self.execute(*args, **kwargs)
        return result

    def close(self):
        if self._conn is not None:
            self._conn.close()


class DatabaseAccess(object):
    """
    A core element of PyIron, which generally deals with accessing the database: getting, sending, changing some data
    to the db.

    Args:
        connection_string (str): SQLalchemy connection string which specifies the database to connect to
                                 typical form: dialect+driver://username:password@host:port/database
                                 example: 'postgresql://scott:tiger@cmcent56.mpie.de/mdb'
        table_name (str): database table name, a simple string like: 'simulation'

    Murat Han Celik
    """
    def __init__(self, connection_string, table_name):
        """
        Initialize the Database connection

        Args:
            connection_string (str): SQLalchemy connection string which specifies the database to connect to
                                     typical form: dialect+driver://username:password@host:port/database
                                     example: 'postgresql://scott:tiger@cmcent56.mpie.de/mdb'
            table_name (str): database table name, a simple string like: 'simulation'
        """
        self.table_name = table_name
        self._keep_connection = False
        self._sql_lite = 'sqlite' in connection_string
        try:
            if not self._sql_lite:
                self._engine = create_engine(connection_string,
                                             connect_args={'connect_timeout': 15},
                                             poolclass=NullPool)
                self.conn = AutorestoredConnection(self._engine)
            else:
                self._engine = create_engine(connection_string)
                self.conn = self._engine.connect()
                self.conn.connection.create_function("like", 2, self.regexp)
                self._keep_connection = True
        except Exception as except_msg:
            raise ValueError("Connection to database failed: " + str(except_msg))

        self.__reload_db()
        self.simulation_table = Table(str(table_name), self.metadata,
                                      Column('id', Integer, primary_key=True, autoincrement=True),
                                      Column('parentid', Integer),
                                      Column('masterid', Integer),
                                      Column('projectpath', String(50)),
                                      Column('project', String(255)),
                                      Column('job', String(50)),
                                      Column('subjob', String(255)),
                                      Column('chemicalformula', String(30)),
                                      Column('status', String(20)),
                                      Column('hamilton', String(20)),
                                      Column('hamversion', String(50)),
                                      Column('username', String(20)),
                                      Column('computer', String(100)),
                                      Column('timestart', DateTime),
                                      Column('timestop', DateTime),
                                      Column('totalcputime', Float),
                                      extend_existing=True)
        self.metadata.create_all()
        self._viewer_mode = False

    @property
    def viewer_mode(self):
        """
        Get viewer_mode - if viewer_mode is enable pyiron has read only access to the database.

        Returns:
            bool: returns TRUE when viewer_mode is enabled
        """
        return self._viewer_mode

    @viewer_mode.setter
    def viewer_mode(self, value):
        """
        Set viewer_mode - if viewer_mode is enable pyiron has read only access to the database.

        Args:
            value (bool): TRUE or FALSE
        """
        if isinstance(value, bool):
            self._viewer_mode = value
        else:
            raise TypeError('Viewmode can only be TRUE or FALSE.')

    # Internal functions
    def __del__(self):
        """
        Close database connection

        Returns:

        """
        self.conn.close()

    def __reload_db(self):
        """
        Reload database

        Returns:

        """
        self.metadata = MetaData(bind=self._engine)
        self.metadata.reflect(self._engine)

    @staticmethod
    def regexp(expr, item):
        """
        Regex function for SQLite
        Args:
            expr: str, regex expression
            item: str, item which needs to be checked

        Returns:

        """
        expr = expr.replace('%', '(.)*')
        expr = expr.replace('_', '.')
        expr = '^' + expr
        if expr[-1] is not '%':
            expr += '$'
        reg = re.compile(expr)
        if item is not None:
            return reg.search(item) is not None

    # Table functions
    def get_table_headings(self, table_name=None):
        """
        Get column names

        Args:
            table_name (str): simple string of a table_name like: 'jobs_username'

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
        if table_name is None:
            table_name = self.table_name
        self.__reload_db()
        try:
            simulation_list = Table(str(table_name), self.metadata, autoload=True, autoload_with=self._engine)
        except Exception:
            raise ValueError(str(table_name) + " does not exist")
        return [column.name for column in iter(simulation_list.columns)]

    def add_column(self, col_name, col_type):
        """
        Add an additional column - required for modification on the database

        Args:
            col_name (str, list): name of the new column, normal string like: 'myColumn'
            col_type (str, list: SQL type of the new column, SQL type like: 'varchar(50)'

        Returns:

        """
        if not self._viewer_mode:
            if isinstance(col_name, list):
                col_name = col_name[-1]
            if isinstance(col_type, list):
                col_type = col_type[-1]
            self._engine.execute('ALTER TABLE %s ADD COLUMN %s %s' % (self.simulation_table.name, col_name, col_type))
        else:
            raise PermissionError('Not avilable in viewer mode.')

    def change_column_type(self, col_name, col_type):
        """
        Modify data type of an existing column - required for modification on the database

        Args:
            col_name (str, list): name of the new column, normal string like: 'myColumn'
            col_type (str, list: SQL type of the new column, SQL type like: 'varchar(50)'

        Returns:

        """
        if not self._viewer_mode:
            if isinstance(col_name, list):
                col_name = col_name[-1]
            if isinstance(col_type, list):
                col_type = col_type[-1]
            self._engine.execute('ALTER TABLE %s ALTER COLUMN %s TYPE %s' % (self.simulation_table.name,
                                                                             col_name,
                                                                             col_type))
        else:
            raise PermissionError('Not avilable in viewer mode.')

    def get_items_sql(self, where_condition=None, sql_statement=None):
        """
        Submit an SQL query to the database

        Args:
            where_condition (str): SQL where query, query like: "project LIKE 'lammps.phonons.Ni_fcc%'"
            sql_statement (str): general SQL query, normal SQL statement

        Returns:
            list: get a list of dictionaries, where each dictionary represents one item of the table like:
                 [{u'chemicalformula': u'BO',
                  u'computer': u'localhost',
                  u'hamilton': u'VAMPS',
                  u'hamversion': u'1.1',
                  u'id': 1,
                  u'job': u'testing',
                  u'masterid': None,
                  u'parentid': 0,
                  u'project': u'database.testing',
                  u'projectpath': u'/TESTING',
                  u'status': u'KAAAA',
                  u'subjob': u'testJob',
                  u'timestart': u'2016-05-02 11:31:04.253377',
                  u'timestop': u'2016-05-02 11:31:04.371165',
                  u'totalcputime': 0.117788,
                  u'username': u'User'},
                 {u'chemicalformula': u'BO',
                  u'computer': u'localhost',
                  u'hamilton': u'VAMPS',
                  u'hamversion': u'1.1',
                  u'id': 2,
                  u'job': u'testing',
                  u'masterid': 0,
                  u'parentid': 0,
                  u'project': u'database.testing',
                  u'projectpath': u'/TESTING',
                  u'status': u'KAAAA',
                  u'subjob': u'testJob',
                  u'timestart': u'2016-05-02 11:31:04.253377',
                  u'timestop': u'2016-05-02 11:31:04.371165',
                  u'totalcputime': 0.117788,
                  u'username': u'User'}.....]
        """

        if where_condition:
            where_condition = where_condition.replace(
                'like', 'similar to') if self._engine.dialect.name is "postgresql" else where_condition
            try:
                query = "select * from " + self.table_name + " where " + where_condition
                query.replace('%', '%%')
                result = self.conn.execute(text(query))
            except Exception as except_msg:
                print("EXCEPTION in get_items_sql: ", except_msg)
                raise ValueError("EXCEPTION in get_items_sql: ", except_msg)
        elif sql_statement:
            sql_statement = sql_statement.replace(
                'like', 'similar to') if self._engine.dialect.name is "postgresql" else sql_statement
            # TODO: make it save against SQL injection
            result = self.conn.execute(text(sql_statement))
        else:
            result = self.conn.execute(text("select * from " + self.table_name))
        row = result.fetchall()
        if not self._keep_connection:
            self.conn.close()

        # change the date of str datatype back into datetime object
        output_list = []
        for col in row:
            # ensures working with db entries, which are camel case
            timestop_index = [item.lower() for item in col.keys()].index('timestop')
            timestart_index = [item.lower() for item in col.keys()].index('timestart')
            tmp_values = col.values()
            if (col.values()[timestop_index] and col.values()[timestart_index]) is not None:
                # changes values
                try:
                    tmp_values[timestop_index] = datetime.strptime(
                        str(tmp_values[timestop_index]), '%Y-%m-%d %H:%M:%S.%f')
                    tmp_values[timestart_index] = datetime.strptime(
                        str(tmp_values[timestart_index]), '%Y-%m-%d %H:%M:%S.%f')
                except ValueError:
                    print("error in: ", str(col))
            output_list += [dict(zip(col.keys(), tmp_values))]
        return output_list

    # Item functions
    def add_item_dict(self, par_dict):
        """
        Create a new database item

        Args:
            par_dict (dict): Dictionary with the item values and column names as keys, like:
                              {'chemicalformula': 'BO',
                             'computer': 'localhost',
                             'hamilton': 'VAMPS',
                             'hamversion': '1.1',
                             'job': 'testing',
                             'subjob' : 'SubJob',
                             'parentid': 0L,
                             'myCol': 'Blubbablub',
                             'project': 'database.testing',
                             'projectpath': '/root/directory/tmp',
                             'status': 'KAAAA',
                             'timestart': datetime(2016, 5, 2, 11, 31, 4, 253377),
                             'timestop': datetime(2016, 5, 2, 11, 31, 4, 371165),
                             'totalcputime': 0.117788,
                             'username': 'Test'}

        Returns:
            int: Database ID of the item created as an int, like: 3
        """
        if not self._viewer_mode:
            try:
                par_dict = dict((key.lower(), value) for key, value in par_dict.items())           # make keys lowercase
                result = self.conn.execute(self.simulation_table.insert(par_dict))
                if not self._keep_connection:
                    self.conn.close()
                return result.inserted_primary_key[-1]
            except Exception as except_msg:
                raise ValueError("Error occurred: " + str(except_msg))
        else:
            raise PermissionError('Not avilable in viewer mode.')

    def __get_items(self, col_name, var):
        """
        Get multiple items from the database

        Args:
            col_name (str): column to query for, like :  'id'
            var (str, int): value of the specific column, like: '2'

        ----> __get_items('id', '2')

        Returns:
            dict: Dictionary where the key is the column name, like:
                    [{'chemicalformula': u'BO',
                      'computer': u'computer',
                      'hamilton': u'VAMPS',
                      'hamversion': u'1.1',
               ------>'id': 2,
                      'job': u'testing',
                      'parentid': 0,
                      'project': u'database.testing',
                      'projectpath': u'/root/directory/tmp',
                      'samucol': None,
                      'status': u'Testing',
                      'timestart': datetime.datetime(2016, 5, 2, 11, 31, 4, 253377),
                      'timestop': datetime.datetime(2016, 5, 2, 11, 31, 4, 371165),
                      'totalcputime': 0.117788,
                      'username': u'Test'}]
        """
        try:
            if type(var) is list:
                var = var[-1]
            query = select([self.simulation_table], self.simulation_table.c[str(col_name)] == var)
        except Exception:
            raise ValueError("There is no Column named: " + col_name)
        try:
            result = self.conn.execute(query)
        except (OperationalError, DatabaseError):
            if not self._sql_lite:
                self.conn = AutorestoredConnection(self._engine)
            else:
                self.conn = self._engine.connect()
                self.conn.connection.create_function("like", 2, self.regexp)
            result = self.conn.execute(query)
        row = result.fetchall()
        if not self._keep_connection:
            self.conn.close()
        return [dict(zip(col.keys(), col.values())) for col in row]

    def item_update(self, par_dict, item_id):
        """
        Modify Item in database

        Args:
            par_dict (dict): Dictionary of the parameters to be modified,, where the key is the column name.
                            {'job' : 'maximize',
                             'subjob' : 'testing',
                             ........}
            item_id (int, list): Database Item ID (Integer) - '38'  can also be [38]

        Returns:

        """
        if not self._viewer_mode:
            if type(item_id) is list:
                item_id = item_id[-1]                                           # sometimes a list is given, make it int
            if np.issubdtype(type(item_id), np.integer):
                item_id = int(item_id)
            # all items must be lower case, ensured here
            par_dict = dict((key.lower(), value) for key, value in par_dict.items())
            query = self.simulation_table.update(self.simulation_table.c['id'] == item_id).values()
            try:
                self.conn.execute(query, par_dict)
            except (OperationalError, DatabaseError):
                if not self._sql_lite:
                    self.conn = AutorestoredConnection(self._engine)
                else:
                    self.conn = self._engine.connect()
                    self.conn.connection.create_function("like", 2, self.regexp)

                self.conn.execute(query, par_dict)
            if not self._keep_connection:
                self.conn.close()
        else:
            raise PermissionError('Not avilable in viewer mode.')

    def delete_item(self, item_id):
        """
        Delete Item from database

        Args:
            item_id (int): Databse Item ID (Integer), like: 38

        Returns:

        """
        if not self._viewer_mode:
            self.conn.execute(self.simulation_table.delete(self.simulation_table.c['id'] == int(item_id)))
            if not self._keep_connection:
                self.conn.close()
        else:
            raise PermissionError('Not avilable in viewer mode.')

    # Shortcut
    def get_item_by_id(self, item_id):
        """
        Get item from database by searching for a specific item Id.

        Args:
            item_id (int): Databse Item ID (Integer), like: 38

        Returns:
            dict: Dictionary where the key is the column name, like:
                    {'chemicalformula': u'BO',
                     'computer': u'localhost',
                     'hamilton': u'VAMPS',
                     'hamversion': u'1.1',
                     'id': 1,
                     'job': u'testing',
                     'masterid': None,
                     'parentid': 0,
                     'project': u'database.testing',
                     'projectpath': u'/root/directory/tmp',
                     'status': u'KAAAA',
                     'subjob': u'SubJob',
                     'timestart': datetime.datetime(2016, 5, 2, 11, 31, 4, 253377),
                     'timestop': datetime.datetime(2016, 5, 2, 11, 31, 4, 371165),
                     'totalcputime': 0.117788,
                     'username': u'Test'}
        """
        # convert item_id to int type
        # needed since psycopg2 gives otherwise an error for np.int64 type (bigint in database)
        if item_id is None:
            return None
        if isinstance(item_id, (str, float)):
            item_id = int(item_id)
        if np.issubdtype(type(item_id), np.integer):
            try:
                return self.__get_items('id', int(item_id))[-1]
            except TypeError as except_msg:
                raise TypeError("Wrong data type given as parameter. item_id has to be Integer or String: ", except_msg)
            except IndexError as except_msg:
                raise IndexError("Error when trying to find elements by given Job ID: ", except_msg)
        else:
            raise TypeError('THE SQL database ID has to be an integer.')

    def query_for_element(self, element):
        return or_(
            *[self.simulation_table.c['chemicalformula'].like('%' + element +
                                                              '[ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]%'),
              self.simulation_table.c['chemicalformula'].like('%' + element)])

    def get_items_dict(self, item_dict, return_all_columns=True):
        """

        Args:
            item_dict (dict): a dict type, which has a certain syntax for this function:
                              a normal dict like {'hamilton': 'VAMPE', 'hamversion': '1.1'} has similarities with a
                              simple query like
                                  select * from table_name where hamilton = 'VAMPE AND hamversion = '1.1'
                              as seen it puts an AND for every key, value combination in the dict and searches for it.

                              another syntax is for an OR statement, simply: {'hamilton': ['VAMPE', 'LAMMPS']}, the
                              query would be:
                                  select * from table_name where hamilton = 'VAMPE' OR hamilton = 'LAMMPS'

                              and lastly for a LIKE statement, simply: {'project': 'database.%'}, the query would be
                                  select * from table_name where project LIKE 'database.%'
                              that means you can simply add the syntax for a like statement like '%' and it will
                              automatically operate a like-search

                              of course you can also use a more complex select method, with everything in use:
                                  {'hamilton': ['VAMPE', 'LAMMPS'],
                                   'project': 'databse%',
                                   'hamversion': '1.1'}
                                  select * from table_name where (hamilton = 'VAMPE' Or hamilton = 'LAMMPS') AND
                                      (project LIKE 'database%') AND hamversion = '1.1'
            return_all_columns (bool): return all columns or only the 'id' - still the format stays the same.

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
        if not isinstance(item_dict, dict):
            raise TypeError("Wrong DataType! Only Dicts are usable!")
        and_statement = []      # list for the whole sqlalchemy statement
        # here we go through all keys and values of item_dict
        for key, value in item_dict.items():
            # if a value of item_dict is a list, we have to make an or statement of it
            if key == 'element_lst':
                part_of_statement = [self.query_for_element(element=element) for element in value]
            elif isinstance(value, list):
                or_statement = [self.simulation_table.c[str(key)] == element
                                if '%' not in element
                                else self.simulation_table.c[str(key)].like(element)
                                for element in value]
                # here we wrap the given values in an sqlalchemy-type or_statement
                part_of_statement = [or_(*or_statement)]
            else:
                if '%' not in str(value):
                    part_of_statement = [self.simulation_table.c[str(key)] == value]
                else:
                    part_of_statement = [self.simulation_table.c[str(key)].like(value)]
            # here all statements are wrapped together for the and statement
            and_statement += part_of_statement
        if return_all_columns:
            query = select([self.simulation_table], and_(*and_statement))
        else:
            query = select([self.simulation_table.columns['id']], and_(*and_statement))
        try:
            result = self.conn.execute(query)
        except (OperationalError, DatabaseError):
            if not self._sql_lite:
                self.conn = AutorestoredConnection(self._engine)
            else:
                self.conn = self._engine.connect()
                self.conn.connection.create_function("like", 2, self.regexp)

            result = self.conn.execute(query)
        row = result.fetchall()
        if not self._keep_connection:
            self.conn.close()
        return [dict(zip(col.keys(), col.values())) for col in row]

