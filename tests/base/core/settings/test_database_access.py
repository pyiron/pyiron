"""
Normally, test-based development is done in reversed order unlike this DatabaseClass.
First the Unittest should be created and then the wanted class which should pass all the tests.
Because of that, this Unittest is based on the DatabaseClass making the executed test cases a bit strange
and specific extra for this class.

Murat Han Celik
"""
import unittest
import os
from datetime import datetime
from random import choice
from string import ascii_uppercase
from pyiron.base.core.settings.database import DatabaseAccess


class TestDatabaseAccess(unittest.TestCase):
    """
    Standard Unittest of the DatabaseAccess class
    """
    @classmethod
    def setUpClass(cls):
        """
        Set up whole class for testing
        Returns:
        """
        # we assume everything working on sqlite, should also work on postgres in sqlalchemy
        cls.database = DatabaseAccess('sqlite:///test_database.db', 'simulation')

    @classmethod
    def tearDownClass(cls):
        """
        Tear down whole class after testing
        Returns:
        """
        cls.database.conn.close()
        os.remove('test_database.db')

    def tearDown(self):
        """
        Deletes all entries after every tested function
        Returns:
        """
        self.database.conn.execute('delete from simulation')

    def test_get_table_headings(self):
        """
        Tests get_table_headings
        Returns:
        """
        heading_list = ['id',
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
        # general headings have to be at least a part of get_table_headings
        for item in heading_list:
            self.assertTrue(item in self.database.get_table_headings())

    def test_get_items_sql(self):
        """
        Tests get_items_sql function
        Returns:
        """
        self.add_items('Blub')
        self.add_items('Bluk')
        # has to return a list
        self.assertIsInstance(self.database.get_items_sql("chemicalformula LIKE 'Blu%'"), list)
        self.assertRaises(Exception, self.database.get_items_sql, "A Wrong where Clause")  # Where clause must be right
        # A valid sqlstatement should also return a valid list
        self.assertIsInstance(self.database.get_items_sql(where_condition="",
                                                          sql_statement="select * from simulation"), list)
        par_dict = self.add_items('BO')
        key = par_dict['id']
        # be sure that get_items_sql returns right result with right statement
        self.assertDictContainsSubset(par_dict, self.database.get_items_sql(where_condition="",
                                                                            sql_statement='select * from simulation '
                                                                                          'where id=%s' % key)[-1])

    def test_get_items_sql_like_regex(self):
        """
        Tests the regex functionality of 'like'
        Returns:
        """
        elem1 = self.add_items('H4Ni2')
        elem2 = self.add_items('H2')
        elem3 = self.add_items('H6')
        elem4 = self.add_items('HeNi2')
        elem5 = self.add_items('H12Ni5')
        elem6 = self.add_items('H12')
        elem7 = self.add_items('He2')

        # H([0-9]*) matches H2, H6 and H12
        self.assertEqual([elem2, elem3, elem6], self.database.get_items_sql(r"chemicalformula like 'H([0-9]*)'"))
        # He(\d)*(Ni)?\d* matches HeNi2, He2
        self.assertEqual([elem4, elem7], self.database.get_items_sql(r"chemicalformula like 'He(\d)*(Ni)?\d*'"))
        # H\d*Ni\d* matches H4Ni2, H12Ni5
        self.assertEqual([elem1, elem5], self.database.get_items_sql(r"chemicalformula like 'H\d*Ni\d*'"))
        # assert that not something random really is in the Database, recommended by Samuel Hartke
        # Murat: 'Just ignore the line!'
        self.assertEqual([], self.database.get_items_sql(r"chemicalformula like 'B\d[a-z]'"))

    def test_add_item_dict(self):
        """
        Tests add_item_dict function
        Returns:
        """
        par_dict = self.add_items('BO')
        key = par_dict['id']
        # list as parameter shall not work
        self.assertRaises(Exception, self.database.add_item_dict, [{'chemicalformula': 'BO'}])
        self.assertIsInstance(key, int)  # returned value must be int
        # added and got dict must be(almost) the same
        self.assertDictContainsSubset(par_dict, self.database.get_item_by_id(key))

    def test_item_update(self):
        """
        Tests item_update function
        Returns:
        """
        par_dict = self.add_items('BO')
        key = par_dict['id']
        # function does only accept a dict, no list
        self.assertRaises(Exception, self.database.item_update, [{'job': 'testing2'}], key)
        try:  # Function works with int, str and list, normally I would test against list, but then our project
            # would not work anymore.
            self.database.item_update({'job': 'testing2'}, key)
            self.database.item_update({'job': 'testing2'}, [key])
            self.database.item_update({'job': 'testing2'}, str(key))
        except TypeError:
            self.fail('Unexpectedly, item_update raises an Error with types of ids which should be usable')

    def test_delete_item(self):
        """
        Tests delete_item function
        Returns:
        """
        par_dict = self.add_items('BO')
        key = par_dict['id']
        self.database.delete_item(key)
        self.assertRaises(Exception, self.database.delete_item, [key])  # use only str or int
        # self.assertRaises(Exception, self.database.get_item_by_id, key)  # ensure item does not exist anymore

    def test_get_item_by_id(self):
        """
        Tests get_item_by_id function
        Returns:
        """
        par_dict = self.add_items('BO')
        key = par_dict['id']
        self.assertRaises(Exception, self.database.get_item_by_id, [key])  # given key must be int or str
        # self.assertRaises(Exception, self.database.get_item_by_id,
        #                   str(key + 1))  # must give Error, if id does not exist
        self.assertIsInstance(self.database.get_item_by_id(key), dict)  # return value has to be a dict
        # added dict must (almost) be same as the got ones
        self.assertDictContainsSubset(par_dict, self.database.get_item_by_id(key))

    def test_get_items_dict_and(self):
        """
        Tests the 'and' functionality of get_items_dict function
        Returns:
        """
        self.add_items('Blub')
        # tests general and statements
        item_dict = {'hamilton': 'VAMPE',
                     'hamversion': '1.1'}
        self.assertEqual(self.database.get_items_dict(item_dict),
                         self.database.get_items_sql("hamilton='VAMPE' and hamversion='1.1'"))

    def test_get_items_dict_project(self):
        """
        Tests whether a query for {'project': 'Projecta%'} gives Projecta, Projecta/b/c , but not Projectas
        Returns:
        """
        par_dict = {'chemicalformula': 'H2',
                    'computer': 'localhost',
                    'hamilton': 'VAMPE',
                    'hamversion': '1.1',
                    'job': 'testing',
                    'parentid': 0,
                    'project': 'Projecta/',
                    'projectpath': '/TESTING',
                    'status': 'KAAAA',
                    'timestart': datetime(2016, 5, 2, 11, 31, 4, 253377),
                    'timestop': datetime(2016, 5, 2, 11, 31, 4, 371165),
                    'totalcputime': 0.117788,
                    'username': 'User',
                    'masterid': 0,
                    'subjob': 'testJob'}
        second_dict = dict(par_dict, project='Projecta/b/c/')
        third_dict = dict(par_dict, project='Projectas')
        par_dict['id'] = self.database.add_item_dict(par_dict)
        second_dict['id'] = self.database.add_item_dict(second_dict)
        third_dict['id'] = self.database.add_item_dict(third_dict)
        self.assertEqual([par_dict, second_dict], self.database.get_items_dict({'project': 'Projecta/%'}))

    def test_get_items_dict_or(self):
        """
        Tests 'or' functionality of get_items_dict function
        Returns:
        """
        self.add_items('Blub')
        self.add_items('Blab')
        # tests an example or statement
        item_dict = {'chemicalformula': ['Blub', 'Blab']}
        # assert that both the sql and non-sql methods give the same result
        sql_db = self.database.get_items_sql("chemicalformula='Blub' or chemicalformula='Blab'")
        dict_db = self.database.get_items_dict(item_dict)
        for item in sql_db:
            self.assertTrue(item in dict_db)

    def test_get_items_dict_like(self):
        """
        Tests 'like' functionality of get_items_dict function
        Returns:
        """
        self.add_items('Blub')
        # tests an example like statement
        item_dict = {'status': '%AA%'}
        # assert that both the sql and non-sql methods give the same result
        sql_db = self.database.get_items_sql("status like '%AA%'")
        dict_db = self.database.get_items_dict(item_dict)
        for item in sql_db:
            self.assertTrue(item in dict_db)

    def test_get_items_dict_datatype(self):
        """
        Tests datatype error functionality of get_items_dict function
        Returns:
        """
        # ensures right datatype
        item_dict = ['test', 'test']
        self.assertRaises(TypeError, self.database.get_items_dict, item_dict)

    def test_z_add_column(self):
        """
        Tests add_column function
        Name includes a z so that it is run last. Altering a table can lead
        the other tests to fail.
        Returns:
        """
        self.add_items('blub')
        column = 'myColumn5'

        if column not in self.database.get_table_headings():
            self.database.add_column(column, 'varchar(50)')
        self.assertRaises(Exception, self.database.add_column, column)  # cannot add column with same name
        self.assertTrue(column in self.database.get_table_headings())  # see whether myColumn has been included
        try:
            # list should be usable, but function will just take last element of lists
            second_column = ''.join(choice(ascii_uppercase) + str(i) for i in range(12))  # random generator for columns
            self.database.add_column([second_column], ['varchar(50)'])
            self.database.add_column([second_column + "2"], 'varchar(50)')
        except TypeError:
            self.fail("Unexpectedly add_column cannot take lists as parameter.")
        self.assertRaises(Exception, self.database.add_column, ['mycolumn'], 10)  # cannot use int as params

    # NOT A TEST #
    def add_items(self, formula):
        """
        Simple generic helper function to add items to DB
        Args:
            formula: string for chemicalformula

        Returns:
        """
        par_dict = {'chemicalformula': formula,
                    'computer': 'localhost',
                    'hamilton': 'VAMPE',
                    'hamversion': '1.1',
                    'job': 'testing',
                    'parentid': 0,
                    'project': 'database.testing/',
                    'projectpath': '/TESTING',
                    'status': 'KAAAA',
                    'timestart': datetime(2016, 5, 2, 11, 31, 4, 253377),
                    'timestop': datetime(2016, 5, 2, 11, 31, 4, 371165),
                    'totalcputime': 0.117788,
                    'username': 'User',
                    'masterid': 0,
                    'subjob': 'testJob'}
        par_dict['id'] = self.database.add_item_dict(par_dict)
        return par_dict


if __name__ == '__main__':
    unittest.main()
