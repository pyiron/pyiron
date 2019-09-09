# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

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
from pyiron.base.database.filetable import FileTable


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
        cls.database = FileTable(os.path.dirname(os.path.abspath(__file__)))

    def test_get_table_headings(self):
        """
        Tests get_table_headings
        Returns:
        """
        heading_list = [
            "id",
            "parentid",
            "masterid",
            "projectpath",
            "project",
            "job",
            "subjob",
            "chemicalformula",
            "status",
            "hamilton",
            "hamversion",
            "username",
            "computer",
            "timestart",
            "timestop",
            "totalcputime",
        ]
        # general headings have to be at least a part of get_table_headings
        for item in heading_list:
            self.assertTrue(item in self.database.get_table_headings())

    def test_add_item_dict(self):
        """
        Tests add_item_dict function
        Returns:
        """
        par_dict = self.add_items("BO")
        key = par_dict["id"]
        # list as parameter shall not work
        self.assertRaises(
            Exception, self.database.add_item_dict, [{"chemicalformula": "BO"}]
        )
        self.assertIsInstance(key, int)  # returned value must be int
        # added and got dict must be(almost) the same
        self.assertDictContainsSubset(par_dict, self.database.get_item_by_id(key))

    def test_item_update(self):
        """
        Tests item_update function
        Returns:
        """
        par_dict = self.add_items("BO")
        key = par_dict["id"]
        # function does only accept a dict, no list
        self.assertRaises(
            Exception, self.database.item_update, [{"job": "testing2"}], key
        )
        try:  # Function works with int, str and list, normally I would test against list, but then our project
            # would not work anymore.
            self.database.item_update({"job": "testing2"}, key)
            self.database.item_update({"job": "testing2"}, [key])
            self.database.item_update({"job": "testing2"}, str(key))
        except TypeError:
            self.fail(
                "Unexpectedly, item_update raises an Error with types of ids which should be usable"
            )

    def test_delete_item(self):
        """
        Tests delete_item function
        Returns:
        """
        par_dict = self.add_items("BO")
        key = par_dict["id"]
        self.database.delete_item(key)
        self.assertRaises(
            Exception, self.database.delete_item, [key]
        )  # use only str or int
        # self.assertRaises(Exception, self.database.get_item_by_id, key)  # ensure item does not exist anymore

    def test_get_item_by_id(self):
        """
        Tests get_item_by_id function
        Returns:
        """
        par_dict = self.add_items("BO")
        key = par_dict["id"]
        self.assertRaises(
            Exception, self.database.get_item_by_id, [key]
        )  # given key must be int or str
        # self.assertRaises(Exception, self.database.get_item_by_id,
        #                   str(key + 1))  # must give Error, if id does not exist
        self.assertIsInstance(
            self.database.get_item_by_id(key), dict
        )  # return value has to be a dict
        # added dict must (almost) be same as the got ones
        self.assertDictContainsSubset(par_dict, self.database.get_item_by_id(key))

    def test_get_items_dict_project(self):
        """
        Tests whether a query for {'project': 'Projecta%'} gives Projecta, Projecta/b/c , but not Projectas
        Returns:
        """
        par_dict = {
            "chemicalformula": "H2",
            "computer": "localhost",
            "hamilton": "VAMPE",
            "hamversion": "1.1",
            "job": "testing",
            "parentid": 0,
            "project": "Projecta/",
            "projectpath": "/TESTING",
            "status": "KAAAA",
            "timestart": datetime(2016, 5, 2, 11, 31, 4, 253377),
            "timestop": datetime(2016, 5, 2, 11, 31, 4, 371165),
            "totalcputime": 0.117788,
            "username": "User",
            "masterid": 0,
            "subjob": "testJob",
        }
        second_dict = dict(par_dict, project="Projecta/b/c/")
        third_dict = dict(par_dict, project="Projectas")
        par_dict["id"] = self.database.add_item_dict(par_dict)
        second_dict["id"] = self.database.add_item_dict(second_dict)
        third_dict["id"] = self.database.add_item_dict(third_dict)
        self.assertEqual(
            [par_dict, second_dict],
            self.database.get_items_dict({"project": "Projecta/%"}),
        )

    def test_get_items_dict_datatype(self):
        """
        Tests datatype error functionality of get_items_dict function
        Returns:
        """
        # ensures right datatype
        item_dict = ["test", "test"]
        self.assertRaises(TypeError, self.database.get_items_dict, item_dict)

    # NOT A TEST #
    def add_items(self, formula):
        """
        Simple generic helper function to add items to DB
        Args:
            formula: string for chemicalformula

        Returns:
        """
        par_dict = {
            "chemicalformula": formula,
            "computer": "localhost",
            "hamilton": "VAMPE",
            "hamversion": "1.1",
            "job": "testing",
            "parentid": 0,
            "project": "database.testing/",
            "projectpath": "/TESTING",
            "status": "KAAAA",
            "timestart": datetime(2016, 5, 2, 11, 31, 4, 253377),
            "timestop": datetime(2016, 5, 2, 11, 31, 4, 371165),
            "totalcputime": 0.117788,
            "username": "User",
            "masterid": 0,
            "subjob": "testJob",
        }
        par_dict["id"] = self.database.add_item_dict(par_dict)
        return par_dict


if __name__ == "__main__":
    unittest.main()
