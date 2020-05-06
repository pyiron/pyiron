# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from __future__ import print_function
from pyiron.base.settings.generic import Settings

"""
Functions to update existing pyiron installations - mainly modify the database columns. 
"""

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


def database():
    """
    Convenience function to update an existing (older) version of the database to the latest version, by modifying the
    database columns. This is only possible if no other pyiron session is accessing the database. Therefore the script
    might take some time to be executed successfully.
    """
    s = Settings()
    s.open_connection()
    db = s.database
    try:
        if "projectPath".lower() not in db.get_table_headings(db.table_name):
            print("add missing column: " + "projectPath")
            db.add_column(col_name="projectPath", col_type="varchar(255)")
        if "subJob".lower() not in db.get_table_headings(db.table_name):
            print("add missing column: " + "subJob")
            db.add_column(col_name="subJob", col_type="varchar(255)")
        else:
            print("change data type of subJob")
            db.change_column_type(col_name="subJob", col_type="varchar(255)")
        if "masterID".lower() not in db.get_table_headings(db.table_name):
            print("add missing column: " + "masterid")
            db.add_column(col_name="masterid", col_type="bigint")

        if "hamversion" in db.get_table_headings(db.table_name):
            print("change data type hamversion")
            db.change_column_type(col_name="hamversion", col_type="varchar(50)")

        if "job" in db.get_table_headings(db.table_name):
            print("change data type job")
            db.change_column_type(col_name="job", col_type="varchar(50)")
        print(db.table_name, " - database successful updated")
    except ValueError:
        print(db.table_name, " - database failed")

    print("database update done")


if __name__ == "__main__":
    database()
