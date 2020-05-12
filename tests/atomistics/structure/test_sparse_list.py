# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import unittest
from pyiron.atomistics.structure.sparse_list import SparseList, SparseArray


class TestSparseList(unittest.TestCase):
    def setUp(self):
        self.aList = SparseList({1: True, 4: True}, length=7)
        self.cList = SparseList({1: 0.2, 3: 0.5}, default=0, length=5)

    def test__len__(self):
        self.assertEqual(len(self.aList), 7)
        self.assertEqual(len(self.cList), 5)

    def test_list(self):
        self.assertEqual(self.aList.list(), [None, True, None, None, True, None, None])
        self.assertEqual(self.cList.list(), [0, 0.2, 0, 0.5, 0])

    def test_keys(self):
        self.assertEqual(list(self.aList.keys()), [1, 4])
        self.assertEqual(list(self.cList.keys()), [1, 3])

    def test_items(self):
        self.assertEqual(list(self.aList.items()), [(1, True), (4, True)])
        self.assertEqual(list(self.cList.items()), [(1, 0.2), (3, 0.5)])

    def test__iter__(self):
        self.assertEqual([a.index for a in self.aList], [1, 4])
        self.assertEqual([c for c in self.cList], [0, 0.2, 0, 0.5, 0])

    def test__getitem__(self):
        self.assertEqual(self.aList[1], True)
        self.assertEqual(list(self.aList[2:4].keys()), [])
        self.assertEqual(list(self.aList[2:5].items()), [(2, True)])
        self.assertEqual(self.aList[-1], None)
        self.assertTrue(isinstance(self.aList[2:4], SparseList))

        self.assertEqual(self.cList[1], 0.2)
        self.assertEqual(list(self.cList[2:4].keys()), [1])
        self.assertEqual(list(self.cList[2:5].items()), [(1, 0.5)])
        self.assertEqual(self.cList[-1], 0)
        self.assertTrue(isinstance(self.cList[2:4], SparseList))

    def test__setitem__(self):
        self.aList[::2] = True
        self.assertEqual(self.aList.list(), [True, True, True, None, True, None, True])
        self.aList[4] = False
        self.assertEqual(self.aList.list(), [True, True, True, None, False, None, True])

    def test__delitem__(self):
        del self.aList[1:3]
        self.assertEqual(self.aList.list(), [None, None, True, None, None])

        del self.aList[:]
        self.assertEqual(self.aList.list(), [])

        del self.cList[2]
        self.assertEqual(self.cList.list(), [0, 0.2, 0.5, 0])

    def test__add__(self):
        bList = self.aList + self.aList
        self.assertEqual(len(bList), 14)
        bList += self.aList
        self.assertEqual(len(bList), 21)
        self.assertEqual(self.aList.list(), [None, True, None, None, True, None, None])

    def test__mul__(self):
        self.assertEqual(len(200 * self.aList), 1400)
        self.assertEqual(len(self.aList * 2), 14)

    def test__str__(self):
        self.assertEqual(self.aList.__str__(), "[(1: True) (4: True)]")


class TestSparseArray(unittest.TestCase):
    def setUp(self):
        aList = SparseList({1: True, 4: True}, length=6)
        cList = SparseList({1: 0.2, 3: 0.5}, default=0, length=6)
        list_1 = [5, 4, 3, 4, 5, 6]
        self.aMatrix = SparseArray(list_1=list_1, aList=aList, cList=cList)

    def test_len(self):
        self.assertEqual(len(self.aMatrix), 6)

    def test_add_tag(self):
        self.aMatrix.add_tag("coordinates")
        self.assertTrue("coordinates" in self.aMatrix.keys())
        self.assertTrue(len(self.aMatrix.coordinates) == 6)
        self.assertTrue(isinstance(self.aMatrix.coordinates, SparseList))

        self.aMatrix.add_tag(rel=[True, True, True])
        self.assertTrue(self.aMatrix.rel[1] == [True, True, True])

        # self.assertRaises(ValueError, self.aMatrix.add_tag, SparseList([], length=5))

    def test__getitem__(self):
        self.assertEqual(len(self.aMatrix[2:4]), 2)
        self.assertEqual(self.aMatrix[2].list_1, 3)

    def test__add__(self):
        bMatrix = self.aMatrix[2:4]
        bMatrix.add_tag(rel=True)
        bMatrix += self.aMatrix
        self.assertEqual(len(bMatrix), 8)
        self.assertEqual(bMatrix.list_1, [3, 4, 5, 4, 3, 4, 5, 6])

        # TODO: introduce SparseElement and test sum

    def test__mul__(self):
        fac = 100
        bMatrix = self.aMatrix * fac
        self.assertEqual(len(bMatrix), 6 * fac)
        self.assertEqual(bMatrix.list_1, fac * [5, 4, 3, 4, 5, 6])

        bMatrix = fac * self.aMatrix
        self.assertEqual(len(bMatrix), 6 * fac)
        self.assertEqual(bMatrix.list_1, fac * [5, 4, 3, 4, 5, 6])

    def test__getattr__(self):
        self.aMatrix.add_tag("spin")
        self.assertEqual(len(self.aMatrix.spin), 6)
        self.aMatrix.spin[1::2] = True
        self.aMatrix.spin[0::2] = False
        self.assertEqual(
            self.aMatrix.spin.list(), [False, True, False, True, False, True]
        )

    def test_keys(self):
        self.assertEqual(
            sorted(list(self.aMatrix.keys())), ["aList", "cList", "list_1"]
        )

    def test_items(self):
        for key, val in self.aMatrix[2:4].items():
            if key == "list_1":
                self.assertEqual(val, [3, 4])


if __name__ == "__main__":
    unittest.main()
