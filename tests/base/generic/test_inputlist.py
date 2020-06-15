from pyiron.base.generic.inputlist import InputList
from pyiron.base.generic.hdfio import ProjectHDFio
from pyiron.base.project.generic import Project
import os
import unittest
import json

class TestInputList(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pl = InputList([
                {'foo': 'bar'},
                2,
                42,
                {'next':
                    [0,
                        {'depth': 23}
                    ]
                }
        ], table_name = 'input')
        cls.pl['tail'] = InputList([2,4,8])

        file_location = os.path.dirname(os.path.abspath(__file__))
        pr = Project(file_location)
        cls.file_name = os.path.join(file_location, "input.h5")
        cls.hdf = ProjectHDFio(project=pr, file_name=cls.file_name,
                               h5_path="/test", mode="a")

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.file_name)

    ## Init tests
    def test_init_none(self):
        pl = InputList()
        self.assertEqual(len(pl), 0, 'not empty after initialized with None')

    def test_init_list(self):
        l = [1, 2, 3, 4]
        pl = InputList(l)
        self.assertEqual(len(pl), len(l),
                'not the same length as source list')
        self.assertEqual(list(pl.values()), l,
                'conversion to list not the same as source list')

    def test_init_tuple(self):
        t = (1, 2, 3, 4)
        pl = InputList(t)
        self.assertEqual(len(pl), len(t),
                'not the same length as source tuple')
        self.assertEqual(tuple(pl.values()), t,
                'conversion to tuple not the same as source tuple')

    def test_init_set(self):
        s = {1, 2, 3, 4}
        pl = InputList(s)
        self.assertEqual(len(pl), len(s),
                'not the same length as source set')
        self.assertEqual(set(pl.values()), s,
                'conversion to set not the same as source set')

    def test_init_dict(self):
        d = {'foo': 23, 'test case': 'bar'}
        pl = InputList(d)
        self.assertEqual(tuple(pl.items()), tuple(d.items()),
                'source dict items not preserved')
        with self.assertRaises(ValueError,
                               msg = 'no ValueError on invalid initializer'):
            pl = InputList({2: 0, 1: 1})


    # access tests
    def test_get_nested(self):
        n = [
                {'foo': 'bar'},
                2,
                42,
                {'next':
                    [0,
                        {'depth': 23}
                    ]
                }
        ]
        pl = InputList(n)
        self.assertEqual(type(pl[0]), InputList,
                'nested dict not converted to InputList')
        self.assertEqual(type(pl['3/next']), InputList,
                'nested list not converted to InputList')
        self.assertEqual(type(pl['0/foo']), str,
                'nested str converted to InputList')

    def test_get_item(self):
        self.assertEqual(self.pl[0], {'foo': 'bar'},
                'index with integer does not give correct element')
        self.assertEqual(self.pl[0]['foo'], 'bar',
                'index with string does not give correct element')

    def test_get_attr(self):
        self.assertEqual(self.pl.tail, InputList([2, 4, 8]),
                'attribute access does not give correct element')
        self.assertEqual(self.pl[3].next, InputList([0, InputList({'depth': 23})]),
                'nested attribute access does not give correct element')

    def test_get_sempath(self):
        self.assertEqual(self.pl['0'], {'foo': 'bar'},
                'decimal string not converted to integer')
        self.assertEqual(self.pl['0/foo'], 'bar',
                'nested access does not give correct element')
        self.assertEqual(self.pl['3/next/1/depth'], 23,
                'nested access does not give correct element')
        self.assertEqual(self.pl['3/next/0'], 0,
                'nested access does not give correct element')
        self.assertEqual(self.pl['3/next/1/depth'],
                         self.pl[3, 'next', 1, 'depth'],
                         'access via semantic path and tuple not the same')

    def test_get_string_int(self):
        self.assertEqual(self.pl[0], self.pl['0'],
                         'access via index and digit-only string not the same')

    def test_set_some_keys(self):
        pl = InputList([1,2])
        pl['end'] = 3
        self.assertEqual(pl, InputList({0: 1, 1: 2, 'end': 3}))

    def test_set_append(self):
        pl = InputList()
        # should not raise and exception
        pl[0] = 1
        pl[1] = 2
        self.assertEqual(pl[0], 1,
                         'append via index broken on empty list')
        self.assertEqual(pl[1], 2,
                         'append via index broken on non-empty list')

    def test_update(self):
        pl = InputList()
        d = self.pl.to_builtin()
        pl.update(d, wrap = True)
        self.assertEqual(pl, self.pl,
                         'update from to_builtin does not restore list')


    # hdf tests
    def test_to_hdf(self):
        self.pl.to_hdf(hdf=self.hdf)
        self.assertEqual(self.hdf["input/NAME"],
                         "InputList")
        self.assertEqual(self.hdf["input/OBJECT"],
                         "InputList")
        self.assertEqual(self.hdf["input/TYPE"],
                         "<class 'pyiron.base.generic.inputlist.InputList'>")
        l = InputList(self.hdf["input/data"])
        self.assertEqual(self.pl, l)

    def test_to_hdf_group(self):
        self.pl.to_hdf(hdf=self.hdf, group_name = "test_group")
        self.assertEqual(self.hdf["test_group/NAME"],
                         "InputList")
        self.assertEqual(self.hdf["test_group/TYPE"],
                         "<class 'pyiron.base.generic.inputlist.InputList'>")
        self.assertEqual(self.hdf["test_group/OBJECT"],
                         "InputList")
        l = InputList(self.hdf["test_group/data"])
        self.assertEqual(self.pl, l)

    def test_from_hdf(self):
        self.pl.to_hdf(hdf=self.hdf)
        l = InputList(table_name = "input")
        l.from_hdf(hdf=self.hdf)
        self.assertEqual(self.pl, l)

    def test_from_hdf_group(self):
        self.pl.to_hdf(hdf=self.hdf, group_name = "test_group")
        l = InputList(table_name = "input")
        l.from_hdf(hdf=self.hdf, group_name = "test_group")
        self.assertEqual(self.pl, l)


if __name__ == '__main__':
    unittest.main()
