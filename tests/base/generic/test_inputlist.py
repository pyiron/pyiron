from pyiron.base.generic.inputlist import InputList
import unittest

class InitTest(unittest.TestCase):

    def test_none(self):
        pl = InputList()
        self.assertEqual(len(pl), 0, 'not empty after initialized with None')

    def test_list(self):
        l = [1, 2, 3, 4]
        pl = InputList(l)
        self.assertEqual(len(pl), len(l),
                'not the same length as source list')
        self.assertEqual(list(pl.values()), l,
                'conversion to list not the same as source list')

    def test_tuple(self):
        t = (1, 2, 3, 4)
        pl = InputList(t)
        self.assertEqual(len(pl), len(t),
                'not the same length as source tuple')
        self.assertEqual(tuple(pl.values()), t,
                'conversion to tuple not the same as source tuple')

    def test_set(self):
        s = {1, 2, 3, 4}
        pl = InputList(s)
        self.assertEqual(len(pl), len(s),
                'not the same length as source set')
        self.assertEqual(set(pl.values()), s,
                'conversion to set not the same as source set')

    def test_dict(self):
        d = {'foo': 23, 'test case': 'bar'}
        pl = InputList(d)
        self.assertEqual(tuple(pl.items()), tuple(d.items()),
                'source dict items not preserved')

    def test_nested(self):
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

class GetterTest(unittest.TestCase):

    def setUp(self):
        self.pl = InputList([
                {'foo': 'bar'},
                2,
                42,
                {'next':
                    [0,
                        {'depth': 23}
                    ]
                }
        ])
        self.pl['tail'] = InputList([2,4,8])

    def test_item(self):
        self.assertEqual(self.pl[0], {'foo': 'bar'},
                'index with integer does not give correct element')
        self.assertEqual(self.pl[0]['foo'], 'bar',
                'index with string does not give correct element')

    def test_attr(self):
        self.assertEqual(self.pl.tail, InputList([2, 4, 8]),
                'attribute access does not give correct element')
        self.assertEqual(self.pl[3].next, InputList([0, InputList({'depth', 23})]),
                'nested attribute access does not give correct element')

    def test_sempath(self):
        self.assertEqual(self.pl['0'], {'foo': 'bar'},
                'decimal string not converted to integer')
        self.assertEqual(self.pl['0/foo'], 'bar',
                'nested access does not give correct element')
        self.assertEqual(self.pl['3/next/1/depth'], 23,
                'nested access does not give correct element')
        self.assertEqual(self.pl['3/next/0'], 0,
                'nested access does not give correct element')

if __name__ == '__main__':
    unittest.main()
