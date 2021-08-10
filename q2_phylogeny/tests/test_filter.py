# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io

import numpy as np
from biom.table import Table
import skbio

from q2_phylogeny import (filter_table, filter_tree)


class FilterTableTests(unittest.TestCase):

    def test_tree_filter_table_some(self):
        rooted_nwk = io.StringIO("(O1:4.5,(d:4,(a:1,b:1):2):0.5);")
        tree = skbio.TreeNode.read(rooted_nwk)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_table(table, tree)
        expected = table.filter(['O1'], axis='observation')
        self.assertEqual(actual, expected)

    def test_tree_filter_table_all(self):
        rooted_nwk = io.StringIO("(c:4.5,(d:4,(a:1,b:1):2):0.5);")
        tree = skbio.TreeNode.read(rooted_nwk)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_table(table, tree)
        expected = table.filter([], axis='observation')
        self.assertEqual(actual, expected)

    def test_tree_filter_table_none(self):
        rooted_nwk = io.StringIO("(O1:4.5,(O2:4,(a:1,b:1):2):0.5);")
        tree = skbio.TreeNode.read(rooted_nwk)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_table(table, tree)
        expected = table.filter(['O1', 'O2'], axis='observation')
        self.assertEqual(actual, expected)


class FilterTreeTests(unittest.TestCase):

    def test_filter_tree_some_intersection(self):
        rooted_nwk = io.StringIO("(O1:4.5,(d:4,(a:1,b:1):2):0.5);")
        tree = skbio.TreeNode.read(rooted_nwk)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_tree(table, tree)
        expected = tree.shear(['O1']).prune()
        self.assertEqual(actual, expected)

    def test_filter_tree_no_intersection(self):
        rooted_nwk = io.StringIO("(c:4.5,(d:4,(a:1,b:1):2):0.5);")
        tree = skbio.TreeNode.read(rooted_nwk)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            filter_tree(table, tree)

    def test_filter_tree_perfect_intersection(self):
        rooted_nwk = io.StringIO("(O1:4.5,(O2:4,(a:1,b:1):2):0.5);")
        tree = skbio.TreeNode.read(rooted_nwk)
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_tree(table, tree)
        expected = tree.shear(['O1', 'O2']).prune()
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
