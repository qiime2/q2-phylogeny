# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io

import numpy as np
import pandas as pd
from biom.table import Table
import skbio

from qiime2 import Metadata

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
    def setUp(self):
        rooted_nwk = io.StringIO("((A:0.1, B:0.2)C:0.3, D:0.4, E:0.5)root;")
        self.tree = skbio.TreeNode.read(rooted_nwk)
        self.metadata = Metadata(pd.DataFrame(
                data=np.array([['Bacteria', '1'],
                               ['Archea', '1']], dtype=object),
                index=pd.Index(['A', 'D'], name='Feature ID'),
                columns=['kingdom', 'keep'],
            ))
        self.table = Table(data=np.array([[0, 1, 2], [2, 2, 2]]),
                           observation_ids=['A', 'D'],
                           sample_ids=['S1', 'S2', 'S3']
                           )
        self.filtered_tree = self.tree.copy().shear(['A', 'D'])
        self.filtered_tree.prune()

    def test_filter_tree_error_no_filter_art(self):
        with self.assertRaises(ValueError) as err:
            filter_tree(self.tree)
        self.assertEqual(
            str(err.exception),
            ('A feature table, sequences or metadata must be provided for '
                'filtering.')
            )

    def test_filter_tree_error_multiple_filter_arts(self):
        with self.assertRaises(ValueError) as err:
            filter_tree(self.tree,
                        table=self.table,
                        metadata=self.metadata)
        self.assertEqual(
            str(err.exception),
            ('Filtering can only be performed using one reference file. '
             'Please choose between filtering with a feature table, '
             'sequences, or metadata.')
        )

    def test_filter_tree_error_where_no_metadata(self):
        with self.assertRaises(ValueError) as err:
            filter_tree(self.tree,
                        table=self.table,
                        where='[kingdom]="Archea"')
        self.assertEqual(
            str(err.exception),
            ("Metadata must be provided if 'where' is specified")
            )

    def test_filter_tree_error_filter_superset(self):
        metadata = Metadata(pd.DataFrame(
            data=np.array([[1, 1, 0]]).T,
            index=pd.Index(['A', 'D', 'F'], name='feature-id'),
            columns=['keep']
        ))
        with self.assertRaises(ValueError) as err:
            filter_tree(self.tree,
                        metadata=metadata)
        self.assertEqual(
            str(err.exception),
            'The ids for filtering must be a subset of the tips in the tree.'
            )

    def test_filter_tree_metadata(self):
        test = filter_tree(self.tree, metadata=self.metadata)
        self.assertEqual(str(test), str(self.filtered_tree))

    def test_filter_tree_table(self):
        test = filter_tree(self.tree, table=self.table)
        self.assertEqual(str(test), str(self.filtered_tree))


if __name__ == "__main__":
    unittest.main()
