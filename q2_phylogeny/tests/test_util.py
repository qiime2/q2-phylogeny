# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io

import skbio

from q2_phylogeny import midpoint_root, robinson_foulds


class MidpointRootTests(unittest.TestCase):

    def test_midpoint_root(self):
        unrooted_nwk = io.StringIO("((a:1,b:1):1,(c:5,d:4):1);")
        rooted_nwk = io.StringIO("(c:4.5,(d:4,(a:1,b:1):2):0.5);")
        unrooted_tree = skbio.TreeNode.read(unrooted_nwk)
        expected = skbio.TreeNode.read(rooted_nwk)
        actual = midpoint_root(unrooted_tree)
        # skbio.TreeNode currently does not have an equality operator
        # (see https://github.com/biocore/scikit-bio/issues/465), so do a
        # few specific comparisons of actual and expected

        # confirm ids are as expected
        actual_tip_tip_distances = actual.tip_tip_distances()
        expected_tip_tip_distances = expected.tip_tip_distances()
        self.assertEqual(sorted(actual_tip_tip_distances.ids),
                         sorted(expected_tip_tip_distances.ids))

        # confirm tip-to-tip distances are as expected
        for id1 in actual_tip_tip_distances.ids:
            for id2 in actual_tip_tip_distances.ids:
                actual_distance = actual_tip_tip_distances[id1, id2]
                expected_distance = expected_tip_tip_distances[id1, id2]
                self.assertEqual(
                 actual_distance, expected_distance,
                 "Distances between %s and %s are not equal. (%f versus %f)." %
                 (id1, id2, actual_distance, expected_distance))

        # confirm tip-to-root distances are as expected
        for id_ in actual_tip_tip_distances.ids:
            self.assertEqual(actual.find(id_).distance(actual.root()),
                             expected.find(id_).distance(expected.root()))


class TestRobinsonFoulds(unittest.TestCase):
    def test_expected(self):
        # Using Felsenstein's test data set from here:
        # http://evolution.genetics.washington.edu/phylip/doc/treedist.html
        test_data_set_1 = [  # name from original source
            ['(A,(B,(H,(D,(J,(((G,E),(F,I)),C))))));'],
            ['(A,(B,(D,((J,H),(((G,E),(F,I)),C)))));'],
            ['(A,(B,(D,(H,(J,(((G,E),(F,I)),C))))));'],
            ['(A,(B,(E,(G,((F,I),((J,(H,D)),C))))));'],
            ['(A,(B,(E,(G,((F,I),(((J,H),D),C))))));'],
            ['(A,(B,(E,((F,I),(G,((J,(H,D)),C))))));'],
            ['(A,(B,(E,((F,I),(G,(((J,H),D),C))))));'],
            ['(A,(B,(E,((G,(F,I)),((J,(H,D)),C)))));'],
            ['(A,(B,(E,((G,(F,I)),(((J,H),D),C)))));'],
            ['(A,(B,(E,(G,((F,I),((J,(H,D)),C))))));'],
            ['(A,(B,(D,(H,(J,(((G,E),(F,I)),C))))));'],
            ['(A,(B,(E,((G,(F,I)),((J,(H,D)),C)))));']
        ]
        trees = [skbio.TreeNode.read(nwk) for nwk in test_data_set_1]
        expected = skbio.DistanceMatrix(
            [
                [0,   4,  2, 10, 10, 10, 10, 10, 10, 10,  2, 10],
                [4,   0,  2, 10,  8, 10,  8, 10,  8, 10,  2, 10],
                [2,   2,  0, 10, 10, 10, 10, 10, 10, 10,  0, 10],
                [10, 10, 10,  0,  2,  2,  4,  2,  4,  0, 10,  2],
                [10,  8, 10,  2,  0,  4,  2,  4,  2,  2, 10,  4],
                [10, 10, 10,  2,  4,  0,  2,  2,  4,  2, 10,  2],
                [10,  8, 10,  4,  2,  2,  0,  4,  2,  4, 10,  4],
                [10, 10, 10,  2,  4,  2,  4,  0,  2,  2, 10,  0],
                [10,  8, 10,  4,  2,  4,  2,  2,  0,  4, 10,  2],
                [10, 10, 10,  0,  2,  2,  4,  2,  4,  0, 10,  2],
                [2,   2,  0, 10, 10, 10, 10, 10, 10, 10,  0, 10],
                [10, 10, 10,  2,  4,  2,  4,  0,  2,  2, 10,  0]
            ], ids=['tree_%d' % d for d in range(1, 13)])

        result = robinson_foulds(trees)

        self.assertEqual(result, expected)

    def test_single_tree_and_label(self):
        trees = [skbio.TreeNode.read(['(A:0.2, B:1.5, C, (E, F));'])]
        expected = skbio.DistanceMatrix([[0]], ids=['foo'])

        result = robinson_foulds(trees, labels=['foo'])

        self.assertEqual(result, expected)

    def test_labels_too_short(self):
        test_data = [
            ['(A,(B,(H,(D,(J,(((G,E),(F,I)),C))))));'],
            ['(A,(B,(D,((J,H),(((G,E),(F,I)),C)))));'],
        ]
        trees = [skbio.TreeNode.read(nwk) for nwk in test_data]

        with self.assertRaisesRegex(ValueError, 'number.*match'):
            robinson_foulds(trees, labels=['A'])

    def test_missing_single_tip(self):
        test_data = [
            ['(A,(B,(H,(D,(J,(((G,E),(F,I)),C))))));'],
            ['((B,(D,((J,H),(((G,E),(F,I)),C)))));'],
        ]
        trees = [skbio.TreeNode.read(nwk) for nwk in test_data]

        with self.assertRaisesRegex(ValueError, "tips.*shared.*'A'"):
            robinson_foulds(trees)

    def test_missing_many_tips(self):
        test_data = [
            ['(A,(B,(H,(D,(J,(((G,E),(F,I,K,L)),C))))));'],
            ['(A,(,((,),(((,),(,)),))));'],
        ]
        trees = [skbio.TreeNode.read(nwk) for nwk in test_data]

        with self.assertRaisesRegex(ValueError, "tips.*shared"):
            robinson_foulds(trees)

    def test_missing_all_tips(self):
        test_data = [
            ['(A,(B,(H,(D,(J,(((G,E),(F,I,K,L)),C))))));'],
            ['(,(,((,),(((,),(,)),))));'],
        ]
        trees = [skbio.TreeNode.read(nwk) for nwk in test_data]

        with self.assertRaisesRegex(ValueError, "No tip names.*shared"):
            robinson_foulds(trees)

    def test_missing_intersect_all(self):
        test_data = [
            ['(A,(B,(H,(J,(((G,E),(F,I)),C)))));'],
            ['(B,((J,H),(((G,E),(F,I)),C)));'],
            ['(B,(D,(H,(J,(((G,E),(F,I)),C)))));'],
        ]
        trees = [skbio.TreeNode.read(nwk) for nwk in test_data]
        expected = skbio.DistanceMatrix([[0, 2, 0], [2, 0, 2], [0, 2, 0]])

        result = robinson_foulds(trees, labels=['0', '1', '2'],
                                 missing_tips='intersect-all')

        self.assertEqual(result, expected)

    def test_invalid_missing(self):
        test_data = [
            ['(A,(B,(H,(D,(J,(((G,E),(F,I,K,L)),C))))));'],
        ]
        trees = [skbio.TreeNode.read(nwk) for nwk in test_data]

        with self.assertRaisesRegex(ValueError, "not-an-option"):
            robinson_foulds(trees, missing_tips="not-an-option")


if __name__ == "__main__":
    unittest.main()
