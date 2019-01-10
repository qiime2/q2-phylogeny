# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io

import skbio

from q2_phylogeny import midpoint_root


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


if __name__ == "__main__":
    unittest.main()
