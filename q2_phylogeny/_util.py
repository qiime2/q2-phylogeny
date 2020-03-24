# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio


def midpoint_root(tree: skbio.TreeNode) -> skbio.TreeNode:
    return tree.root_at_midpoint()


def robinson_foulds(trees: skbio.TreeNode, labels: str=None) \
        -> skbio.DistanceMatrix:
    if labels is None:
        labels = ['tree_%d' % d for d in range(1, len(trees) + 1)]
    elif len(trees) != len(labels):
        raise ValueError("The number of trees and labels must match.")

    return skbio.DistanceMatrix.from_iterable(
        trees, metric=lambda a, b: a.compare_rfd(b),
        keys=labels, validate=False)
