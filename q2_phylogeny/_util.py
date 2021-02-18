# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio


def midpoint_root(tree: skbio.TreeNode) -> skbio.TreeNode:
    return tree.root_at_midpoint()


def robinson_foulds(trees: skbio.TreeNode, labels: str = None,
                    missing_tips: str = 'error') -> skbio.DistanceMatrix:
    if labels is None:
        labels = ['tree_%d' % d for d in range(1, len(trees) + 1)]
    elif len(trees) != len(labels):
        raise ValueError("The number of trees and labels must match.")

    tips = [{t.name for t in tree.tips()} for tree in trees]
    shared_tips = set.intersection(*tips)
    if not shared_tips:
        raise ValueError("No tip names are shared between these trees.")
    if missing_tips == 'intersect-all':
        trees = [t.shear(shared_tips) for t in trees]
    elif missing_tips == 'error':
        all_tips = set.union(*tips)
        if shared_tips != all_tips:
            SKIP = 10
            miss = list(all_tips - shared_tips)
            missing_repr = ", ".join(map(repr, miss[:SKIP]))
            if len(miss) > SKIP:
                missing_repr += ", ...<%d ommitted>" % (len(miss) - SKIP)

            raise ValueError("Not all tips are shared between trees: "
                             + missing_repr)
    else:
        raise ValueError("Unknown argument for missing_tips=%r"
                         % missing_tips)

    return skbio.DistanceMatrix.from_iterable(
        trees, metric=skbio.TreeNode.compare_rfd,
        keys=labels, validate=False)
