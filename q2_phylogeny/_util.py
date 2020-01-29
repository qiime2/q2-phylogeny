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
