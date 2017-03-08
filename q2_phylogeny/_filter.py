# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import skbio


def filter_table(table: biom.Table, tree: skbio.TreeNode) -> biom.Table:
    """ Filter table to remove feature ids that are not tip ids in tree
    """
    tip_ids = set([t.name for t in tree.tips()])
    feature_ids = set(table.ids(axis='observation'))
    # ids_to_keep can only include ids that are in table
    ids_to_keep = tip_ids & feature_ids
    table.filter(ids_to_keep, axis='observation', inplace=True)
    return table
