# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import numpy as np
import pandas as pd
import qiime2
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


def filter_tree(tree: skbio.TreeNode,
                table: biom.Table = None,
                metadata: qiime2.Metadata = None,
                where: str = None,
                ) -> skbio.TreeNode:
    """
    Prunes a phylogenetic tree to match the input ids
    """
    # Checks the input metadata
    if ((table is None) & (metadata is None)):
        raise ValueError('A feature table, sequences or metadata must be '
                         'provided for filtering.')
    filter_refs = [table, metadata]
    if np.sum([(ref is not None) for ref in filter_refs]).sum() > 1:
        raise ValueError('Filtering can only be performed using one reference'
                         ' file. Please choose between filtering with a '
                         'feature table, sequences, or metadata.')
    if (where is not None) & (metadata is None):
        raise ValueError("Metadata must be provided if 'where' is specified")

    # Gets the list of IDs to keep
    if table is not None:
        ids_to_keep = table.ids(axis='observation')
    if metadata is not None:
        ids_to_keep = metadata.get_ids(where)

    # Gets the list of tips
    tip_ids = set([t.name for t in tree.tips()])

    # Checks for an intersection between ids
    if not set(tip_ids).issuperset(set(ids_to_keep)):
        raise ValueError('The ids for filtering must be a subset of '
                         'the tips in the tree.')
    sub_tree = tree.shear(ids_to_keep)
    sub_tree.prune()

    return sub_tree
