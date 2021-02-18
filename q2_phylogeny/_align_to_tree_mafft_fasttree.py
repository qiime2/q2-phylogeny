# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def align_to_tree_mafft_fasttree(ctx, sequences, n_threads=1,
                                 mask_max_gap_frequency=1.0,
                                 mask_min_conservation=0.40,
                                 parttree=False):
    mafft = ctx.get_action('alignment', 'mafft')
    mask = ctx.get_action('alignment', 'mask')
    fasttree = ctx.get_action('phylogeny', 'fasttree')
    midpoint_root = ctx.get_action('phylogeny', 'midpoint_root')

    aligned_seq, = mafft(sequences=sequences, n_threads=n_threads,
                         parttree=parttree)
    masked_seq, = mask(alignment=aligned_seq,
                       max_gap_frequency=mask_max_gap_frequency,
                       min_conservation=mask_min_conservation)
    unrooted_tree, = fasttree(alignment=masked_seq, n_threads=n_threads)
    rooted_tree, = midpoint_root(tree=unrooted_tree)

    return (aligned_seq, masked_seq, unrooted_tree, rooted_tree)
