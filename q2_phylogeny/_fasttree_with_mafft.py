# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def fasttree_with_mafft(ctx, sequences, n_threads=1, max_gap_frequency=1.0,
                        min_conservation=0.40):
    mafft = ctx.get_action('alignment', 'mafft')
    mask = ctx.get_action('alignment', 'mask')
    fasttree = ctx.get_action('phylogeny', 'fasttree')
    midpoint_root = ctx.get_action('phylogeny', 'midpoint_root')

    aligned_seq, = mafft(sequences=sequences, n_threads=n_threads)
    masked_seq, = mask(alignment=aligned_seq,
                       max_gap_frequency=max_gap_frequency,
                       min_conservation=min_conservation)
    unrooted_tree, = fasttree(alignment=masked_seq, n_threads=n_threads)
    rooted_tree, = midpoint_root(tree=unrooted_tree)

    return (aligned_seq, masked_seq, unrooted_tree, rooted_tree)
