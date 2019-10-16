# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def align_to_tree_mafft_iqtree(ctx, sequences, n_threads=1,
                               mask_max_gap_frequency=1.0,
                               mask_min_conservation=0.40,
                               substitution_model='MFP',
                               fast=False, alrt=None,
                               seed=None):
    mafft = ctx.get_action('alignment', 'mafft')
    mask = ctx.get_action('alignment', 'mask')
    iqtree = ctx.get_action('phylogeny', 'iqtree')
    midpoint_root = ctx.get_action('phylogeny', 'midpoint_root')

    aligned_seq, = mafft(sequences=sequences, n_threads=n_threads)
    masked_seq, = mask(alignment=aligned_seq,
                       max_gap_frequency=mask_max_gap_frequency,
                       min_conservation=mask_min_conservation)
    unrooted_tree, = iqtree(alignment=masked_seq, n_cores=n_threads,
                            fast=fast, alrt=alrt, seed=seed,
                            substitution_model=substitution_model)
    rooted_tree, = midpoint_root(tree=unrooted_tree)

    return (aligned_seq, masked_seq, unrooted_tree, rooted_tree)
