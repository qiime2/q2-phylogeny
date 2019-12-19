# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def align_to_tree_mafft_raxml(ctx, sequences, n_threads=1,
                              mask_max_gap_frequency=1.0,
                              mask_min_conservation=0.40,
                              parttree=False,
                              substitution_model='GTRGAMMA',
                              seed=None, raxml_version='Standard'):
    mafft = ctx.get_action('alignment', 'mafft')
    mask = ctx.get_action('alignment', 'mask')
    raxml = ctx.get_action('phylogeny', 'raxml')
    midpoint_root = ctx.get_action('phylogeny', 'midpoint_root')

    aligned_seq, = mafft(sequences=sequences, n_threads=n_threads,
                         parttree=parttree)
    masked_seq, = mask(alignment=aligned_seq,
                       max_gap_frequency=mask_max_gap_frequency,
                       min_conservation=mask_min_conservation)
    unrooted_tree, = raxml(alignment=masked_seq, n_threads=n_threads,
                           substitution_model=substitution_model,
                           seed=seed, raxml_version=raxml_version)
    rooted_tree, = midpoint_root(tree=unrooted_tree)

    return (aligned_seq, masked_seq, unrooted_tree, rooted_tree)
