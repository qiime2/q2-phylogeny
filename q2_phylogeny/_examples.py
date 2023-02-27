# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


rep_seqs_url = ('https://data.qiime2.org/usage-examples/'
                'moving-pictures/rep-seqs-dada2.qza')


def phylogeny_align_to_tree_mafft_fasttree(use):
    rep_seqs = use.init_artifact_from_url('rep_seqs', rep_seqs_url)

    alignment, masked_alignment, tree, rooted_tree = use.action(
        use.UsageAction('phylogeny', 'align_to_tree_mafft_fasttree'),
        use.UsageInputs(sequences=rep_seqs),
        use.UsageOutputNames(
            alignment='aligned-rep-seqs',
            masked_alignment='masked-aligned-rep-seqs',
            tree='unrooted-tree',
            rooted_tree='rooted-tree',
        )
    )

    alignment.assert_output_type('FeatureData[AlignedSequence]')
    masked_alignment.assert_output_type('FeatureData[AlignedSequence]')
    tree.assert_output_type('Phylogeny[Unrooted]')
    rooted_tree.assert_output_type('Phylogeny[Rooted]')
