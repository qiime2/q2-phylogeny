# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact


class AlignToTreeMafftFasttreePipelineTest(TestPluginBase):
    package = 'q2_phylogeny.tests'

    def setUp(self):
        super().setUp()
        self.align_to_tree_mafft_iqtree = self.plugin.pipelines[
                                    'align_to_tree_mafft_iqtree']

        input_sequences_fp = self.get_data_path('dna-sequences-1.fasta')
        self.input_sequences = Artifact.import_data('FeatureData[Sequence]',
                                                    input_sequences_fp)

    def test_execution(self):
        # Does it run?
        self.align_to_tree_mafft_iqtree(self.input_sequences)
        self.assertTrue(True)

    def test_outputs(self):
        result = self.align_to_tree_mafft_iqtree(self.input_sequences)
        self.assertEqual(4, len(result))
        aligned_seq, masked_seq, unrooted_tree, rooted_tree = result
        self.assertEqual('FeatureData[AlignedSequence]', str(aligned_seq.type))
        self.assertEqual('FeatureData[AlignedSequence]', str(masked_seq.type)),
        self.assertEqual('Phylogeny[Unrooted]', str(unrooted_tree.type)),
        self.assertEqual('Phylogeny[Rooted]', str(rooted_tree.type))


if __name__ == '__main__':
    unittest.main()
