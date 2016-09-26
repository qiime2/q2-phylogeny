# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import skbio

from q2_types.testing import TestPluginBase
from q2_types.feature_data import AlignedDNAFASTAFormat

from q2_phylogeny import fasttree


class FastTreeTests(TestPluginBase):

    package = 'q2_phylogeny.tests'

    def test_fasttree(self):
        input_fp = self.get_data_path('aligned-dna-sequences-1.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        obs = fasttree(input_sequences)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids (the branch lengths can vary with
        # different versions of FastTree)
        obs_tree = skbio.io.read(obs.open(), format='newick',
                                 into=skbio.TreeNode)
        tips = list(obs_tree.tips())
        self.assertEqual(len(tips), 2)
        tip_names = [t.name for t in tips]
        tip_names.sort()
        self.assertEqual(tip_names, ['seq1', 'seq2'])


if __name__ == "__main__":
    unittest.main()
