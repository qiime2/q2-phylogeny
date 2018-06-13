# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest
import skbio
import tempfile

from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio
from q2_types.feature_data import AlignedDNAFASTAFormat

from q2_phylogeny import iqtree
from q2_phylogeny._raxml import run_command


class IqtreeTests(TestPluginBase):

    package = 'q2_phylogeny.tests'

    def test_iqtree(self):
        # Test that output tree is made.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = iqtree(input_sequences)
        obs_tree = skbio.TreeNode.read(str(obs))
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA002142615',
                              'GCA000686145', 'GCA001950115', 'GCA001971985',
                              'GCA900007555']))

    def test_iqtree_underscore_ids(self):
        # Test that output tree is made with underscores in tip IDs.
        # Some programs and python wrappers may strip underscores.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-4.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = iqtree(input_sequences)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA_001510755_1', 'GCA_001045515_1',
                              'GCA_000454205_1', 'GCA_000473545_1',
                              'GCA_000196255_1', 'GCA_002142615_1',
                              'GCA_000686145_1', 'GCA_001950115_1',
                              'GCA_001971985_1', 'GCA_900007555_1']))



if __name__ == "__main__":
    unittest.main()
