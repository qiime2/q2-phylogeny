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
import subprocess

from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio
from q2_types.feature_data import AlignedDNAFASTAFormat
from q2_types.tree import NewickFormat

from q2_phylogeny import raxml
from q2_phylogeny._raxml import run_command


class RaxmlTests(TestPluginBase):

    package = 'q2_phylogeny.tests'


    def test_raxml(self):
        """Test that output tree is made.
        Reads tree output and compares tip labels to expected labels.
        """
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, seed=1723)
        obs_tree = skbio.TreeNode.read(str(obs))
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names), set(['GCA001510755','GCA001045515',
        'GCA000454205', 'GCA000473545','GCA000196255','GCA002142615',
        'GCA000686145','GCA001950115','GCA001971985','GCA900007555']))


    def test_raxml_underscore_ids(self):
        """Test that output tree is made with underscores in tip IDs.
        Some programs and python wrappers may strip underscores.
        Reads tree output and compares tip labels to expected labels.
        """
        input_fp = self.get_data_path('aligned-dna-sequences-4.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, seed=1723)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names), set(['GCA_001510755_1',
        'GCA_001045515_1', 'GCA_000454205_1', 'GCA_000473545_1',
        'GCA_000196255_1', 'GCA_002142615_1', 'GCA_000686145_1',
        'GCA_001950115_1','GCA_001971985_1', 'GCA_900007555_1']))


    def test_raxml_n_threads(self):
        """Test that an output tree is made when invoking threads."""
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, seed=1723, n_threads=2)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names), set(['GCA001510755','GCA001045515',
	    'GCA000454205', 'GCA000473545','GCA000196255','GCA002142615',
	    'GCA000686145', 'GCA001950115','GCA001971985','GCA900007555']))


    def test_raxml_with_seed(self):
        """Test tip-to-tip dists are identical to manually run RAxML output.
        This test is comparing an ordered series of tip-to-tip distances
        to a tree output from a manual run of the default command:
            raxmlHPC -m GTRGAMMA -p 1723 -s aligned-dna-sequences-3.fasta
                -w tmp_dir -n q2
        """
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, seed=1723)
            obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
            obs_series = set(obs_tree.tip_tip_distances().to_series())

        exp_tree = skbio.TreeNode.read(self.get_data_path('test.tre'))
        exp_series = set(exp_tree.tip_tip_distances().to_series())

        self.assertEqual(obs_series, exp_series)


    def test_raxml_model_choice(self):
        """Tip to tip dists should NOT be identical under different models.
        Default is GTRGAMMA, we'll compare ouput to GRTGAMMAI & GTRCAT.

        This test is comparing an ordered series of tip-to-tip distances.
        """
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        # default GTRGAMMA
        with redirected_stdio(stderr=os.devnull):
            gtrg = raxml(input_sequences, seed=1723)
            gtrg_tree = skbio.TreeNode.read(str(gtrg), 
	                convert_underscores=False)
            gtrg_td = set(gtrg_tree.tip_tip_distances().to_series())

        # set GTRGAMMAI
        with redirected_stdio(stderr=os.devnull):
            gtrgi = raxml(input_sequences, seed=1723, 
	            substitution_model='GTRGAMMAI')
            gtrgi_tree = skbio.TreeNode.read(str(gtrgi), 
	                    convert_underscores=False)
            gtrgi_td = set(gtrgi_tree.tip_tip_distances().to_series())

        # set GTRCAT
        with redirected_stdio(stderr=os.devnull):
            gtrcat = raxml(input_sequences, seed=1723, 
	                substitution_model='GTRCAT')
            gtrcat_tree = skbio.TreeNode.read(str(gtrcat), 
	                    convert_underscores=False)
            gtrcat_td = set(gtrcat_tree.tip_tip_distances().to_series())

        # test pairs are not equivalent
        self.assertNotEqual(gtrg_td, gtrgi_td)
        self.assertNotEqual(gtrg_td, gtrcat_td)
        self.assertNotEqual(gtrgi_td, gtrcat_td)


if __name__ == "__main__":
    unittest.main()
