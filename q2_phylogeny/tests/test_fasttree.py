# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
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

from q2_phylogeny import fasttree
from q2_phylogeny._fasttree import run_command


class FastTreeTests(TestPluginBase):

    package = 'q2_phylogeny.tests'

    def test_fasttree(self):
        input_fp = self.get_data_path('aligned-dna-sequences-1.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = fasttree(input_sequences)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids (the branch lengths can vary with
        # different versions of FastTree)
        obs_tree = skbio.TreeNode.read(str(obs))
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        tip_names.sort()
        self.assertEqual(tip_names, ['seq1', 'seq2'])

    def test_fasttree_underscore_ids(self):
        input_fp = self.get_data_path('aligned-dna-sequences-2.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = fasttree(input_sequences)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids (the branch lengths can vary with
        # different versions of FastTree)
        obs_tree = skbio.TreeNode.read(str(obs))
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        tip_names.sort()
        self.assertEqual(tip_names, ['_s_e_q_1_', '_s_e_q_2_'])


class RunCommandTests(TestPluginBase):

    package = 'q2_phylogeny.tests'

    def test_failed_run(self):
        input_fp = self.get_data_path('aligned-dna-sequences-1.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        result = NewickFormat()
        aligned_fp = str(input_sequences)
        tree_fp = str(result)

        cmd = ['FastTree', '-nt', '-not-a-real-parameter', aligned_fp]
        with self.assertRaises(subprocess.CalledProcessError):
            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, tree_fp)

    def test_failed_run_not_verbose(self):
        input_fp = self.get_data_path('aligned-dna-sequences-1.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        result = NewickFormat()
        aligned_fp = str(input_sequences)
        tree_fp = str(result)

        cmd = ['FastTree', '-nt', '-not-a-real-parameter', aligned_fp]
        with self.assertRaises(subprocess.CalledProcessError):
            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, tree_fp, verbose=False)


if __name__ == "__main__":
    unittest.main()
