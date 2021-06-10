# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
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

    def test_fasttree_n_tips(self):
        input_fp = self.get_data_path('aligned-dna-sequences-1.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = fasttree(input_sequences, n_threads='auto')
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids (the branch lengths can vary with
        # different versions of FastTree, and threading can produce
        # non-deterministic trees)
        obs_tree = skbio.TreeNode.read(str(obs))
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        tip_names.sort()
        self.assertEqual(tip_names, ['seq1', 'seq2'])


# We are using pytest vs. unittest for this test
# The capfd arg below is a built-in pytest fixture that allows for
# captured stderr from all subprocesses that directly write
# to operating system level output
def test_fasttree_num_threads(capfd):
    # Rather than providing an actual filepath, we just run a 'help' command
    # Ideally, we supply --help in this param, but not available in FastTree
    # So we are using -expert instead, so that FastTree doesn't fail
    # ######################
    # Output preview/example:
    # ######################
    # Detailed usage for FastTree 2.1.10 Double precision (No SSE3):
    # FastTree [-nt] [-n 100] [-quote] [-pseudo | -pseudo 1.0]
    #            [-boot 1000 | -nosupport]
    #            [-intree starting_trees_file | -intree1 starting_tree_file]
    #            [-quiet | -nopr]
    #            [-nni 10] [-spr 2] [-noml | -mllen | -mlnni 10]
    #            [-mlacc 2] [-cat 20 | -nocat] [-gamma]
    #            [-slow | -fastest] [-2nd | -no2nd] [-slownni] [-seed 1253]
    #            [-top | -notop] [-topm 1.0 [-close 0.75] [-refresh 0.8]]
    #            [-matrix Matrix | -nomatrix] [-nj | -bionj]
    #            [-lg] [-wag] [-nt] [-gtr] [-gtrrates ac ag at cg ct gt]
    #            [-gtrfreq A C G T]
    #            [ -constraints constraintAlignment
    #            [ -constraintWeight 100.0 ] ]
    #            [-log logfile]
    #          [ alignment_file ]
    #         [ -out output_newick_file | > newick_tree]
    fasttree('-expert', n_threads=1)
    captured = capfd.readouterr()
    assert 'FastTree' in captured.err

    fasttree('-expert', n_threads=20)
    captured = capfd.readouterr()
    assert 'OpenMP (20 threads)' in captured.err

    # This test case ensures that when a user enters 'auto', the n_threads
    # var will still be set to the max available on their machine, even if
    # the OMP_NUM_THREADS env var has been set on their machine
    os.environ['OMP_NUM_THREADS'] = '2560'
    fasttree('-expert', n_threads='auto')
    captured = capfd.readouterr()
    assert 'OpenMP (2560 threads)' not in captured.err


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
