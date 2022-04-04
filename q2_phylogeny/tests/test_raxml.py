# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import shutil
import unittest
import skbio
import tempfile

from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio
from q2_types.feature_data import AlignedDNAFASTAFormat

from q2_phylogeny import raxml, raxml_rapid_bootstrap
from q2_phylogeny._raxml import (run_command, _build_rapid_bootstrap_command,
                                 _set_raxml_version)


class RaxmlTests(TestPluginBase):

    package = 'q2_phylogeny.tests'

    @classmethod
    def setUpClass(cls):
        super(TestPluginBase, cls).setUpClass()
        tmpdir = tempfile.mkdtemp()
        src = pkg_resources.resource_filename(cls.package, 'data')
        dst = os.path.join(tmpdir, 'data')
        shutil.copytree(src, dst)
        cls.data_dir = dst

    @classmethod
    def tearDownClass(cls):
        super(TestPluginBase, cls).setUpClass()
        shutil.rmtree(cls.data_dir)

    def get_data_path(self, filename):
        # Override TestPluginBase.get_data_path so that it returns paths to
        # temporary copies of test data.
        return os.path.join(self.data_dir, filename)

    def test_raxml(self):
        # Test that output tree is made.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences)
        obs_tree = skbio.TreeNode.read(str(obs))
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA000686145',
                              'GCA001950115', 'GCA001971985', 'GCA900007555']))

    def test_raxml_underscore_ids(self):
        # Test that output tree is made with underscores in tip IDs.
        # Some programs and python wrappers may strip underscores.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-4.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences)
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

    def test_set_raxml_version(self):
        obs_stand_1 = _set_raxml_version(raxml_version='Standard',
                                         n_threads=1)
        self.assertTrue('raxmlHPC' in str(obs_stand_1[0]))
        self.assertTrue(len(obs_stand_1) == 1)

        obs_sse3_1 = _set_raxml_version(raxml_version='SSE3', n_threads=1)
        self.assertTrue('raxmlHPC-SSE3' in str(obs_sse3_1[0]))
        self.assertTrue(len(obs_sse3_1) == 1)

        obs_avx2_1 = _set_raxml_version(raxml_version='AVX2', n_threads=1)
        self.assertTrue('raxmlHPC-AVX2' in str(obs_avx2_1[0]))
        self.assertTrue(len(obs_avx2_1) == 1)

        obs_stand_4 = _set_raxml_version(raxml_version='Standard',
                                         n_threads=4)
        self.assertTrue('raxmlHPC-PTHREADS' in str(obs_stand_4[0]))
        self.assertTrue('4' in str(obs_stand_4[1]))
        self.assertTrue(len(obs_stand_4) == 2)

        obs_sse3_4 = _set_raxml_version(raxml_version='SSE3', n_threads=4)
        self.assertTrue('raxmlHPC-PTHREADS-SSE3' in str(obs_sse3_4[0]))
        self.assertTrue('4' in str(obs_sse3_4[1]))
        self.assertTrue(len(obs_sse3_4) == 2)

        obs_avx2_4 = _set_raxml_version(raxml_version='AVX2', n_threads=4)
        self.assertTrue('raxmlHPC-PTHREADS-AVX2' in str(obs_avx2_4[0]))
        self.assertTrue('4' in str(obs_avx2_4[1]))
        self.assertTrue(len(obs_avx2_4) == 2)

    def test_raxml_version(self):
        # Test that an output tree is made when invoking threads.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, raxml_version='SSE3')
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)

        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]

        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA000686145',
                              'GCA001950115', 'GCA001971985', 'GCA900007555']))

    def test_raxml_n_threads(self):
        # Test that an output tree is made when invoking threads.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, n_threads=2)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)

        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]

        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA000686145',
                              'GCA001950115', 'GCA001971985', 'GCA900007555']))

    def test_raxml_with_seed(self):
        # Test tip-to-tip dists are identical to manually run RAxML output.
        # This test is comparing an ordered series of tip-to-tip distances
        # to a tree output from a manual run of the default command:
        # raxmlHPC -m GTRGAMMA -p 1723 -s aligned-dna-sequences-3.fasta -n q2
        # NOTE: I cleanly rounded the tip-to-tip dists (i.e. `%.4f`) as RAxML
        # may return slightly different rounding errors on different
        # systems.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, seed=1723)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        obs_tl = list(obs_tree.tip_tip_distances().to_series())
        obs_series = set(['%.4f' % e for e in obs_tl])

        exp_tree = skbio.TreeNode.read(self.get_data_path('test.tre'))
        exp_tl = list(exp_tree.tip_tip_distances().to_series())
        exp_series = set(['%.4f' % e for e in exp_tl])

        self.assertEqual(obs_series, exp_series)

    def test_raxml_model_choice(self):
        # Tip to tip dists should NOT be identical under different models.
        # Default is GTRGAMMA, we'll compare ouput to GRTGAMMAI & GTRCAT.
        # This test is comparing an ordered series of tip-to-tip distances.
        # Take note, that for this comparison to work, all must have the same
        # seed value set.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        # default GTRGAMMA
        with redirected_stdio(stderr=os.devnull):
            gtrg = raxml(input_sequences, seed=1723)
            gtrg_tree = skbio.TreeNode.read(
                        str(gtrg), convert_underscores=False)
            gtrg_td = set(gtrg_tree.tip_tip_distances().to_series())

        # set GTRGAMMAI
        with redirected_stdio(stderr=os.devnull):
            gtrgi = raxml(input_sequences, seed=1723,
                          substitution_model='GTRGAMMAI')
            gtrgi_tree = skbio.TreeNode.read(
                         str(gtrgi), convert_underscores=False)
            gtrgi_td = set(gtrgi_tree.tip_tip_distances().to_series())

        # set GTRCAT
        with redirected_stdio(stderr=os.devnull):
            gtrcat = raxml(input_sequences, seed=1723,
                           substitution_model='GTRCAT')
            gtrcat_tree = skbio.TreeNode.read(
                          str(gtrcat), convert_underscores=False)
            gtrcat_td = set(gtrcat_tree.tip_tip_distances().to_series())

        # test pairs are not equivalent
        self.assertNotEqual(gtrg_td, gtrgi_td)
        self.assertNotEqual(gtrg_td, gtrcat_td)
        self.assertNotEqual(gtrgi_td, gtrcat_td)

    def test_raxml_num_searches(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = raxml(input_sequences, seed=1723, n_searches=5)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        obs_tl = list(obs_tree.tip_tip_distances().to_series())
        obs_series = set(['%.4f' % e for e in obs_tl])

        exp_tree = skbio.TreeNode.read(self.get_data_path('test3.tre'))
        exp_tl = list(exp_tree.tip_tip_distances().to_series())
        exp_series = set(['%.4f' % e for e in exp_tl])
        self.assertEqual(obs_series, exp_series)

    def test_rapid_bootstrap_command(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with tempfile.TemporaryDirectory() as temp_dir:
            with redirected_stdio(stderr=os.devnull):
                obs = _build_rapid_bootstrap_command(input_sequences, 1723,
                                                     8752, 15, 'GTRGAMMA',
                                                     temp_dir, 'bs')
        self.assertTrue(str(input_sequences) in str(obs[11]))
        self.assertTrue('1723' in obs[5])
        self.assertTrue('8752' in obs[7])
        self.assertTrue('15' in obs[9])
        self.assertTrue('GTRGAMMA' in obs[3])
        self.assertTrue(str(temp_dir) in obs[13])
        self.assertTrue('bs' in obs[15])

    def test_raxml_rapid_bootstrap(self):
        # Test that output tree is made.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = raxml_rapid_bootstrap(input_sequences)
        obs_tree = skbio.TreeNode.read(str(obs))
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA000686145',
                              'GCA001950115', 'GCA001971985', 'GCA900007555']))

    def test_raxml_rapid_bootstrap_n_threads(self):
        # Test that an output tree is made when invoking threads.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = raxml_rapid_bootstrap(input_sequences, n_threads=2)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)

        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]

        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA000686145',
                              'GCA001950115', 'GCA001971985', 'GCA900007555']))

    def test_raxml_rapid_bootstrap_with_seed(self):
        # Test tip-to-tip dists are identical to manually run RAxML output.
        # This test is comparing an ordered series of tip-to-tip distances
        # to a tree output from a manual run of the default command:
        #     raxmlHPC -f a -m GTRGAMMA -p 1723 -x 3871 -N 10
        #         -s aligned-dna-sequences-3.fasta -n q2
        # NOTE: I cleanly rounded the tip-to-tip dists (i.e. `%.4f`) as RAxML
        # may return slightly different rounding errors on different
        # systems (and at times, between conda environments).
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        # test that branchlengths are identical
        with redirected_stdio(stderr=os.devnull):
            obs = raxml_rapid_bootstrap(input_sequences, seed=1723,
                                        rapid_bootstrap_seed=3871,
                                        bootstrap_replicates=10)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        # sometimes we lose the last set of numbers on long floats
        obs_tl = list(obs_tree.tip_tip_distances().to_series())
        obs_series = set(['%.4f' % e for e in obs_tl])

        exp_tree = skbio.TreeNode.read(self.get_data_path('test2.tre'),
                                       convert_underscores=True)
        exp_tl = list(exp_tree.tip_tip_distances().to_series())
        exp_series = set(['%.4f' % e for e in exp_tl])
        self.assertEqual(obs_series, exp_series)

        # test that bootstrap supports are identical
        obs_bs = [node.name for node in obs_tree.non_tips()].sort()
        exp_bs = [node.name for node in exp_tree.non_tips()].sort()
        self.assertEqual(obs_bs, exp_bs)

    def test_run_not_verbose(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        aligned_fp = str(input_sequences)

        with tempfile.TemporaryDirectory() as temp_dir:
            cmd = ['raxmlHPC',
                   '-m', 'GTRGAMMA',
                   '-p', '1723',
                   '-s', aligned_fp,
                   '-w', temp_dir,
                   '-n', 'q2']

            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, verbose=False)
            obs_tree_fp = os.path.join(temp_dir, 'RAxML_bestTree.q2')
            obs_tree = skbio.TreeNode.read(str(obs_tree_fp),
                                           convert_underscores=False)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515',
                              'GCA000454205', 'GCA000473545',
                              'GCA000196255', 'GCA000686145',
                              'GCA001950115', 'GCA001971985',
                              'GCA900007555']))

    def test_run_rapid_bs_not_verbose(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        aligned_fp = str(input_sequences)

        with tempfile.TemporaryDirectory() as temp_dir:
            cmd = ['raxmlHPC',
                   '-m', 'GTRGAMMA',
                   '-p', '1723',
                   '-s', aligned_fp,
                   '-w', temp_dir,
                   '-n', 'q2',
                   '-f', 'a',
                   '-x', '9834',
                   '-N', '10']

            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, verbose=False)

            obs_tree_fp = os.path.join(temp_dir, 'RAxML_bipartitions.q2')
            obs_tree = skbio.TreeNode.read(str(obs_tree_fp),
                                           convert_underscores=False)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515',
                              'GCA000454205', 'GCA000473545',
                              'GCA000196255', 'GCA000686145',
                              'GCA001950115', 'GCA001971985',
                              'GCA900007555']))


if __name__ == "__main__":
    unittest.main()
