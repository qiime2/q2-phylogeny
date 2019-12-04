# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
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

from q2_phylogeny import iqtree, iqtree_ultrafast_bootstrap
from q2_phylogeny._raxml import run_command
from q2_phylogeny._iqtree import (_build_iqtree_command,
                                  _build_iqtree_ufbs_command)


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

    def test_iqtree_safe_allnni(self):
        # Same as `test_iqtree` but testing the `-safe` and `-allnni `flags
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = iqtree(input_sequences, safe='True', allnni='True')
        obs_tree = skbio.TreeNode.read(str(obs))
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

    def test_iqtree_n_cores(self):
        # Test that an output tree is made when invoking multiple cores.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = iqtree(input_sequences, n_cores=2)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)

        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]

        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA002142615',
                              'GCA000686145', 'GCA001950115', 'GCA001971985',
                              'GCA900007555']))

    def test_iqtree_n_cores_auto(self):
        # Test that an output tree is made when invoking automatic threads.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = iqtree(input_sequences, n_cores='auto')
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)

        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]

        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515', 'GCA000454205',
                              'GCA000473545', 'GCA000196255', 'GCA002142615',
                              'GCA000686145', 'GCA001950115', 'GCA001971985',
                              'GCA900007555']))

    def test_iqtree_with_seed(self):
        # Test tip-to-tip dists are identical to manually run IQ-TREE output.
        # This test is comparing an ordered series of tip-to-tip distances
        # to a tree output from a manual run of the default command:
        #     iqtree -seed 1723 -m HKY -s aligned-dna-sequences-3.fasta
        #            -nt 1 -pre q2iqtree
        # NOTE: I cleanly rounded the tip-to-tip dists (i.e. `%.4f`) as
        # IQ-TREE may return slightly different rounding errors on different
        # systems.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = iqtree(input_sequences, seed=1723,
                         substitution_model='HKY')
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        obs_tl = list(obs_tree.tip_tip_distances().to_series())
        obs_series = set(['%.4f' % e for e in obs_tl])

        exp_tree = skbio.TreeNode.read(self.get_data_path('test4.tre'))
        exp_tl = list(exp_tree.tip_tip_distances().to_series())
        exp_series = set(['%.4f' % e for e in exp_tl])

        self.assertEqual(obs_series, exp_series)

    def test_iqtree_model_choice(self):
        # Tip to tip dists should NOT be identical under different models.
        # Default is MFP (auto select substitution model). We'll compare ouput
        # of the GTR+G and HKY models.
        # This test is comparing an ordered series of tip-to-tip distances.
        # Take note, that for this comparison to work, all must have the same
        # seed value set.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        # default GTR+G
        with redirected_stdio(stderr=os.devnull):
            gtrg = iqtree(input_sequences, seed=1723,
                          substitution_model='GTR+G')
            gtrg_tree = skbio.TreeNode.read(
                        str(gtrg), convert_underscores=False)
            gtrg_td = set(gtrg_tree.tip_tip_distances().to_series())

        # set HKY
        with redirected_stdio(stderr=os.devnull):
            hky = iqtree(input_sequences, seed=1723,
                         substitution_model='HKY')
            hky_tree = skbio.TreeNode.read(
                         str(hky), convert_underscores=False)
            hky_td = set(hky_tree.tip_tip_distances().to_series())

        # test pairs are not equivalent
        self.assertNotEqual(gtrg_td, hky_td)

    def test_build_iqtree_command(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with tempfile.TemporaryDirectory() as temp_dir:
            run_prefix = os.path.join(temp_dir, 'q2iqtree')
            with redirected_stdio(stderr=os.devnull):
                obs = _build_iqtree_command(input_sequences,
                                            seed=1723,
                                            n_cores='auto',
                                            n_runs=2,
                                            substitution_model='MFP',
                                            run_prefix=run_prefix,
                                            dtype='DNA',
                                            safe='True',
                                            fast='True',
                                            alrt=1000,
                                            abayes=True,
                                            lbp=1000,
                                            n_init_pars_trees=200,
                                            n_top_init_trees=30,
                                            n_best_retain_trees=10,
                                            n_iter=80,
                                            stop_iter=300,
                                            perturb_nni_strength=0.55,
                                            spr_radius=8,
                                            allnni='True')
        self.assertTrue('DNA' in obs[2])
        self.assertTrue('2' in obs[4])
        self.assertTrue(str(input_sequences) in str(obs[6]))
        self.assertTrue('MFP' in obs[8])
        self.assertTrue(str(run_prefix) in obs[10])
        self.assertTrue('AUTO' in obs[12])
        self.assertTrue('1723' in obs[14])
        self.assertTrue(str('-safe') in obs[15])
        self.assertTrue(str('-fast') in obs[16])
        self.assertTrue(str('1000') in obs[18])
        self.assertTrue(str('-abayes') in obs[19])
        self.assertTrue(str('1000') in obs[21])
        self.assertTrue(str('-allnni') in obs[22])
        self.assertTrue(str('200') in obs[24])
        self.assertTrue(str('30') in obs[26])
        self.assertTrue(str('10') in obs[28])
        self.assertTrue(str('80') in obs[30])
        self.assertTrue(str('300') in obs[32])
        self.assertTrue(str('0.55') in obs[34])
        self.assertTrue(str('8') in obs[36])

    def test_build_iqtree_ufbs_command(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with tempfile.TemporaryDirectory() as temp_dir:
            run_prefix = os.path.join(temp_dir, 'q2iqtreeufboot')
            with redirected_stdio(stderr=os.devnull):
                obs = _build_iqtree_ufbs_command(input_sequences,
                                                 seed=1723,
                                                 n_cores='auto',
                                                 n_runs=5,
                                                 bootstrap_replicates=2000,
                                                 substitution_model='MFP',
                                                 run_prefix=run_prefix,
                                                 dtype='DNA',
                                                 safe='True',
                                                 allnni='True',
                                                 alrt=500,
                                                 abayes=True,
                                                 lbp=400,
                                                 bnni=True,
                                                 n_init_pars_trees=200,
                                                 n_top_init_trees=30,
                                                 n_best_retain_trees=10,
                                                 stop_iter=300,
                                                 perturb_nni_strength=0.55,
                                                 spr_radius=8,
                                                 n_max_ufboot_iter=600,
                                                 n_ufboot_steps=80,
                                                 min_cor_ufboot=0.66,
                                                 ep_break_ufboot=0.51)
        self.assertTrue('2000' in obs[2])
        self.assertTrue('DNA' in obs[4])
        self.assertTrue('5' in obs[6])
        self.assertTrue(str(input_sequences) in str(obs[8]))
        self.assertTrue('MFP' in obs[10])
        self.assertTrue(str(run_prefix) in obs[12])
        self.assertTrue('AUTO' in obs[14])
        self.assertTrue('1723' in obs[16])
        self.assertTrue('-safe' in obs[17])
        self.assertTrue('-allnni' in obs[18])
        self.assertTrue('500' in obs[20])
        self.assertTrue('-abayes' in obs[21])
        self.assertTrue('400' in obs[23])
        self.assertTrue('-bnni' in obs[24])
        self.assertTrue('200' in obs[26])
        self.assertTrue(str('30') in obs[28])
        self.assertTrue(str('10') in obs[30])
        self.assertTrue(str('300') in obs[32])
        self.assertTrue(str('0.55') in obs[34])
        self.assertTrue(str('8') in obs[36])
        self.assertTrue(str('600') in obs[38])
        self.assertTrue(str('80') in obs[40])
        self.assertTrue(str('0.66') in obs[42])
        self.assertTrue(str('0.51') in obs[44])

    def test_iqtree_ultrafast_bootstrap(self):
        # Test that output tree is made.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = iqtree_ultrafast_bootstrap(input_sequences)
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

    def test_iqtree_ultrafast_bootstrap_singlebranch_methods(self):
        # Comparing branch support to manually constructed tree
        # using the following command:
        #  iqtree -s aligned-dna-sequences-3.fasta -alrt 1500 -lbp 1500
        #         -abayes -bb 1500 -m 'HKY' -seed 1723
        # Here I am simply checking if the support values are identical
        # to the manual run. Also check for number of values.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            obs = iqtree_ultrafast_bootstrap(input_sequences, seed=1723,
                                             substitution_model='HKY',
                                             alrt=1500, lbp=1500,
                                             abayes=True,
                                             bootstrap_replicates=1500)
        obs_tree = skbio.TreeNode.read(str(obs), convert_underscores=False)
        obs_supp = [node.name for node in obs_tree.non_tips()]

        exp_tree = skbio.TreeNode.read(self.get_data_path('test5.tre'))
        exp_supp = [node.name for node in exp_tree.non_tips()]

        self.assertEqual(set(obs_supp), set(exp_supp))
        self.assertEqual(len(obs_supp[0].split('/')), 4)  # should be 4 values
        self.assertEqual(len(exp_supp[0].split('/')), 4)  # should be 4 values

    def test_iqtree_ultrafast_bootstrap_safe_allnni(self):
        # Test that output tree is made.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = iqtree_ultrafast_bootstrap(input_sequences, safe='True',
                                             allnni='True')
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

    def test_iqtree_ultrafast_bootstrap_auto_threads(self):
        # Test that output tree is made after auto optimizing threads.
        # Reads tree output and compares tip labels to expected labels.
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        with redirected_stdio(stderr=os.devnull):
            obs = iqtree_ultrafast_bootstrap(input_sequences, n_cores='auto')
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


class IqtreeRunCommandTests(TestPluginBase):

    package = 'q2_phylogeny.tests'

    def test_run_not_verbose(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        aligned_fp = str(input_sequences)

        with tempfile.TemporaryDirectory() as temp_dir:
            run_prefix = os.path.join(temp_dir, 'q2iqtree')
            cmd = ['iqtree',
                   '-m', 'HKY',
                   '-seed', '1723',
                   '-s', aligned_fp,
                   '-pre', run_prefix,
                   '-nt', '2']

            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, verbose=False)
            obs_tree_fp = run_prefix + '.treefile'
            obs_tree = skbio.TreeNode.read(str(obs_tree_fp),
                                           convert_underscores=False)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515',
                              'GCA000454205', 'GCA000473545',
                              'GCA000196255', 'GCA002142615',
                              'GCA000686145', 'GCA001950115',
                              'GCA001971985', 'GCA900007555']))

    def test_run_ultrafast_bs_not_verbose(self):
        input_fp = self.get_data_path('aligned-dna-sequences-3.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        aligned_fp = str(input_sequences)

        with tempfile.TemporaryDirectory() as temp_dir:
            run_prefix = os.path.join(temp_dir, 'q2iqtreeufboot')
            cmd = ['iqtree',
                   '-m', 'HKY',
                   '-seed', '1723',
                   '-bb', '1000',
                   '-s', aligned_fp,
                   '-pre', run_prefix,
                   '-nt', '2']

            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, verbose=False)
            obs_tree_fp = run_prefix + '.treefile'
            obs_tree = skbio.TreeNode.read(str(obs_tree_fp),
                                           convert_underscores=False)
        # load the resulting tree and test that it has the right number of
        # tips and the right tip ids
        tips = list(obs_tree.tips())
        tip_names = [t.name for t in tips]
        self.assertEqual(set(tip_names),
                         set(['GCA001510755', 'GCA001045515',
                              'GCA000454205', 'GCA000473545',
                              'GCA000196255', 'GCA002142615',
                              'GCA000686145', 'GCA001950115',
                              'GCA001971985', 'GCA900007555']))


if __name__ == "__main__":
    unittest.main()
