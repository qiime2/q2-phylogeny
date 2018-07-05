# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

from random import randint

from q2_types.feature_data import AlignedDNAFASTAFormat
from q2_types.tree import NewickFormat
from q2_phylogeny._raxml import run_command


def _build_iqtree_command(alignment, seed,
                          n_threads=1,
                          substitution_model='MFP',
                          run_prefix='q2iqtree',
                          dtype='DNA',
                          n_init_pars_trees=None,
                          n_top_init_trees=None,
                          n_best_retain_trees=None,
                          n_iter=None,
                          stop_iter=None,
                          perturb_nni_strength=None,
                          spr_radius=None,
                          allnni=False,
                          safe=False):

    cmd = ['iqtree', '-nt', '%i' % n_threads]

    if seed is None:
        seed = randint(1000, 10000)

    cmd += ['-seed', str(seed),
            '-st', str(dtype),
            '-s', str(alignment),
            '-m', str(substitution_model),
            '-pre', str(run_prefix)]

    if safe:
        cmd += ['-safe']

    if allnni:
        cmd += ['-allnni']

    if n_init_pars_trees:
        cmd += ['-ninit', str(n_init_pars_trees)]

    if n_top_init_trees:
        cmd += ['-ntop', str(n_top_init_trees)]

    if n_best_retain_trees:
        cmd += ['-nbest', str(n_best_retain_trees)]

    if n_iter:
        cmd += ['-n', str(n_iter)]

    if stop_iter:
        cmd += ['-nstop', str(stop_iter)]

    if perturb_nni_strength:
        cmd += ['-pers', str(perturb_nni_strength)]

    if spr_radius:
        cmd += ['-sprrad', str(spr_radius)]

    return cmd


def iqtree(alignment: AlignedDNAFASTAFormat,
           seed: int=None,
           n_threads: int=1,
           substitution_model: str='MFP',
           n_init_pars_trees: int=None,
           n_top_init_trees: int=None,
           n_best_retain_trees: int=None,
           n_iter: int=None,
           stop_iter: int=None,
           perturb_nni_strength: float=None,
           spr_radius: int=None,
           allnni: bool=False,
           safe: bool=False) -> NewickFormat:
    result = NewickFormat()

    with tempfile.TemporaryDirectory() as temp_dir:
        run_prefix = os.path.join(temp_dir, 'q2iqtree')
        cmd = _build_iqtree_command(alignment, seed,
                                    n_threads=n_threads,
                                    substitution_model=substitution_model,
                                    run_prefix=run_prefix,
                                    n_init_pars_trees=n_init_pars_trees,
                                    n_top_init_trees=n_top_init_trees,
                                    n_best_retain_trees=n_best_retain_trees,
                                    n_iter=n_iter,
                                    stop_iter=stop_iter,
                                    perturb_nni_strength=perturb_nni_strength,
                                    spr_radius=spr_radius,
                                    allnni=allnni,
                                    safe=safe)
        run_command(cmd)

        tree_tmp_fp = os.path.join(temp_dir, '%s.treefile' % run_prefix)
        os.rename(tree_tmp_fp, str(result))

    return result


def _build_iqtree_ufbs_command(alignment, seed,
                               n_threads=1,
                               substitution_model='MFP',
                               bootstrap_replicates=1000,
                               run_prefix='q2iqtree',
                               dtype='DNA',
                               n_init_pars_trees=None,
                               n_top_init_trees=None,
                               n_best_retain_trees=None,
                               stop_iter=None,
                               perturb_nni_strength=None,
                               spr_radius=None,
                               n_max_ufboot_iter=None,
                               n_ufboot_steps=None,
                               min_cor_ufboot=None,
                               ep_break_ufboot=None,
                               allnni=False,
                               safe=False):
    # This is a separate command becuase there are several
    # bootstrap specific options.
    cmd = ['iqtree', '-nt', '%i' % n_threads, '-bb',
           '%i' % bootstrap_replicates]

    if seed is None:
        seed = randint(1000, 10000)

    cmd += ['-seed', str(seed),
            '-st', str(dtype),
            '-s', str(alignment),
            '-m', str(substitution_model),
            '-bb', str(bootstrap_replicates),
            '-pre', str(run_prefix)]

    if safe:
        cmd += ['-safe']

    if allnni:
        cmd += ['-allnni']

    if n_init_pars_trees:
        cmd += ['-ninit', str(n_init_pars_trees)]

    if n_top_init_trees:
        cmd += ['-ntop', str(n_top_init_trees)]

    if n_best_retain_trees:
        cmd += ['-nbest', str(n_best_retain_trees)]

    if stop_iter:
        cmd += ['-nstop', str(stop_iter)]

    if perturb_nni_strength:
        cmd += ['-pers', str(perturb_nni_strength)]

    if spr_radius:
        cmd += ['-sprrad', str(spr_radius)]

    if n_max_ufboot_iter:
        cmd += ['-nm', '%i' % n_max_ufboot_iter]

    if n_ufboot_steps:
        cmd += ['-nstep', '%i' % n_ufboot_steps]

    if min_cor_ufboot:
        cmd += ['-bcor', '%f' % min_cor_ufboot]

    if ep_break_ufboot:
        cmd += ['-beps', '%f' % ep_break_ufboot]

    return cmd


def iqtree_ultrafast_bootstrap(alignment: AlignedDNAFASTAFormat,
                               seed: int=None,
                               n_threads: int=1,
                               substitution_model: str='MFP',
                               bootstrap_replicates: int=1000,
                               n_init_pars_trees: int=None,
                               n_top_init_trees: int=None,
                               n_best_retain_trees: int=None,
                               stop_iter: int=None,
                               perturb_nni_strength: float=None,
                               spr_radius: int=None,
                               n_max_ufboot_iter: int=None,
                               n_ufboot_steps: int=None,
                               min_cor_ufboot: float=None,
                               ep_break_ufboot: float=None,
                               allnni: bool=False,
                               safe: bool=False) -> NewickFormat:
    # NOTE: the IQ-TREE command `-n` (called as `n_iter` in the `iqtree`
    # method) is not compatable with ultrafast_bootstrap `-bb`
    result = NewickFormat()

    with tempfile.TemporaryDirectory() as temp_dir:
        run_prefix = os.path.join(temp_dir, 'q2iqtreeufboot')
        cmd = _build_iqtree_ufbs_command(
                    alignment, seed,
                    n_threads=n_threads,
                    substitution_model=substitution_model,
                    bootstrap_replicates=bootstrap_replicates,
                    run_prefix=run_prefix,
                    n_init_pars_trees=n_init_pars_trees,
                    n_top_init_trees=n_top_init_trees,
                    n_best_retain_trees=n_best_retain_trees,
                    stop_iter=stop_iter,
                    perturb_nni_strength=perturb_nni_strength,
                    spr_radius=spr_radius,
                    n_max_ufboot_iter=n_max_ufboot_iter,
                    n_ufboot_steps=n_ufboot_steps,
                    min_cor_ufboot=min_cor_ufboot,
                    ep_break_ufboot=ep_break_ufboot,
                    allnni=allnni,
                    safe=safe)
        run_command(cmd)
        tree_tmp_fp = os.path.join(temp_dir, '%s.treefile' % run_prefix)
        os.rename(tree_tmp_fp, str(result))

    return result
