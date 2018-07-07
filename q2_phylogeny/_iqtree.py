# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

from q2_types.feature_data import AlignedDNAFASTAFormat
from q2_types.tree import NewickFormat
from q2_phylogeny._raxml import run_command

_iqtree_defaults = {
    'seed': None,
    'n_cores': 1,
    'substitution_model': 'MFP',
    'run_prefix': 'q2iqtree',
    'dtype': 'DNA',
    'n_init_pars_trees': None,
    'n_top_init_trees': None,
    'n_best_retain_trees': None,
    'n_iter': None,
    'stop_iter': None,
    'perturb_nni_strength': None,
    'spr_radius': None,
    'allnni': False,
    'safe': False,
    'bootstrap_replicates': 1000,
    'n_max_ufboot_iter': None,
    'n_ufboot_steps': None,
    'min_cor_ufboot': None,
    'ep_break_ufboot': None,
}


def _build_iqtree_command(
        alignment,
        seed: int=_iqtree_defaults['seed'],
        n_cores: int=_iqtree_defaults['n_cores'],
        substitution_model: str=_iqtree_defaults['substitution_model'],
        run_prefix: str=_iqtree_defaults['run_prefix'],
        dtype: str=_iqtree_defaults['dtype'],
        n_init_pars_trees: int=_iqtree_defaults['n_init_pars_trees'],
        n_top_init_trees: int=_iqtree_defaults['n_top_init_trees'],
        n_best_retain_trees: int=_iqtree_defaults['n_best_retain_trees'],
        n_iter: int=_iqtree_defaults['n_iter'],
        stop_iter: int=_iqtree_defaults['stop_iter'],
        perturb_nni_strength: float=_iqtree_defaults['perturb_nni_strength'],
        spr_radius: int=_iqtree_defaults['spr_radius'],
        allnni: bool=_iqtree_defaults['allnni'],
        safe: bool=_iqtree_defaults['safe']):

    cmd = ['iqtree']

    cmd += ['-st', str(dtype),
            '-s', str(alignment),
            '-m', str(substitution_model),
            '-pre', str(run_prefix)]

    if n_cores == 0:
        cmd += ['-nt', 'AUTO']
    else:
        cmd += ['-nt', str(n_cores)]

    if seed:
        cmd += ['-seed', str(seed)]

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


def iqtree(
    alignment: AlignedDNAFASTAFormat,
    seed: int=_iqtree_defaults['seed'],
    n_cores: int=_iqtree_defaults['n_cores'],
    substitution_model: str=_iqtree_defaults['substitution_model'],
    n_init_pars_trees: int=_iqtree_defaults['n_init_pars_trees'],
    n_top_init_trees: int=_iqtree_defaults['n_top_init_trees'],
    n_best_retain_trees: int=_iqtree_defaults['n_best_retain_trees'],
    n_iter: int=_iqtree_defaults['n_iter'],
    stop_iter: int=_iqtree_defaults['stop_iter'],
    perturb_nni_strength: float=_iqtree_defaults['perturb_nni_strength'],
    spr_radius: int=_iqtree_defaults['spr_radius'],
    allnni: bool=_iqtree_defaults['allnni'],
    safe: bool=_iqtree_defaults['safe'],
            ) -> NewickFormat:
    result = NewickFormat()

    with tempfile.TemporaryDirectory() as temp_dir:
        run_prefix = os.path.join(temp_dir, 'q2iqtree')
        cmd = _build_iqtree_command(alignment,
                                    seed=seed,
                                    n_cores=n_cores,
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


def _build_iqtree_ufbs_command(
        alignment,
        seed: int=_iqtree_defaults['seed'],
        n_cores: int=_iqtree_defaults['n_cores'],
        substitution_model: str=_iqtree_defaults['substitution_model'],
        bootstrap_replicates: int=_iqtree_defaults['bootstrap_replicates'],
        run_prefix: str=_iqtree_defaults['run_prefix'],
        dtype: str=_iqtree_defaults['dtype'],
        n_init_pars_trees: int=_iqtree_defaults['n_init_pars_trees'],
        n_top_init_trees: int=_iqtree_defaults['n_top_init_trees'],
        n_best_retain_trees: int=_iqtree_defaults['n_best_retain_trees'],
        stop_iter: int=_iqtree_defaults['stop_iter'],
        perturb_nni_strength: float=_iqtree_defaults['perturb_nni_strength'],
        spr_radius: int=_iqtree_defaults['spr_radius'],
        n_max_ufboot_iter: int=_iqtree_defaults['n_max_ufboot_iter'],
        n_ufboot_steps: int=_iqtree_defaults['n_ufboot_steps'],
        min_cor_ufboot: float=_iqtree_defaults['min_cor_ufboot'],
        ep_break_ufboot: float=_iqtree_defaults['ep_break_ufboot'],
        allnni: bool=_iqtree_defaults['allnni'],
        safe: bool=_iqtree_defaults['safe']):
    # This is a separate command becuase there are several
    # bootstrap specific options.

    cmd = ['iqtree', '-bb', '%i' % bootstrap_replicates]

    cmd += ['-st', str(dtype),
            '-s', str(alignment),
            '-m', str(substitution_model),
            '-pre', str(run_prefix)]

    if n_cores == 0:
        cmd += ['-nt', 'AUTO']
    else:
        cmd += ['-nt', str(n_cores)]

    if seed:
        cmd += ['-seed', str(seed)]

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


def iqtree_ultrafast_bootstrap(
    alignment: AlignedDNAFASTAFormat,
    seed: int=_iqtree_defaults['seed'],
    n_cores: int=_iqtree_defaults['n_cores'],
    substitution_model: str=_iqtree_defaults['substitution_model'],
    bootstrap_replicates: int=_iqtree_defaults['bootstrap_replicates'],
    n_init_pars_trees: int=_iqtree_defaults['n_init_pars_trees'],
    n_top_init_trees: int=_iqtree_defaults['n_top_init_trees'],
    n_best_retain_trees: int=_iqtree_defaults['n_best_retain_trees'],
    stop_iter: int=_iqtree_defaults['stop_iter'],
    perturb_nni_strength: float=_iqtree_defaults['perturb_nni_strength'],
    spr_radius: int=_iqtree_defaults['spr_radius'],
    n_max_ufboot_iter: int=_iqtree_defaults['n_max_ufboot_iter'],
    n_ufboot_steps: int=_iqtree_defaults['n_ufboot_steps'],
    min_cor_ufboot: float=_iqtree_defaults['min_cor_ufboot'],
    ep_break_ufboot: float=_iqtree_defaults['ep_break_ufboot'],
    allnni: bool=_iqtree_defaults['allnni'],
    safe: bool=_iqtree_defaults['safe']
                                ) -> NewickFormat:
    # NOTE: the IQ-TREE command `-n` (called as `n_iter` in the `iqtree`
    # method) is not compatable with ultrafast_bootstrap `-bb`
    result = NewickFormat()

    with tempfile.TemporaryDirectory() as temp_dir:
        run_prefix = os.path.join(temp_dir, 'q2iqtreeufboot')
        cmd = _build_iqtree_ufbs_command(
                    alignment,
                    seed=seed,
                    n_cores=n_cores,
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
