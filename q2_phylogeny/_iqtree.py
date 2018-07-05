# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import subprocess

from random import randint

from q2_types.feature_data import AlignedDNAFASTAFormat
from q2_types.tree import NewickFormat

from q2_phylogeny._raxml import run_command

def _build_iqtree_command(alignment, seed,
                          n_threads=1,
                          substitution_model='MFP',
                          run_prefix='q2iqtree',
                          dtype='DNA',
                          safe=False):

    cmd = ['iqtree','-nt', '%i' % n_threads]

    if seed is None:
        seed = randint(1000, 10000)

    cmd += ['-seed', str(seed),
        '-st', str(dtype),
        '-s', str(alignment),
        '-m', str(substitution_model),
        '-pre', str(run_prefix)]

    if safe:
        cmd += ['-safe']

    return cmd


def iqtree(alignment: AlignedDNAFASTAFormat,
          seed: int=None,
          n_threads: int=1,
          substitution_model: str='MFP',
          safe: bool=False) -> NewickFormat:
    result = NewickFormat()

    with tempfile.TemporaryDirectory() as temp_dir:
        run_prefix = os.path.join(temp_dir, 'q2iqtree')
        cmd = _build_iqtree_command(alignment, seed,
                                     n_threads=n_threads,
                                     substitution_model=substitution_model,
                                     run_prefix=run_prefix, safe=safe)

        run_command(cmd)

        tree_tmp_fp = os.path.join(temp_dir, '%s.treefile' % run_prefix)
        os.rename(tree_tmp_fp, str(result))

    return result


def _build_iqtree_ultrafast_bootstrap_command(alignment, seed,
                          n_threads=1,
                          substitution_model='MFP',
                          bootstrap_replicates=1000,
                          run_prefix='q2iqtree',
                          dtype='DNA',
                          safe=False):
    # This is a separate command becuase there are many other
    # bootstrap specific option we may want to add later.
    cmd = ['iqtree','-nt', '%i' % n_threads, '-bb', '%i' % bootstrap_replicates]

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

    return cmd


def iqtree_ultrafast_bootstrap(alignment: AlignedDNAFASTAFormat,
          seed: int=None,
          n_threads: int=1,
          substitution_model: str='MFP',
          bootstrap_replicates: int=1000,
          safe: bool=False) -> NewickFormat:
    result = NewickFormat()

    with tempfile.TemporaryDirectory() as temp_dir:
        run_prefix = os.path.join(temp_dir, 'q2iqtreeufboot')
        cmd = _build_iqtree_ultrafast_bootstrap_command(alignment, seed,
                                    n_threads=n_threads,
                                    substitution_model=substitution_model,
                                    bootstrap_replicates=bootstrap_replicates,
                                    run_prefix=run_prefix, safe=safe)
        run_command(cmd)

        tree_tmp_fp = os.path.join(temp_dir, '%s.treefile' % run_prefix)
        os.rename(tree_tmp_fp, str(result))

    return result
