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
#from q2_phylogeny._iqtree import _IQTREE_DNA_MODELS

_iqtree_defaults = {
    'n_threads': 1,
    'substitution_model': 'MFP',
    'run_prefix': 'q2-iqtree',
    'dtype': 'DNA',
    'safe': False
}


def _build_iqtree_command(alignment, seed,
                          n_threads=_iqtree_defaults['n_threads'],
                          substitution_model=_iqtree_defaults[
                                             'substitution_model'],
                          run_prefix=_iqtree_defaults['run_prefix'],
                          dtype=_iqtree_defaults['dtype'],
                          safe=_iqtree_defaults['safe']):
    cmd = [
        'iqtree',
        '-m', str(substitution_model),
        '-nt', str(n_threads),
        '-s', str(alignment),
        '-st', str(dtype),
        '-pre', str(run_prefix),
        '-seed', str(seed)
            ]

    if safe:
        cmd += ['-safe']

    #test if substitution_model not in _IQTREE_DNA_MODELS:
    return cmd


def iqtree(alignment: AlignedDNAFASTAFormat,
          seed: int=randint(1000, 10000),
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
