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

_IQTREE_DNA_MODELS = ['JC', 'JC+I', 'JC+G', 'JC+I+G', 'JC+R2', 'JC+R3', 'JC+R4',
'JC+R5', 'JC+R6', 'JC+R7', 'JC+R8', 'JC+R9', 'JC+R10', 'F81', 'F81+I',
'F81+G', 'F81+I+G', 'F81+R2', 'F81+R3', 'F81+R4', 'F81+R5', 'F81+R6',
'F81+R7', 'F81+R8', 'F81+R9', 'F81+R10', 'K80', 'K80+I', 'K80+G', 'K80+I+G',
'K80+R2', 'K80+R3', 'K80+R4', 'K80+R5', 'K80+R6', 'K80+R7', 'K80+R8',
'K80+R9', 'K80+R10', 'HKY', 'HKY+I', 'HKY+G', 'HKY+I+G', 'HKY+R2', 'HKY+R3',
'HKY+R4', 'HKY+R5', 'HKY+R6', 'HKY+R7', 'HKY+R8', 'HKY+R9', 'HKY+R10', 'TNe',
'TNe+I', 'TNe+G', 'TNe+I+G', 'TNe+R2', 'TNe+R3', 'TNe+R4', 'TNe+R5', 'TNe+R6',
'TNe+R7', 'TNe+R8', 'TNe+R9', 'TNe+R10', 'TN', 'TN+I', 'TN+G', 'TN+I+G',
'TN+R2', 'TN+R3', 'TN+R4', 'TN+R5', 'TN+R6', 'TN+R7', 'TN+R8', 'TN+R9',
'TN+R10', 'K81', 'K81+I', 'K81+G', 'K81+I+G', 'K81+R2', 'K81+R3', 'K81+R4',
'K81+R5', 'K81+R6', 'K81+R7', 'K81+R8', 'K81+R9', 'K81+R10', 'K81u', 'K81u+I',
'K81u+G', 'K81u+I+G', 'K81u+R2', 'K81u+R3', 'K81u+R4', 'K81u+R5', 'K81u+R6',
'K81u+R7', 'K81u+R8', 'K81u+R9', 'K81u+R10', 'TPM2', 'TPM2+I', 'TPM2+G',
'TPM2+I+G', 'TPM2+R2', 'TPM2+R3', 'TPM2+R4', 'TPM2+R5', 'TPM2+R6', 'TPM2+R7',
'TPM2+R8', 'TPM2+R9', 'TPM2+R10', 'TPM2u', 'TPM2u+I', 'TPM2u+G', 'TPM2u+I+G',
'TPM2u+R2', 'TPM2u+R3', 'TPM2u+R4', 'TPM2u+R5', 'TPM2u+R6', 'TPM2u+R7',
'TPM2u+R8', 'TPM2u+R9', 'TPM2u+R10', 'TPM3', 'TPM3+I', 'TPM3+G', 'TPM3+I+G',
'TPM3+R2', 'TPM3+R3', 'TPM3+R4', 'TPM3+R5', 'TPM3+R6', 'TPM3+R7', 'TPM3+R8',
'TPM3+R9', 'TPM3+R10', 'TPM3u', 'TPM3u+I', 'TPM3u+G', 'TPM3u+I+G', 'TPM3u+R2',
'TPM3u+R3', 'TPM3u+R4', 'TPM3u+R5', 'TPM3u+R6', 'TPM3u+R7', 'TPM3u+R8',
'TPM3u+R9', 'TPM3u+R10', 'TIMe', 'TIMe+I', 'TIMe+G', 'TIMe+I+G', 'TIMe+R2',
'TIMe+R3', 'TIMe+R4', 'TIMe+R5', 'TIMe+R6', 'TIMe+R7', 'TIMe+R8', 'TIMe+R9',
'TIMe+R10', 'TIM', 'TIM+I', 'TIM+G', 'TIM+I+G', 'TIM+R2', 'TIM+R3', 'TIM+R4',
'TIM+R5', 'TIM+R6', 'TIM+R7', 'TIM+R8', 'TIM+R9', 'TIM+R10', 'TIM2e',
'TIM2e+I', 'TIM2e+G', 'TIM2e+I+G', 'TIM2e+R2', 'TIM2e+R3', 'TIM2e+R4',
'TIM2e+R5', 'TIM2e+R6', 'TIM2e+R7', 'TIM2e+R8', 'TIM2e+R9', 'TIM2e+R10',
'TIM2', 'TIM2+I', 'TIM2+G', 'TIM2+I+G', 'TIM2+R2', 'TIM2+R3', 'TIM2+R4',
'TIM2+R5', 'TIM2+R6', 'TIM2+R7', 'TIM2+R8', 'TIM2+R9', 'TIM2+R10', 'TIM3e',
'TIM3e+I', 'TIM3e+G', 'TIM3e+I+G', 'TIM3e+R2', 'TIM3e+R3', 'TIM3e+R4',
'TIM3e+R5', 'TIM3e+R6', 'TIM3e+R7', 'TIM3e+R8', 'TIM3e+R9', 'TIM3e+R10',
'TIM3', 'TIM3+I', 'TIM3+G', 'TIM3+I+G', 'TIM3+R2', 'TIM3+R3', 'TIM3+R4',
'TIM3+R5', 'TIM3+R6', 'TIM3+R7', 'TIM3+R8', 'TIM3+R9', 'TIM3+R10', 'TVMe',
'TVMe+I', 'TVMe+G', 'TVMe+I+G', 'TVMe+R2', 'TVMe+R3', 'TVMe+R4', 'TVMe+R5',
'TVMe+R6', 'TVMe+R7', 'TVMe+R8', 'TVMe+R9', 'TVMe+R10', 'TVM', 'TVM+I',
'TVM+G', 'TVM+I+G', 'TVM+R2', 'TVM+R3', 'TVM+R4', 'TVM+R5', 'TVM+R6',
'TVM+R7', 'TVM+R8', 'TVM+R9', 'TVM+R10', 'SYM', 'SYM+I', 'SYM+G', 'SYM+I+G',
'SYM+R2', 'SYM+R3', 'SYM+R4', 'SYM+R5', 'SYM+R6', 'SYM+R7', 'SYM+R8',
'SYM+R9', 'SYM+R10', 'GTR', 'GTR+I', 'GTR+G', 'GTR+I+G', 'GTR+R2', 'GTR+R3',
'GTR+R4', 'GTR+R5', 'GTR+R6', 'GTR+R7', 'GTR+R8', 'GTR+R9', 'GTR+R10', 'MFP',
'TEST']

def _build_iqtree_command(alignment, seed,
                          n_threads=1,
                          substitution_model='MFP',
                          run_prefix='q2iqtree',
                          dtype='DNA',
                          safe=False):
    cmd = [
        'iqtree',
        '-m', str(substitution_model),
        '-nt', str(n_threads),
        '-s', str(alignment),
        '-st', str(dtype),
        '-pre', str(run_prefix),
            ]

    if safe:
        cmd += ['-safe']

    if seed is None:
        cmd += ['-seed', str(randint(1000, 10000))]

    # if substitution_model not in _IQTREE_DNA_MODELS:
    #     print("\'%s\' is not one of the allowed models."
    #             %(substitution_model))
    #     print("Allowed substitution models are:\n%s"
    #             %('\n'.join(_IQTREE_DNA_MODELS)))

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
