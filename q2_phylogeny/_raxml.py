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


def run_command(cmd, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def raxml(alignment: AlignedDNAFASTAFormat,
          seed: int=None,
          n_searches: int=1,
          n_threads: int=1,
          substitution_model: str='GTRGAMMA') -> NewickFormat:
    result = NewickFormat()

    if n_threads == 1:
        cmd = ['raxmlHPC']
    else:
        cmd = ['raxmlHPC-PTHREADS', '-T %i' % n_threads]

    if seed is None:
        seed = randint(1000, 10000)

    runname = 'q2'
    with tempfile.TemporaryDirectory() as temp_dir:
        cmd += ['-m', str(substitution_model),
                '-p', str(seed),
                '-N', str(n_searches),
                '-s', str(alignment),
                '-w', temp_dir,
                '-n', runname]
        run_command(cmd)

        tree_tmp_fp = os.path.join(temp_dir, 'RAxML_bestTree.%s' % runname)
        os.rename(tree_tmp_fp, str(result))

    return result


def _build_rapid_bootstrap_command(alignment, seed, rapid_bootstrap_seed,
                                   bootstrap_replicates, substitution_model,
                                   temp_dir, runname):
    cmd = ['-f', 'a',  # always set, rapid bootstrapping
           '-m', str(substitution_model),
           '-p', str(seed),
           '-x', str(rapid_bootstrap_seed),
           '-N', str(bootstrap_replicates),
           '-s', str(alignment),
           '-w', temp_dir,
           '-n', runname]
    return cmd


def raxml_rapid_bootstrap(alignment: AlignedDNAFASTAFormat,
                          seed: int=None, rapid_bootstrap_seed: int=None,
                          bootstrap_replicates: int=100, n_threads: int=1,
                          substitution_model: str='GTRGAMMA') -> NewickFormat:
    result = NewickFormat()

    if n_threads == 1:
        cmd = ['raxmlHPC']
    else:
        cmd = ['raxmlHPC-PTHREADS', '-T %i' % n_threads]

    if seed is None:
        seed = randint(1000, 10000)

    if rapid_bootstrap_seed is None:
        rapid_bootstrap_seed = randint(1000, 10000)

    runname = 'q2bootstrap'
    with tempfile.TemporaryDirectory() as temp_dir:
        cmd += _build_rapid_bootstrap_command(alignment, seed,
                                              rapid_bootstrap_seed,
                                              bootstrap_replicates,
                                              substitution_model, temp_dir,
                                              runname)
        run_command(cmd)

        tree_tmp_fp = os.path.join(temp_dir, 'RAxML_bipartitions.%s' % runname)
        os.rename(tree_tmp_fp, str(result))

    return result
