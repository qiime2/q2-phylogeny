# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess

from q2_types.feature_data import AlignedDNAFASTAFormat
from q2_types.tree import NewickFormat


def run_command(cmd, output_fp, verbose=True, env=None):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')

    with open(output_fp, 'w') as output_f:
        subprocess.run(cmd, stdout=output_f, check=True, env=env)


def fasttree(alignment: AlignedDNAFASTAFormat,
             n_threads: int=1) -> NewickFormat:
    result = NewickFormat()
    aligned_fp = str(alignment)
    tree_fp = str(result)

    env = None
    if n_threads == 1:
        cmd = ['FastTree']
    else:
        env = os.environ.copy()
        env.update({'OMP_NUM_THREADS': str(n_threads)})
        cmd = ['FastTreeMP']

    cmd.extend(['-quote', '-nt', aligned_fp])
    run_command(cmd, tree_fp, env=env)
    return result
