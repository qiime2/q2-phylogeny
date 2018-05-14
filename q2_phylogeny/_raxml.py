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

import skbio

from q2_types.feature_data import AlignedDNAFASTAFormat
from q2_types.tree import NewickFormat


def run_command(cmd, verbose=True, env=None):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')

    subprocess.run(cmd, check=True, env=env)


def raxml(alignment: AlignedDNAFASTAFormat,
             n_threads: int=1, seed: int=1723,
	     substitution_model: str='GTRGAMMA') -> NewickFormat:
    aligned_fp = str(alignment)
    result = NewickFormat()

    if n_threads == 1:
        cmd = ['raxmlHPC']
    else:
        cmd = ['raxmlHPC-PTHREADS', '-T %s' % n_threads]
        #cmd = ['raxmlHPC-PTHREADS']

    runname = 'q2'
    with tempfile.TemporaryDirectory() as temp_dir:
        cmd += ['-m', str(substitution_model),
	        '-p', str(seed),
	        '-s', aligned_fp,
	        '-w', temp_dir,
	        '-n', runname
	       ]
        run_command(cmd)

        tree_tmp_fp = os.path.join(temp_dir, 'RAxML_bestTree.%s' % runname)
        os.rename(tree_tmp_fp, str(result))

    return result
