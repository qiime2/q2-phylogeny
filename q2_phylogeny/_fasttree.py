# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

from q2_types.feature_data import AlignedDNAFASTAFormat
from q2_types.tree import NewickFormat


def fasttree(alignment: AlignedDNAFASTAFormat) -> NewickFormat:
    result = NewickFormat()
    aligned_fp = str(alignment)
    tree_fp = str(result)
    cmd = ['FastTree', '-nt', '-quiet', aligned_fp]
    # construct phylogenetic tree and write the output file
    with open(tree_fp, 'w') as tree_f:
        subprocess.run(cmd, stdout=tree_f)
    return result
