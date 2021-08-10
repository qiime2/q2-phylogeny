# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._util import midpoint_root, robinson_foulds
from ._fasttree import fasttree
from ._raxml import raxml, raxml_rapid_bootstrap
from ._iqtree import iqtree, iqtree_ultrafast_bootstrap
from ._filter import filter_table, filter_tree
from ._version import get_versions
from ._align_to_tree_mafft_fasttree import align_to_tree_mafft_fasttree
from ._align_to_tree_mafft_iqtree import align_to_tree_mafft_iqtree
from ._align_to_tree_mafft_raxml import align_to_tree_mafft_raxml

__version__ = get_versions()['version']
del get_versions

__all__ = ["midpoint_root", "fasttree", "align_to_tree_mafft_fasttree",
           "raxml", "raxml_rapid_bootstrap", "iqtree", "filter_table",
           "iqtree_ultrafast_bootstrap", "align_to_tree_mafft_iqtree",
           "align_to_tree_mafft_raxml", "robinson_foulds", 'filter_tree']
