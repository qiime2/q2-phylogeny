# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from ._util import midpoint_root
from ._fasttree import fasttree
from ._filter import filter_table


__version__ = pkg_resources.get_distribution('q2-phylogeny').version

__all__ = ["midpoint_root", "fasttree", "filter_table"]
