# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin
from q2_types.tree import Phylogeny, Unrooted, Rooted
from q2_types.feature_data import FeatureData, AlignedSequence
from q2_types.feature_table import FeatureTable, Frequency
from qiime2.core.type import Int

import q2_phylogeny

plugin = Plugin(
    name='phylogeny',
    version=q2_phylogeny.__version__,
    website='https://github.com/qiime2/q2-phylogeny',
    package='q2_phylogeny',
    citation_text=("FastTree 2 – Approximately Maximum-Likelihood Trees for "
                   "Large Alignments. Price MN, Dehal PS, Arkin AP (2010) "
                   "PLOS ONE 5(3): e9490. doi: 10.1371/journal.pone.0009490"),
    description=('This QIIME 2 plugin supports generating and manipulating '
                 'phylogenetic trees.'),
    short_description='Plugin for generating and manipulating phylogenies.'
)

plugin.methods.register_function(
    function=q2_phylogeny.midpoint_root,
    inputs={'tree': Phylogeny[Unrooted]},
    parameters={},
    outputs=[('rooted_tree', Phylogeny[Rooted])],
    input_descriptions={'tree': 'The phylogenetic tree to be rooted.'},
    parameter_descriptions={},
    output_descriptions={'rooted_tree': 'The rooted phylogenetic tree.'},
    name='Midpoint root an unrooted phylogenetic tree.',
    description=("Midpoint root an unrooted phylogenetic tree.")
)

plugin.methods.register_function(
    function=q2_phylogeny.fasttree,
    inputs={'alignment': FeatureData[AlignedSequence]},
    parameters={'threads': Int},
    outputs=[('tree', Phylogeny[Unrooted])],
    input_descriptions={
        'alignment': ('Aligned sequences to be used for phylogenetic '
                      'reconstruction.')
    },
    parameter_descriptions={},
    output_descriptions={'tree': 'The resulting phylogenetic tree.'},
    name='Construct a phylogenetic tree with FastTree.',
    description=("Construct a phylogenetic tree with FastTree.")
)

plugin.methods.register_function(
    function=q2_phylogeny.filter_table,
    inputs={'table': FeatureTable[Frequency],
            'tree': Phylogeny[Rooted | Unrooted]},
    parameters={},
    outputs=[('filtered_table', FeatureTable[Frequency])],
    input_descriptions={
        'table': 'Feature table that features should be filtered from.',
        'tree': ('Tree where tip identifiers are the feature identifiers that '
                 'should be retained in the table.')
    },
    parameter_descriptions={},
    output_descriptions={'filtered_table': 'The resulting feature table.'},
    name="Remove features from table if they're not present in tree.",
    description=("Remove features from a feature table if their identifiers "
                 "are not tip identifiers in tree.")
)
