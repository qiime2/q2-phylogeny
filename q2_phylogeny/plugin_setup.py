# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin, Citations, Int, Range, Str, Choices
from q2_types.tree import Phylogeny, Unrooted, Rooted
from q2_types.feature_data import FeatureData, AlignedSequence
from q2_types.feature_table import FeatureTable, Frequency

import q2_phylogeny

_RAXML_MODEL_OPT = ['GTRGAMMA', 'GTRGAMMAI', 'GTRCAT', 'GTRCATI']

citations = Citations.load('citations.bib', package='q2_phylogeny')
plugin = Plugin(
    name='phylogeny',
    version=q2_phylogeny.__version__,
    website='https://github.com/qiime2/q2-phylogeny',
    package='q2_phylogeny',
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
    parameters={'n_threads': Int % Range(-1, None)},
    outputs=[('tree', Phylogeny[Unrooted])],
    input_descriptions={
        'alignment': ('Aligned sequences to be used for phylogenetic '
                      'reconstruction.')
    },
    parameter_descriptions={
        'n_threads': 'The number of threads. Using more than one thread '
                     'runs the non-deterministic variant of `FastTree` '
                     '(`FastTreeMP`), and may result in a different tree than '
                     'single-threading. See '
                     'http://www.microbesonline.org/fasttree/#OpenMP for '
                     'details. (Use -1 to automatically use all available '
                     'cores)'
    },
    output_descriptions={'tree': 'The resulting phylogenetic tree.'},
    name='Construct a phylogenetic tree with FastTree.',
    description=("Construct a phylogenetic tree with FastTree."),
    citations=[citations['price2010fasttree']]
)

plugin.methods.register_function(
    function=q2_phylogeny.raxml,
    inputs={
            'alignment': FeatureData[AlignedSequence]},
    parameters={
            'seed': Int,
            'n_threads': Int % Range(1, None),
            'substitution_model': Str %
		                      Choices(_RAXML_MODEL_OPT),
		},
    outputs=[('tree', Phylogeny[Unrooted])],
    input_descriptions={
        'alignment': ('Aligned sequences to be used for phylogenetic '
                      'reconstruction.'),
    },
    parameter_descriptions={
        'n_threads': ('The number of threads. Using more than one thread '
                     'will enable the PTHREADS version of RAxML'),
        'substitution_model': ('Model of Nucleotide Substitution'),
        'seed': ('Random number seed for the parsimony starting tree.'
	         'This allows you to reproduce your results.'),
    },
    output_descriptions={'tree': 'The resulting phylogenetic tree.'},
    name='Construct a phylogenetic tree with RAxML.',
    description=("Construct a phylogenetic tree with RAxML."),
    citations=[citations['Stamatakis2014raxml']]
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
