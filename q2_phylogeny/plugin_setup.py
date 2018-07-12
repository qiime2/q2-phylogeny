# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Citations, Int, Range, Str, Choices, Bool,
                           Float)
from q2_types.tree import Phylogeny, Unrooted, Rooted
from q2_types.feature_data import FeatureData, AlignedSequence, Sequence
from q2_types.feature_table import FeatureTable, Frequency

import q2_phylogeny

_RAXML_MODEL_OPT = ['GTRGAMMA', 'GTRGAMMAI', 'GTRCAT', 'GTRCATI']
_IQTREE_DNA_MODELS = ['JC', 'JC+I', 'JC+G', 'JC+I+G', 'JC+R2', 'JC+R3',
                      'JC+R4', 'JC+R5', 'JC+R6', 'JC+R7', 'JC+R8', 'JC+R9',
                      'JC+R10', 'F81', 'F81+I', 'F81+G', 'F81+I+G', 'F81+R2',
                      'F81+R3', 'F81+R4', 'F81+R5', 'F81+R6', 'F81+R7',
                      'F81+R8', 'F81+R9', 'F81+R10', 'K80', 'K80+I', 'K80+G',
                      'K80+I+G', 'K80+R2', 'K80+R3', 'K80+R4', 'K80+R5',
                      'K80+R6', 'K80+R7', 'K80+R8', 'K80+R9', 'K80+R10',
                      'HKY', 'HKY+I', 'HKY+G', 'HKY+I+G', 'HKY+R2', 'HKY+R3',
                      'HKY+R4', 'HKY+R5', 'HKY+R6', 'HKY+R7', 'HKY+R8',
                      'HKY+R9', 'HKY+R10', 'TNe', 'TNe+I', 'TNe+G',
                      'TNe+I+G', 'TNe+R2', 'TNe+R3', 'TNe+R4', 'TNe+R5',
                      'TNe+R6', 'TNe+R7', 'TNe+R8', 'TNe+R9', 'TNe+R10', 'TN',
                      'TN+I', 'TN+G', 'TN+I+G', 'TN+R2', 'TN+R3', 'TN+R4',
                      'TN+R5', 'TN+R6', 'TN+R7', 'TN+R8', 'TN+R9', 'TN+R10',
                      'K81', 'K81+I', 'K81+G', 'K81+I+G', 'K81+R2', 'K81+R3',
                      'K81+R4', 'K81+R5', 'K81+R6', 'K81+R7', 'K81+R8',
                      'K81+R9', 'K81+R10', 'K81u', 'K81u+I', 'K81u+G',
                      'K81u+I+G', 'K81u+R2', 'K81u+R3', 'K81u+R4', 'K81u+R5',
                      'K81u+R6', 'K81u+R7', 'K81u+R8', 'K81u+R9', 'K81u+R10',
                      'TPM2', 'TPM2+I', 'TPM2+G', 'TPM2+I+G', 'TPM2+R2',
                      'TPM2+R3', 'TPM2+R4', 'TPM2+R5', 'TPM2+R6', 'TPM2+R7',
                      'TPM2+R8', 'TPM2+R9', 'TPM2+R10', 'TPM2u', 'TPM2u+I',
                      'TPM2u+G', 'TPM2u+I+G', 'TPM2u+R2', 'TPM2u+R3',
                      'TPM2u+R4', 'TPM2u+R5', 'TPM2u+R6', 'TPM2u+R7',
                      'TPM2u+R8', 'TPM2u+R9', 'TPM2u+R10', 'TPM3', 'TPM3+I',
                      'TPM3+G', 'TPM3+I+G', 'TPM3+R2', 'TPM3+R3', 'TPM3+R4',
                      'TPM3+R5', 'TPM3+R6', 'TPM3+R7', 'TPM3+R8', 'TPM3+R9',
                      'TPM3+R10', 'TPM3u', 'TPM3u+I', 'TPM3u+G', 'TPM3u+I+G',
                      'TPM3u+R2', 'TPM3u+R3', 'TPM3u+R4', 'TPM3u+R5',
                      'TPM3u+R6', 'TPM3u+R7', 'TPM3u+R8', 'TPM3u+R9',
                      'TPM3u+R10', 'TIMe', 'TIMe+I', 'TIMe+G', 'TIMe+I+G',
                      'TIMe+R2', 'TIMe+R3', 'TIMe+R4', 'TIMe+R5', 'TIMe+R6',
                      'TIMe+R7', 'TIMe+R8', 'TIMe+R9', 'TIMe+R10', 'TIM',
                      'TIM+I', 'TIM+G', 'TIM+I+G', 'TIM+R2', 'TIM+R3',
                      'TIM+R4', 'TIM+R5', 'TIM+R6', 'TIM+R7', 'TIM+R8',
                      'TIM+R9', 'TIM+R10', 'TIM2e', 'TIM2e+I', 'TIM2e+G',
                      'TIM2e+I+G', 'TIM2e+R2', 'TIM2e+R3', 'TIM2e+R4',
                      'TIM2e+R5', 'TIM2e+R6', 'TIM2e+R7', 'TIM2e+R8',
                      'TIM2e+R9', 'TIM2e+R10', 'TIM2', 'TIM2+I', 'TIM2+G',
                      'TIM2+I+G', 'TIM2+R2', 'TIM2+R3', 'TIM2+R4', 'TIM2+R5',
                      'TIM2+R6', 'TIM2+R7', 'TIM2+R8', 'TIM2+R9', 'TIM2+R10',
                      'TIM3e', 'TIM3e+I', 'TIM3e+G', 'TIM3e+I+G', 'TIM3e+R2',
                      'TIM3e+R3', 'TIM3e+R4', 'TIM3e+R5', 'TIM3e+R6',
                      'TIM3e+R7', 'TIM3e+R8', 'TIM3e+R9', 'TIM3e+R10', 'TIM3',
                      'TIM3+I', 'TIM3+G', 'TIM3+I+G', 'TIM3+R2', 'TIM3+R3',
                      'TIM3+R4', 'TIM3+R5', 'TIM3+R6', 'TIM3+R7', 'TIM3+R8',
                      'TIM3+R9', 'TIM3+R10', 'TVMe', 'TVMe+I', 'TVMe+G',
                      'TVMe+I+G', 'TVMe+R2', 'TVMe+R3', 'TVMe+R4', 'TVMe+R5',
                      'TVMe+R6', 'TVMe+R7', 'TVMe+R8', 'TVMe+R9', 'TVMe+R10',
                      'TVM', 'TVM+I', 'TVM+G', 'TVM+I+G', 'TVM+R2', 'TVM+R3',
                      'TVM+R4', 'TVM+R5', 'TVM+R6', 'TVM+R7', 'TVM+R8',
                      'TVM+R9', 'TVM+R10', 'SYM', 'SYM+I', 'SYM+G', 'SYM+I+G',
                      'SYM+R2', 'SYM+R3', 'SYM+R4', 'SYM+R5', 'SYM+R6',
                      'SYM+R7', 'SYM+R8', 'SYM+R9', 'SYM+R10', 'GTR', 'GTR+I',
                      'GTR+G', 'GTR+I+G', 'GTR+R2', 'GTR+R3', 'GTR+R4',
                      'GTR+R5', 'GTR+R6', 'GTR+R7', 'GTR+R8', 'GTR+R9',
                      'GTR+R10', 'MFP', 'TEST']

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
            'n_searches': Int,
            'n_threads': Int % Range(1, None),
            'substitution_model': Str % Choices(_RAXML_MODEL_OPT)},
    outputs=[('tree', Phylogeny[Unrooted])],
    input_descriptions={
        'alignment': ('Aligned sequences to be used for phylogenetic '
                      'reconstruction.'),
    },
    parameter_descriptions={
        'n_searches': ('The number of independent maximum likelihood '
                       'searches to perform. The single best scoring '
                       'tree is returned.'),
        'n_threads': ('The number of threads to use for multithreaded '
                      'processing. Using more than one thread '
                      'will enable the PTHREADS version of RAxML.'),
        'substitution_model': ('Model of Nucleotide Substitution.'),
        'seed': ('Random number seed for the parsimony starting tree. '
                 'This allows you to reproduce tree results. '
                 'If not supplied then one will be randomly chosen.')},
    output_descriptions={'tree': 'The resulting phylogenetic tree.'},
    name='Construct a phylogenetic tree with RAxML.',
    description=("Construct a phylogenetic tree with RAxML."),
    citations=[citations['Stamatakis2014raxml']]
)

plugin.methods.register_function(
    function=q2_phylogeny.raxml_rapid_bootstrap,
    inputs={
            'alignment': FeatureData[AlignedSequence]},
    parameters={
            'seed': Int,
            'rapid_bootstrap_seed': Int,
            'bootstrap_replicates': Int % Range(10, None),
            'n_threads': Int % Range(1, None),
            'substitution_model': Str % Choices(_RAXML_MODEL_OPT)},
    outputs=[('tree', Phylogeny[Unrooted])],
    input_descriptions={
        'alignment': ('Aligned sequences to be used for phylogenetic '
                      'reconstruction.'),
    },
    parameter_descriptions={
        'n_threads': ('The number of threads to use for multithreaded '
                      'processing. Using more than one thread '
                      'will enable the PTHREADS version of RAxML.'),
        'substitution_model': ('Model of Nucleotide Substitution'),
        'seed': ('Random number seed for the parsimony starting tree. '
                 'This allows you to reproduce tree results. '
                 'If not supplied then one will be randomly chosen.'),
        'rapid_bootstrap_seed': ('Specify a random seed for rapid '
                                 'bootstrapping. This allows you to reproduce '
                                 'rapid bootstrap results. If not supplied '
                                 'then one will be randomly chosen.'),
        'bootstrap_replicates': ('The number of bootstrap searches to '
                                 'perform.')},
    output_descriptions={'tree': 'The resulting phylogenetic tree.'},
    name='Construct a phylogenetic tree with bootstrap supports using RAxML.',
    description=('Construct a phylogenetic tree with RAxML with the addition '
                 'of rapid bootstrapping support values.'),
    citations=[citations['Stamatakis2014raxml'],
               citations['Stamatakis2008raxml']]
)

plugin.methods.register_function(
    function=q2_phylogeny.iqtree,
    inputs={'alignment': FeatureData[AlignedSequence]},
    parameters={
            'seed': Int,
            'n_cores': Int,
            'substitution_model': Str % Choices(_IQTREE_DNA_MODELS),
            'n_init_pars_trees': Int % Range(1, None),
            'n_top_init_trees': Int % Range(1, None),
            'n_best_retain_trees': Int % Range(1, None),
            'n_iter': Int % Range(1, None),
            'stop_iter': Int % Range(1, None),
            'perturb_nni_strength': Float % Range(0.01, 99),
            'spr_radius': Int % Range(1, None),
            'allnni': Bool,
            'safe': Bool},
    outputs=[('tree', Phylogeny[Unrooted])],
    input_descriptions={
        'alignment': ('Aligned sequences to be used for phylogenetic '
                      'reconstruction.'),
    },
    parameter_descriptions={
        'n_cores': ('The number of cores to use for parallel '
                    'processing. Use \'0\' to let IQ-TREE automatically '
                    'determine the optimal number of cores to use.'),
        'substitution_model': ('Model of Nucleotide Substitution. '
                               'If not provided, IQ-TREE will determine the '
                               'best fit substitution model automatically.'),
        'seed': ('Random number seed. If not set, program defaults will be '
                 'used. See IQ-TREE manual for details.'),
        'n_init_pars_trees': ('Number of initial parsimony trees. If not '
                              'set, program defaults will be used. See '
                              'IQ-TREE manual for details.'),
        'n_top_init_trees': ('Number of top initial trees. If not set, '
                             'program defaults will be used. See IQ-TREE '
                             'manual for details.'),
        'n_best_retain_trees': ('Number of best trees retained during '
                                'search. If not set, program defaults will '
                                'be used. See IQ-TREE manual for details.'),
        'n_iter': ('Fix number of iterations to stop. If not set, program '
                   'defaults will be used. See IQ-TREE manual for details.'),
        'stop_iter': ('Number of unsuccessful iterations to stop. If not '
                      'set, program defaults will be used. See IQ-TREE '
                      'manual for details.'),
        'perturb_nni_strength': ('Perturbation strength for randomized NNI. '
                                 'If not set, program defaults will be used. '
                                 'See IQ-TREE manual for details.'),
        'spr_radius': ('Radius for parsimony SPR search. If not set, '
                       'program defaults will be used. See IQ-TREE manual '
                       'for details.'),
        'allnni': ('Perform more thorough NNI search.'),
        'safe': ('Safe likelihood kernel to avoid numerical underflow')},
    output_descriptions={'tree': 'The resulting phylogenetic tree.'},
    name='Construct a phylogenetic tree with IQ-TREE.',
    description=('Construct a phylogenetic tree using IQ-TREE '
                 '(http://www.iqtree.org/) with automatic model selection.'),
    citations=[citations['Nguyen2015iqtree'],
               citations['Kalyaanamoorthy2017modelfinder']]
)

plugin.methods.register_function(
    function=q2_phylogeny.iqtree_ultrafast_bootstrap,
    inputs={'alignment': FeatureData[AlignedSequence]},
    parameters={
            'seed': Int,
            'n_cores': Int,
            'substitution_model': Str % Choices(_IQTREE_DNA_MODELS),
            'n_init_pars_trees': Int % Range(1, None),
            'n_top_init_trees': Int % Range(1, None),
            'n_best_retain_trees': Int % Range(1, None),
            'stop_iter': Int % Range(1, None),
            'perturb_nni_strength': Float % Range(0.01, 99),
            'spr_radius': Int % Range(1, None),
            'bootstrap_replicates': Int % Range(1000, None),
            'n_max_ufboot_iter': Int % Range(1, None),
            'n_ufboot_steps': Int % Range(1, None),
            'min_cor_ufboot': Float % Range(0.51, 0.99),
            'ep_break_ufboot': Float % Range(0.01, 0.99),
            'allnni': Bool,
            'safe': Bool},
    outputs=[('tree', Phylogeny[Unrooted])],
    input_descriptions={
        'alignment': ('Aligned sequences to be used for phylogenetic '
                      'reconstruction.'),
    },
    parameter_descriptions={
        'n_cores': ('The number of cores to use for parallel '
                    'processing. Use \'0\' to let IQ-TREE automatically '
                    'determine the optimal number of cores to use.'),
        'substitution_model': ('Model of Nucleotide Substitution.'
                               'If not provided, IQ-TREE will determine the '
                               'best fit substitution model automatically. '),
        'seed': ('Random number seed. If not set, program defaults will be '
                 'used. See IQ-TREE manual for details.'),
        'bootstrap_replicates': ('The number of bootstrap searches to '
                                 'perform. '),
        'n_init_pars_trees': ('Number of initial parsimony trees. If not '
                              'set, program defaults will be used. See '
                              'IQ-TREE manual for details.'),
        'n_top_init_trees': ('Number of top initial trees. If not set, '
                             'program defaults will be used. See IQ-TREE '
                             'manual for details.'),
        'n_best_retain_trees': ('Number of best trees retained during '
                                'search. If not set, program defaults will '
                                'be used. See IQ-TREE manual for details.'),
        'stop_iter': ('Number of unsuccessful iterations to stop. If not '
                      'set, program defaults will be used. See IQ-TREE '
                      'manual for details.'),
        'perturb_nni_strength': ('Perturbation strength for randomized NNI. '
                                 'If not set, program defaults will be used. '
                                 'See IQ-TREE manual for details.'),
        'spr_radius': ('Radius for parsimony SPR search. If not set, '
                       'program defaults will be used. See IQ-TREE manual '
                       'for details.'),
        'n_max_ufboot_iter': ('Maximum number of iterations. If not set, '
                              'program defaults will be used. See IQ-TREE '
                              'manual for details.'),
        'n_ufboot_steps': ('Number of iterations for UFBoot stopping rule. '
                           'If not set, program defaults will be used.'
                           'See IQ-TREE manual for details.'),
        'min_cor_ufboot': ('Minimum correlation coefficient. '
                           'If not set, program defaults will be used.'
                           'See IQ-TREE manual for details.'),
        'ep_break_ufboot': ('Epsilon value to break tie. If not set, program '
                            'defaults will be used. See IQ-TREE manual for '
                            'details.'),
        'allnni': ('Perform more thorough NNI search.'),
        'safe': ('Safe likelihood kernel to avoid numerical underflow')},
    output_descriptions={'tree': 'The resulting phylogenetic tree.'},
    name=('Construct a phylogenetic tree with IQ-TREE with bootstrap '
          'supports.'),
    description=('Construct a phylogenetic tree using IQ-TREE '
                 '(http://www.iqtree.org/) with automatic '
                 'model selection and bootstrap supports.'),
    citations=[citations['Nguyen2015iqtree'],
               citations['Kalyaanamoorthy2017modelfinder'],
               citations['Minh2013ultrafastbootstrap']]
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

plugin.pipelines.register_function(
    function=q2_phylogeny.phylogeny_from_mafft,
    inputs={
        'sequences': FeatureData[Sequence],
    },
    parameters={
        'n_threads': Int % Range(1, None),
        'max_gap_frequency': Float % Range(0, 1, inclusive_end=True),
        'min_conservation': Float % Range(0, 1, inclusive_end=True)
    },
    outputs=[
        ('alignment', FeatureData[AlignedSequence]),
        ('masked_alignment', FeatureData[AlignedSequence]),
        ('tree', Phylogeny[Unrooted]),
        ('rooted_tree', Phylogeny[Rooted]),
    ],
    input_descriptions={
        'sequences': 'The sequences to be used for creating a '
                     'fasttree based rooted phylogenetic tree.'
    },
    parameter_descriptions={
        'n_threads': 'The number of threads. (Use -1 to automatically use all '
                     'available cores) '
                     'This value is used when aligning the sequences and '
                     'creating the tree with fasttree.',
        'max_gap_frequency':  'The maximum relative frequency of gap '
                              'characters in a column for the column to be '
                              'retained. This relative frequency must be a '
                              'number between 0.0 and 1.0 (inclusive), where '
                              '0.0 retains only those columns without gap '
                              'characters, and 1.0 retains all columns '
                              'regardless of gap character frequency.'
                              'This value is used when masking the aligned '
                              'sequences.',
        'min_conservation':  'The minimum relative frequency '
                             'of at least one non-gap character in a '
                             'column for that column to be retained. This '
                             'relative frequency must be a number between 0.0 '
                             'and 1.0 (inclusive). For example, if a value of '
                             '0.4 is provided, a column will only be retained '
                             'if it contains at least one character that is '
                             'present in at least 40% of the sequences.'
                             'This value is used when masking the aligned '
                             'sequences.'
    },
    output_descriptions={
        'alignment': 'The aligned sequences.',
        'masked_alignment': 'The masked alignment.',
        'tree': 'The unrooted phylogenetic tree.',
        'rooted_tree': 'The rooted phylogenetic tree.',
    },
    name='Build a phylogenetic tree using fasttree and mafft alignment',
    description=('This pipeline generates a rooted phylogenetic tree by '
                 'wrapping up methods from q2-alignment and q2-phylogeny.'
                 'The methods include mafft and mask from q2-alignment and '
                 'fasttree and midpoint-root from q2-phylogeny. This pipeline '
                 'returns outputs produced by those methods.'
                 )
)
