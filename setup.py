# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-phylogeny",
    # TODO stop duplicating version string
    version="0.0.7.dev0",
    packages=find_packages(),
    install_requires=['scikit-bio', 'qiime >= 2.0.6', 'q2-types >= 0.0.6'],
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Create and work with Phylogenetic trees in QIIME 2.",
    license="BSD",
    url="http://www.qiime.org",
    entry_points={
        'qiime.plugins': ['q2-phylogeny=q2_phylogeny.plugin_setup:plugin']
    },
    package_data={'q2_phylogeny.tests': ['data/*']}
)
