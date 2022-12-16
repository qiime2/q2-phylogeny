# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase


class TestUsageExample(TestPluginBase):
    package = 'q2_phylogeny.tests'

    def test_usage_examples(self):
        self.execute_examples()
