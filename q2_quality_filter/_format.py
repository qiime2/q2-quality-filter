# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model


class QualityFilterStatsFmt(model.TextFileFormat):
    def sniff(self):
        line = open(str(self)).readline()
        hdr = line.strip().split(',')
        expected = ['sample-id', 'total-input-reads',
                    'total-retained-reads',
                    'reads-truncated',
                    'reads-too-short-after-truncation',
                    'reads-exceeding-maximum-ambiguous-bases']
        return hdr == expected


QualityFilterStatsDirFmt = model.SingleFileDirectoryFormat(
    'QualityFilterStatsDirFmt', 'stats.csv', QualityFilterStatsFmt)
