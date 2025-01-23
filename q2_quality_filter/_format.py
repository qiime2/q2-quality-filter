# ----------------------------------------------------------------------------
# Copyright (c) 2017-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Union

import qiime2.plugin.model as model

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)


class QualityFilterStatsFmt(model.TextFileFormat):
    def _validate_(self, level):
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


_ReadDirectionTypes = Union[
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt
]


class _ReadDirectionUnion:
    format: _ReadDirectionTypes
