# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np
import qiime2

from .plugin_setup import plugin
from ._format import QualityFilterStatsFmt


@plugin.register_transformer
def _1(data: pd.DataFrame) -> QualityFilterStatsFmt:
    ff = QualityFilterStatsFmt()
    data.to_csv(str(ff))
    return ff


_stats_column_dtypes = {
    'sample-id': np.str,
    'total-input-reads': np.number,
    'total-retained-reads': np.number,
    'reads-truncated': np.number,
    'reads-too-short-after-truncation': np.number,
    'reads-exceeding-maximum-ambiguous-bases': np.number
}


def _stats_to_df(ff):
    # https://github.com/pandas-dev/pandas/issues/9435
    df = pd.read_csv(str(ff), dtype=_stats_column_dtypes)
    df.set_index('sample-id', inplace=True)
    return df


@plugin.register_transformer
def _2(ff: QualityFilterStatsFmt) -> pd.DataFrame:
    return _stats_to_df(ff)


@plugin.register_transformer
def _3(ff: QualityFilterStatsFmt) -> qiime2.Metadata:
    return qiime2.Metadata(_stats_to_df(ff))
