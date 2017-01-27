# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import pandas as pd
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_quality_filter', 'assets')


def visualize_stats(output_dir: str, filter_stats: pd.DataFrame) -> None:
    sums = filter_stats.sum()
    sums.name = 'Totals'
    filter_stats = filter_stats.append(sums)

    filter_stats.sort_values('total-input-reads', inplace=True,
                             ascending=False)

    total_retained = filter_stats['total-retained-reads']
    total_input = filter_stats['total-input-reads']
    filter_stats['fraction-retained'] = total_retained / total_input

    # reorder such that retained fraction follows total-input-reads and
    # total-retained-reads
    columns = list(filter_stats.columns)[:-1]
    columns.insert(2, 'fraction-retained')
    filter_stats = filter_stats[columns]

    html = filter_stats.to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    index = os.path.join(TEMPLATES, 'index.html')
    context = {
        'result': html
    }

    q2templates.render(index, output_dir, context=context)
