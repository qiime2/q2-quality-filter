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
    filter_stats.sort_values('total-input-reads', inplace=True,
                             ascending=False)

    html = filter_stats.to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    index = os.path.join(TEMPLATES, 'index.html')
    context = {
        'result': html
    }

    q2templates.render(index, output_dir, context=context)
