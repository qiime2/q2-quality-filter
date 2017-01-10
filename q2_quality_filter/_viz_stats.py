import os
import pkg_resources

import pandas as pd
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_quality_filter', 'assets')

def visualize_stats(output_dir: str, data: pd.DataFrame) -> None:
    data.sort_values('total-input-reads', inplace=True, ascending=False)
    
    html = data.to_html(classes='table table-striped table-hover')
    html = html.replace('border="1"', 'border="0"')
    index = os.path.join(TEMPLATES, 'index.html')
    context = {
        'result': html
    }

    q2templates.render(index, output_dir, context=context)
