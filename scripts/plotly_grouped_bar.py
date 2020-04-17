# plotly_grouped_bar.py
# uses plotly to make a grouped bar chart

import sys
import plotly.express as px
import pandas as pd
import numpy as np

data = pd.read_csv(sys.argv[1])
fig = px.bar(data, x = 'rxn_count', y = 'occurrences', color = 'bitstring')

fig.update_layout(
    barmode = 'group', showlegend = True, legend_orientation = 'h')
fig.show()
