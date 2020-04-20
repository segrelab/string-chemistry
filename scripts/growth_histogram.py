# growth_histogram.py
# histogram of growth fluxes for many networks pruned with different biomass
# reactions in different environments

import pandas as pd
import plotly.express as px

env_growth_df = pd.read_csv(
    'data/ab_5_5ins_3outs_10000envs_10bms.tsv',
    index_col = 0,
    sep = '\t'
) 

fig = px.histogram(
    env_growth_df,
    x = 'growth',
    color = 'biomass',
    nbins = 50
)

fig.update_layout(
    barmode = 'group', showlegend = True, legend_orientation = 'h')
fig.show()
