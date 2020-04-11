# regression_figure.py
# plots describing the regression model for predicting growth from envrionment

import pandas as pd
import plotly.express as px

env_growth_df = pd.read_csv(
    'data/ab_5_5ins_3outs_10000envs_10orgs.csv',
    index_col = 0
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
