#!/usr/bin/env

import cufflinks as cf
import plotly.express as px
import pandas as pd
import numpy as np

filepath="/home/dermoi/Schreibtisch/Segre_Rotation/"
duke=pd.read_csv(filepath+"bitstrings.csv")
fig=px.bar(duke, x='rxn_count',y='occurences',color='bitstring')

fig.update_layout(barmode='group',showlegend=True, legend_orientation='h')
fig.show()
