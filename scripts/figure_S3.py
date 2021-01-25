# figure_S3.py
# makes violin plots to compare sizes of networks pruned with and without
# export reactions

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# read input file
data = pd.read_csv('data/ab_5_2ins_5outs_sizes.csv', index_col = 0)

# make two plots, one for metabolite counts and one for reaction counts
fig, ax = plt.subplots(ncols = 2)
sns.violinplot(data = data, x = 'export', y = 'met_count', ax = ax[0])
sns.violinplot(data = data, x = 'export', y = 'rxn_count', ax = ax[1])

# make the plots aesthetically appealing
ax[0].set_xlabel('Export Reactions')
ax[0].set_ylabel('Number of Metabolites')
ax[0].text(
    # position the panel label above the plot on the left side
    -0.15, 1.05, 'a', va = 'top', ha = 'right',
    # no idea what this does; got it from StackOverflow
    transform = ax[0].transAxes,
    # make the letter large and bold
    fontsize = 14, fontweight = 'bold' 
)

ax[1].set_xlabel('Export Reactions')
ax[1].set_ylabel('Number of Reactions')
ax[1].text(
    # position the panel label above the plot on the left side
    -0.15, 1.05, 'b', va = 'top', ha = 'right',
    # no idea what this does; got it from StackOverflow
    transform = ax[1].transAxes,
    # make the letter large and bold
    fontsize = 14, fontweight = 'bold' 
)

# make sure the two panels don't overlap or anything else unsightly
plt.tight_layout()
plt.savefig('data/figure_S2.png', dpi = 600)
