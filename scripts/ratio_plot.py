# ratio_plot.py
# makes a histogram showing reaction-to-metabolite ratios for many networks
# pruned from the same-sized universal network

import sys
import pandas as pd
from matplotlib import pyplot as plt

in_file = sys.argv[1]

# parse title to figure out how large the universal network is
bits = in_file.lstrip('data').rstrip('.csv').split('_')

data = pd.read_csv(in_file, header=None)
data.hist()
plt.xlabel('Number of reactions/Number of metabolites')
plt.ylabel('Frequency of networks with this ratio')
plt.title(
    f'Reaction/Metabolite Ratios for Universal Network {bits[0][1:]} {bits[1]}'
)
plt.show()
