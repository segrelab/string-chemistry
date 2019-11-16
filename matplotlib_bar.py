# matplotlib_bar.py
# make a bar chart showing how many times each network was found with the 
# networks ordered by commonness and colored by reaction count
# start by making a column ranking the networks by frequency

import sys
import pandas as pd
import matplotlib.pyplot as plt

bitstring_df = pd.read_csv(sys.argv[1])

bitstring_df['freq_rank'] = bitstring_df['occurrences'].rank(
    ascending = False
)
# then sort the dataframe by those ranks
bitstring_df.sort_values('freq_rank', inplace = True)

# then make the bar plot
bitstring_df.plot.bar(
    x = 'freq_rank', y = 'occurrences', legend = False
)
# print the number of reactions in each network on top of each bar
xlocs, xlabs = plt.xticks()
for idx,v in enumerate(bitstring_df['occurrences']):
    plt.text(
        xlocs[idx] - 0.1, # x coord of label
        v + 0.025, # y coord of label
        str(bitstring_df['rxn_count'].iloc[idx]) # label
    )
plt.title(f'Frequency of Observing Each Network After {reps} Prunes')
plt.xlabel('Rank of Network by Number of Times It Was Seen')
plt.ylabel('Number of Times Network Was Seen')
plt.show()
