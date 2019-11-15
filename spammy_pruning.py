# spammy_pruning.py
# makes one network, picks one biomass reaction, and then runs random_prune
# for several different choices of n food source metabolites

import sys
import string_chem_net as scn
import random
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# count number of 1s in a bitstring
def count_bitstring(bitstring):
    count = 0
    for bit in [int(bit) for bit in list(bitstring)]:
        if bit == 1:
            count += 1
    return(count)

# get command-line arguments
try:
    (monos, max_pol, ins, outs, reps, big_reps) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments: monomers, max polymer length, number of food sources, \
number of biomass precursors, number of times to prune each network, number of \
times to reselect food sources.')

# create the reference network and pick some food sources and a biomass rxn
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))
scn.choose_inputs(cobra_model, int(ins))
bm_rxn = scn.choose_bm_mets(cobra_model, int(outs))
cobra_model.objective = bm_rxn

for i in range(int(big_reps)):
    print(f'On set {i} of food sources')
    # choose some new food sources (remove existing ones)
    cobra_model.remove_reactions(cobra_model.boundary)
    scn.choose_inputs(cobra_model, int(ins))
    # see if this choice of metabolites can produce the biomass on this network
    solution = cobra_model.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    if solution.status == 'infeasible' or bm_rxn_flux < 1e-10:
        print('There were no feasible solutions.')
        continue
    # if there's at least one feasible solution, proceed with the pruning

    # will hold bitstrings of all unique reactions and the count of times each
    # one came up
    random_pruned_dict = dict()
    # will hold all the unique networks found by random_prune after reps runs
    random_pruned_nets = list()
    for j in range(1, int(reps)):
        if j % 100 == 0:
            print(j)
        pruned_net = scn.random_prune(cobra_model, bm_rxn)
        # in order to know whether we've seen this model before, we can't just 
        # compare models, since no two models are ever 'equal', so we'll compare
        # reaction presence bitstrings. 
        # We also want to keep track of how many times we see each model, so we 
        # will make a dict with the bitstrings as keys
        bitstring = scn.make_bitstring(cobra_model, pruned_net)
        if bitstring not in random_pruned_dict.keys():
            # make sure all reaction lists are sorted so that all isomorphic
            # networks have the same reaction list
            random_pruned_dict[bitstring] = 1
            random_pruned_nets.append(pruned_net)
        else:
            # if we already found this network once, then increment the
            # appropriate counter by 1
            random_pruned_dict[bitstring] += 1

    # make a dataframe of all the bitstrings we've generated
    # turn the dictionary into a list of lists before coercion
    bitstring_df = pd.DataFrame(list(map(list, random_pruned_dict.items())))
    bitstring_df.columns = ['bitstring', 'occurrences']
    # add in a column for the number of reactions in each network 
    bitstring_df['rxn_count'] = list(map(
        count_bitstring, bitstring_df.bitstring
    ))

    # make a bar chart showing how many times each network was found with the 
    # networks ordered by commonness and colored by reaction count
    # start by making a column ranking the networks by frequency
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
    plt.savefig(
        f'{monos}_{max_pol}_{ins}ins_{outs}outs_{i}of{big_reps}reps.png'
    )
