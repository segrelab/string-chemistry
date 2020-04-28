# multiple_env_random_prune.py
# makes one network, picks one biomass reaction, and then runs random_prune
# for several different choices of n food source metabolites

import sys
import string_chem_net as scn
import random
import pandas as pd
from scipy.stats import chisquare

# count number of 1s in a bitstring
def count_bitstring(bitstring):
    count = 0
    for bit in [int(bit) for bit in list(bitstring)]:
        if bit == 1:
            count += 1
    return(count)

# get command-line arguments
try:
    (monos, max_pol, ins, envs, outs, reps) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments: monomers, max polymer length, number of food ' +
    'sources in each environment, number of environments, number of ' + 
    'biomass precursors, number of times to prune each network')

# create the reference network and pick some food sources and a biomass rxn
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.choose_inputs(int(ins), cobra_model)
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
print(f'Biomass reaction: {bm_rxn.id}')
cobra_model.objective = bm_rxn

i = 0
# record all of the food metabolites that were used
food_mets = list()
# record the p-values from chi-squared tests against the null hypothesis that
# each network was equally likely to be observed
p_vals = list()
# use a while loop and not a for loop so we can go back on occasion
while i < int(envs):
    i = i + 1
    in_rxns = [
        rxn for rxn in cobra_model.boundary
        # only get input reactions
        if rxn.id.startswith('->')
    ]
    foods_string = ' '.join([rxn.metabolites[0] for rxn in in_rxns])
    print(f'Food source group {i}: {foods_string}')
    # choose some new food sources (remove existing ones)
    cobra_model.remove_reactions(in_rxns)
    scn.choose_inputs(int(ins), cobra_model, bm_rxn)
    # see if this choice of metabolites can produce the biomass on this network
    solution = cobra_model.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    if solution.status == 'infeasible' or bm_rxn_flux < 1e-10:
        print('There were no feasible solutions; reselecting food sources.')
        # redo this iteration of the loop
        i = i - 1
        continue
    # record the metabolites that worked before proceeding
    else:
        food_mets.append([met.id for met in cobra_model.boundary])

    # will hold bitstrings of all unique reactions and the count of times each
    # one came up
    random_pruned_dict = dict()
    # will hold all the unique networks found by random_prune after reps runs
    random_pruned_nets = list()
    for j in range(1, int(reps)):
        if j % 100 == 0:
            print(f'Pruned {j} times')
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
    bitstring_df.to_csv(
        f'data/multiple_env_random_prune_{monos}_{max_pol}_{ins}ins_{outs}outs_{i}of{envs}.csv'
    )
    chisq_results = chisquare(bitstring_df.occurrences)
    p_vals.append(chisq_results[1])

# write the list of food sources that were tried to a file
with open(f'data/multiple_env_random_prune{monos}_{max_pol}_{ins}ins_' + 
    f'{outs}outs_foods.csv', 'w') as out:
    # once all the lines are created, concatenate them with newline characters
    output = '\n'.join([
        # add the appropriate p-value to the concatenated food lists
        ','.join([p_val].append(
            # flatten all of the food source lists
            ['-'.join(sublist) for sublist in food_mets]
        )) for p_val in p_vals
    ])
    out.write(output)
