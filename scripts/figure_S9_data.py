# figure_S9_data.py
'''
Make the (2,5) universal string chemistry network and two COBRApy models of it
with and without export reactions, then prune each of those networks using the 
biomass-impact pruner with 100 different biomass reactions and 100 different 
environments per biomass reaction. Save the reaction-inclusion vectors and 
biomass fluxes for the pruned networks in files that figure_S9_plot.py will 
read to make plots
'''

import sys
import string_chem_net as scn
import pandas as pd
import itertools as it
import multiprocessing as mp

# prune one biomass reaction in the specified number of environments
def prune_many_times(arglist):
    # probably a more elegant way to do this but I'm currently new to mp.map()
    full_model, ins, outs, envs = arglist
    # start by making a copy of the original model so we don't have to remove
    # the biomass reaction each time
    model = full_model.copy()
    # add a biomass reaction and set it as the objective
    bm_rxn = scn.choose_bm_mets(outs, model)
    model.objective = bm_rxn
    # keep lists of the environments used, the reaction-inclusion vectors of
    # the pruned networks and the growth rates on the pruned networks
    food_mets = list()
    rxn_incl_vecs = list()
    pruned_growths = list()
    # use a while loop and not a for loop so we can go back on occasion
    i = 0
    # counter for how many times it had to reselct the environment to get a
    # feasible solution with the full network
    j = 0
    while i < envs:
        i +=  1 
        # remove existing input reactions
        in_rxns = [rxn for rxn in model.boundary if rxn.id.startswith('->')]
        model.remove_reactions(in_rxns)
        # choose new input reactions
        scn.choose_inputs(ins, model, bm_rxn)
        in_rxns = [rxn for rxn in model.boundary if rxn.id.startswith('->')]
        foods_string = ' '.join([
            # getting the metabolite IDs out of a reaction is annoying
            list(rxn.metabolites.keys())[0].id for rxn in in_rxns
        ])
        # see if this choice of metabolites can produce the biomass on this network
        solution = model.optimize()
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        if solution.status == 'infeasible' or bm_rxn_flux < 1e-10:
            # redo this iteration of the loop
            i -= 1
            # increment the counter of redos
            j += 1 
            continue
        # record the metabolites that worked and prune the network
        else:
            if i % 100 == 0:
                print(f'On environment {i}')
            # reset the reselection counter
            j = 0
            # get the list of food source metabolites
            food_mets.append('-'.join([
                met.id 
                for rxn in in_rxns 
                for met in rxn.metabolites
            ]))
            # prune the network
            pruned_net = scn.bm_impact_prune(model, bm_rxn)
            rxn_incl = scn.make_rxn_incl(model, pruned_net)
            rxn_incl_vecs.append(rxn_incl)
            # get the growth rate on the pruned network
            solution = pruned_net.optimize()
            pruned_growth = solution.fluxes.get(key = bm_rxn.id)
            pruned_growths.append(pruned_growth)

    # make a dataframe out of the lists and add it to the larger dataframe
    data = pd.DataFrame(list(zip(
        food_mets, rxn_incl_vecs, pruned_growths
    )))
    data.columns = ['env','rxn_incl', 'growth']
    # add a column with the biomass components
    data['biomass'] = list(it.repeat(
        '-'.join([met.id for met in bm_rxn.metabolites]),
        len(food_mets)
    ))
    # reorder columns 
    data = data[['biomass', 'env', 'rxn_incl', 'growth']]
    return(data)

try:
    threads = int(sys.argv[1])
except IndexError:
    sys.exit('Specify a number of threads to run on')

# define parameters
monos = 'ab'
max_pol = 5
ins = 2 # number of input metabolites in each set
envs = 100 # number of different sets of ins to choose
outs = 5 # number of biomass precursors in each biomass reaction
orgs = 100 # number of different sets of outs to choose

# create the universal network
SCN = scn.CreateNetwork(monos, max_pol)
full_model = scn.make_cobra_model(
    SCN.met_list, 
    SCN.rxn_list, 
    allow_export = True
)

# just in case we're trying a particularly large number of biomass reactions,
# run the function in parallel since each biomass reaction can be handled
# completely independently of the others
pool = mp.Pool(threads)
data_bits = pool.map(
    prune_many_times,
    # same arguments every time for orgs times
    [[full_model, ins, outs, envs] for bm in range(orgs)]
)
# concatenate all the dataframes and write to output
all_data = pd.concat(data_bits)
all_data.to_csv(f'data/figure_S9_data.csv')
