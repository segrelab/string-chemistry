# growth_data_generator.py
# does min flux pruning on a single universal scale network and then finds 
# growth in a bunch of different environments

import sys
import string_chem_net as scn
import random
import cobra
import pandas as pd

# get command-line arguments
try:
    (monos, max_pol, ins, outs, env_count, bm_count) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments:\nmonomers\nmax polymer length\nnumber of food ' +
        'sources in each environment\nnumber of biomass precursors\n' +
        'number of environments\nnumber of different biomass reactions ' +
        '(i.e. different pruned networks)' 
    )

# create the universal network
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))

# store data as list of lists with biomass reaction, environment and growth
env_growth_lists = list()

# loop over number of different biomass reactions to use
for bm_trial in range(int(bm_count)):
    # make sure that there's at least one feasible solution before trying to prune
    # have to have some inputs to test this biomass reaction
    scn.choose_inputs(int(ins), cobra_model)
    bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
    cobra_model.objective = bm_rxn

    solution = cobra_model.optimize()
    while solution.status == 'infeasible' or (solution.fluxes == 0).all():
        # remove existing biomass and input reactions
        cobra_model.remove_reactions([bm_rxn])
        cobra_model.remove_reactions(cobra_model.boundary)
        # choose new ones and see if those yield a solvable network
        bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
        scn.choose_inputs(int(ins), cobra_model, bm_rxn)
        cobra_model.objective = bm_rxn
        solution = cobra_model.optimize()

    # now that we've chosen a biomass reaction, make a string describing it
    bm_string = '+'.join([met.id for met in bm_rxn.metabolites])

    # do min flux pruning
    min_flux_pruned = scn.min_flux_prune(cobra_model)

    # get growth flux in as many environments as specified
    for env in range(int(env_count)):
        # start by removing existing input reactions
        min_flux_pruned.remove_reactions(min_flux_pruned.boundary)
        # get random new environment, making sure none of the food sources are also
        # biomass precursors
        scn.choose_inputs(int(ins), min_flux_pruned, bm_rxn)
        # find growth in this environment
        solution = min_flux_pruned.optimize()
        growth = solution.objective_value
        # prepare string describing environment
        in_rxns = min_flux_pruned.boundary
        in_mets = [met.id for rxn in in_rxns for met in rxn.metabolites]
        env_string = '+'.join(in_mets)
        # add to list of lists
        env_growth_lists.append([bm_string, env_string, growth])

    # remove input reactions in preparation for next round of pruning
    cobra_model.remove_reactions(cobra_model.boundary)

# make the list of lists into a DataFrame
env_growth_df = pd.DataFrame(
    env_growth_lists,
    columns = ['biomass', 'environment', 'growth']
)

env_growth_df.to_csv(
    f'data/{monos}_{max_pol}_{ins}ins_{outs}outs_{env_count}envs_{bm_count}orgs.csv'
)
