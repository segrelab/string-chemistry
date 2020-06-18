# multiple_env_fba.py
# runs normal FBA many times on the same universal network with different
# biomass reactions and environments (but such that each biomass reaction is
# run with the same environments)

import sys
import string_chem_net as scn
import pandas as pd
import itertools as it

# get command-line arguments
try:
    (monos, max_pol, ins, envs, outs, orgs, export) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments:\nmonomers\nmax polymer length\n' +
        'number of food sources\nnumber of times to reselect food sources\n' +
        'number of biomass precursors\nnumber of times to reselect biomass\n' +
        'should there be an export reaction for every metabolite? (yes/no)')

if export == 'yes':
    allow_export = True
elif export == 'no':
    allow_export = False
else:
    sys.exit('The last argument must be either "yes" or "no"')

# create the universal network
SCN = scn.CreateNetwork(monos, int(max_pol))
untouched_model = scn.make_cobra_model(
    SCN.met_list, 
    SCN.rxn_list, 
    allow_export = allow_export
)
# make a dataframe to store information about the pruned networks
all_data = pd.DataFrame(columns = ['env', 'rxn_incl', 'biomass'])
# loop over the different biomass reactions
for bm in range(int(orgs)):
    print(f'On biomass reaction {bm}')
    # start by making a copy of the original model so we don't have to remove
    # the biomass reaction each time
    model = untouched_model.copy()
    # add a biomass reaction and set it as the objective
    bm_rxn = scn.choose_bm_mets(int(outs), model)
    model.objective = bm_rxn
    # keep lists of the environments used and the reaction-inclusion vectors of
    # the pruned networks
    food_mets = list()
    rxn_incl_vecs = list()
    # use a while loop and not a for loop so we can go back on occasion
    i = 0
    # counter for how many times it had to reselct the environment to get a
    # feasible solution with the full network
    j = 0
    while i < int(envs):
        i +=  1 
        # remove existing input reactions
        in_rxns = [rxn for rxn in model.boundary if rxn.id.startswith('->')]
        model.remove_reactions(in_rxns)
        # choose new input reactions
        scn.choose_inputs(int(ins), model, bm_rxn)
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
        # record these metabolites and the reactions that had flux
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
            # remove all reactions without flux
            no_flux_rxns = solution.fluxes[solution.fluxes == 0].index
            flux_only = model.copy()
            flux_only.remove_reactions(no_flux_rxns)
            rxn_incl = scn.make_rxn_incl(model, flux_only)
            rxn_incl_vecs.append(rxn_incl)

    # make a dataframe out of the two lists and add it to the larger dataframe
    more_data = pd.DataFrame(list(zip(food_mets, rxn_incl_vecs)))
    more_data.columns = ['env','rxn_incl']
    # add a column with the biomass components
    more_data['biomass'] = list(it.repeat(
        '-'.join([met.id for met in bm_rxn.metabolites]),
        len(food_mets)
    ))
    all_data = all_data.append(more_data)

all_data.to_csv(
    f'data/multiple_env_fba_{monos}_{max_pol}_{ins}ins_{envs}envs_' +
    f'{outs}outs_{orgs}orgs_{export}exp.csv'
)
