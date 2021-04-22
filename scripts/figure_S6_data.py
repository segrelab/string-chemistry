# figure_S6_data.py
'''
Similar to figure_4_data.py, but in addition to varying which reactions
produce nutrients, also vary the stoichiometric coefficients of the biomass
reactions
'''

import string_chem_net as scn
import multiprocessing as mp
import random
import pandas as pd
import cobra
import itertools as it

def prune_many_times(arglist):
    '''
    Given:
    - COBRApy model representing a full/complete/un-pruned string chemistry
    - A number of nutrient sources
    - A number of biomass precursors
    - A number of different sets of nutrient sources
    - A number of times to change stoichiometric coefficients in the biomass
      reaction
    Do:
    - Create a biomass reaction with the designated number of reactants
    - Create the designated number of variants on that reaction with randomized 
      stoichiometric coefficients (each reactant's coefficient is assigned to
      a random integer between 1 and 10)
    - Choose the designated number of sets of the designated number of nutrient
      sources
    - Prune the network once for each combination of biomass reaction and
      set of nutrients
    Return a Dataframe containing:
    - Binary vector indicating which reactions were kept in each pruned network
    - List of biomass precursors used
    - List of nutrient sources used
    - Maximum achievable flux through biomass reaction
    '''
    # probably a more elegant way to do this but I'm currently new to mp.map()
    (full_model, ins, outs, envs, combos) = arglist
    # add a biomass reaction but remove it from the model immediately
    bm_rxn = scn.choose_bm_mets(outs, full_model)
    full_model.remove_reactions([bm_rxn])
    # loop over vaariants of the biomass reaction with different coefficients
    # but identical reactants
    for combo in range(combos):
        if (combo + 1) % 10 == 0:
            print(f'On coefficient set {combo+1} of {combos}')
        # make a new biomass reaction
        new_bm = cobra.Reaction('varied_bm_rxn')
        new_bm.add_metabolites(
            {m: -random.randint(1, 10) for m in bm_rxn.metabolites}
        )
        # make a copy of the model before adding the new biomass reaction
        model = full_model.copy()
        model.add_reaction(new_bm)
        model.objective = new_bm
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
            scn.choose_inputs(ins, model, new_bm)
            in_rxns = [rxn for rxn in model.boundary if rxn.id.startswith('->')]
            foods_string = ' '.join([
                # getting the metabolite IDs out of a reaction is annoying
                list(rxn.metabolites.keys())[0].id for rxn in in_rxns
            ])
            # see if this choice of metabolites can produce the biomass on this network
            solution = model.optimize()
            bm_rxn_flux = solution.fluxes.get(key = new_bm.id)
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
                pruned_net = scn.min_flux_prune(model, new_bm)
                rxn_incl = scn.make_rxn_incl(model, pruned_net)
                rxn_incl_vecs.append(rxn_incl)
                # get the growth rate on the pruned network
                solution = pruned_net.optimize()
                pruned_growth = solution.fluxes.get(key = new_bm.id)
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

monos = 'ab' # characters to use as monomers
max_len = 5  # maximum length of each string chemical
ins = 2      # number of food sources / environmental nutrients
envs = 50    # number of different sets of food sources per biomass reaction
outs = 5     # number of biomass precursors
orgs = 10    # number of different biomass reactions
combos = 50  # number of times to perturb coefficients per biomass reaction
threads = 4  # threads to use when pruning in parallel

SCN = scn.CreateNetwork(monos, max_len)
full_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)

# prune network using many biomass reactions and environments
pool = mp.Pool(1)
data_bits = pool.map(
    prune_many_times,
    # same arguments every time for orgs times
    [[full_model, ins, outs, envs, combos] for bm in range(orgs)]
)
data = pd.concat(data_bits)
data.to_csv('data/figure_S6_data.csv', index = False)
