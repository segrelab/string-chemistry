# var_bm_coefs.py
'''
Make a network with a particular set of input reactions and a biomass reaction
and vary the coefficients in that biomass reaction while keeping the identities
of the nutrients and biomass precursors constant
'''

import string_chem_net as scn
import random
import pandas as pd

monos = 'ab' # characters to use as monomers
max_len = 5  # maximum length of each string chemical
ins = 2      # number of food sources / environmental nutrients
outs = 5     # number of biomass precursors
combos = 100 # number of times to perturb biomass coefficients

SCN = scn.CreateNetwork(monos, max_len)
full_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)

# make some reactions to import nutrient metabolites and represent creation of
# biomass precursors
bm_rxn = scn.choose_bm_mets(outs, full_model)
scn.choose_inputs(ins, full_model, bm_rxn)
full_model.objective = bm_rxn

# while that combination of nutrients and biomass precursors might result in a
# feasible network, it frequently doesn't, so here's the way I make sure that
# the pair I've chosen is viable before doing anything else
solution = full_model.optimize()
# can't just check solution.status because sometimes it's feasible but the
# flux through the biomass reaction is vanishingly small
bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
while solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
    # if the solution isn't feasible, pick a different environment
    in_rxns = [
        # don't want to remove all boundary reactions because that would
        # also remove all of the export reactions
        rxn for rxn in full_model.boundary if rxn.id.startswith('->')
    ]
    full_model.remove_reactions(in_rxns)
    scn.choose_inputs(int(ins), full_model, bm_rxn)
    solution = full_model.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)

# randomly reassign coefficients in biomass reaction to numbers between 1 and
# 10 a bunch of times
rxn_incl_vecs = list()
for combo in range(combos):
    print(f'On round {combo} of {combos}')
    # copy the biomass reaction so the changes to coefficients don't add up
    new_bm = bm_rxn.copy()
    new_bm.add_metabolites(
        {m: -random.randint(1, 9) for m in bm_rxn.metabolites}
    )
    # also copy the model and replace the old biomass reaction with this one
    model = full_model.copy()
    print('Removing previous biomass reaction')
    model.remove_reactions([bm_rxn])
    print('Adding new biomass reaction')
    model.add_reaction(new_bm)
    model.objective = new_bm
    # now prune
    pruned_model = scn.min_flux_prune(model, new_bm)
    # record which reactions were kept
    rxn_incl_vecs.append(scn.make_rxn_incl(full_model, pruned_model))

data = pd.DataFrame({'Reactions' : rxn_incl_vecs})
data.to_csv('data/varied_bm_coefs.csv', index = False)
