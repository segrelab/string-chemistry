# pruning_test.py
# just runs one or both of the pruning algorithms once to see what they're
# doing for mostly debugging purposes

import string_chem_net as scn

# generate a network of some arbitrary size
monos = 'ab'
max_len = 5
ins = 2
outs = 5

SCN = scn.CreateNetwork('ab', 5)
cobra_model = scn.make_cobra_model(
    SCN.met_list,
    SCN.rxn_list,
    allow_export = True
)
bm_rxn = scn.choose_bm_mets(outs, cobra_model)
scn.choose_inputs(ins, cobra_model, bm_rxn)
cobra_model.objective = bm_rxn

# see if this environment supports growth with this biomass reaction
solution = cobra_model.optimize()
while solution.status == 'infeasible' or (solution.fluxes == 0).all():
    # if not, remove existing biomass and input reactions
    cobra_model.remove_reactions([bm_rxn])
    in_rxns = [rxn for rxn in cobra_model.boundary if rxn.id.startswith('->')]
    cobra_model.remove_reactions(in_rxns)
    # choose new ones and see if those yield a solvable network
    bm_rxn = scn.choose_bm_mets(outs, cobra_model)
    scn.choose_inputs(ins, cobra_model, bm_rxn)
    cobra_model.objective = bm_rxn
    solution = cobra_model.optimize()

# prune using the random pruner
pruned = scn.random_prune(cobra_model, bm_rxn)
print(scn.make_rxn_incl(cobra_model, pruned))
pruned = scn.random_prune(cobra_model, bm_rxn)
print(scn.make_rxn_incl(cobra_model, pruned))
pruned = scn.random_prune(cobra_model, bm_rxn)
print(scn.make_rxn_incl(cobra_model, pruned))
pruned = scn.random_prune(cobra_model, bm_rxn)
print(scn.make_rxn_incl(cobra_model, pruned))
pruned = scn.random_prune(cobra_model, bm_rxn)
print(scn.make_rxn_incl(cobra_model, pruned))
pruned = scn.random_prune(cobra_model, bm_rxn)
print(scn.make_rxn_incl(cobra_model, pruned))

