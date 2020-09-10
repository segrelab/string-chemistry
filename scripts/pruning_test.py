# pruning_test.py
# prune one network once and get detailed information on the output (mostly for
# debugging purposes)

import sys
import string_chem_net as scn

# get specs governing network size from command-line arguments
try:
    (monos, max_pol, ins, outs) = sys.argv[1:]
except ValueError:
    sys.exit(
        'Arguments:\nmonomers\nmax polymer length\nnumber of food ' + 
        'sources\nnumber of biomass precursors'
    )

SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(
    SCN.met_list,
    SCN.rxn_list,
    allow_export = True
)
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
scn.choose_inputs(int(ins), cobra_model, bm_rxn)
cobra_model.objective = bm_rxn

# see if this environment supports growth with this biomass reaction
solution = cobra_model.optimize()
while solution.status == 'infeasible' or (solution.fluxes == 0).all():
    # if not, remove existing biomass and input reactions
    cobra_model.remove_reactions([bm_rxn])
    in_rxns = [rxn for rxn in cobra_model.boundary if rxn.id.startswith('->')]
    cobra_model.remove_reactions(in_rxns)
    # choose new ones and see if those yield a solvable network
    bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
    scn.choose_inputs(int(ins), cobra_model, bm_rxn)
    cobra_model.objective = bm_rxn
    solution = cobra_model.optimize()

# prune the network and print some summary statistics about it
pruned = scn.min_flux_prune(cobra_model, bm_rxn)
# metabolites aren't automatically removed when all their associated reactions
# are removed, so count the number of metabolites left with reactions
met_count = len([m for m in pruned.metabolites if len(m.reactions) > 0])
ratio = len(pruned.reactions) / met_count
print(f'Number of metabolites: {met_count}')
print(f'Number of reactions: {len(pruned.reactions)}')
print(f'Reaction-to-metabolite ratio: {ratio}')
