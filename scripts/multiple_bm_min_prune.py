# multiple_bm_min_prune.py
# makes one network, picks one environment, and prunes the universal network
# with that environment and many randomly-chosen biomass reactions using the
# minimum flux pruner

import sys
import string_chem_net as scn
import random
import cobra

# count number of 1s in a bitstring
def count_bitstring(bitstring):
    count = 0
    for bit in [int(bit) for bit in list(bitstring)]:
        if bit == 1:
            count += 1
    return(count)

# get command-line arguments
try:
    (monos, max_pol, ins, outs, reps) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments:\nmonomers\nmax polymer length\nnumber of food ' +
        'sources in environment\nnumber of biomass precursors\nnumber of ' + 
        'different biomass reactions to use'
    )

# create the reference network and pick an environment
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.choose_inputs(int(ins), cobra_model)

# reaction bitstrings as keys, biomass components as values
pruned_nets = dict()
i = 0
while i < int(reps):
    # pick a new biomass reaction and set it as the objective
    bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
    cobra_model.objective = bm_rxn
    # see if these biomass precursors can be produced on this environment
    # before we bother pruning
    solution = cobra_model.optimize()
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    if solution.status != 'infeasible' and bm_rxn_flux > 10e-10:
        # only increment the counter if we've chosen a usable biomass reaction
        i += 1
        print(f'On biomass reaction {i}')
        # run the minimum flux pruner
        pruned_net = scn.min_flux_prune(cobra_model, bm_rxn)
        # remove the biomass reaction before making the reaction bitstring
        # can't just remove the reaction because somehow in the pruning process
        # the biomass reaction became different in some subtle way
        pruned_bm_rxn = pruned_net.reactions.get_by_id(bm_rxn.id)
        pruned_net.remove_reactions([pruned_bm_rxn])
        bitstring = scn.make_bitstring(cobra_model, pruned_net)
        pruned_nets[bitstring] = [met.id for met in bm_rxn.metabolites]
    # remove this biomass reaction from the full network regardless of whether
    # it worked or not
    cobra_model.remove_reactions([bm_rxn])

# print a bunch of info but also write it out to a tsv
with open(
        f'../data/{monos}_{max_pol}_{ins}env_{reps}x{outs}outs.tsv', 'w'
    ) as out:
    out.write('bitstring\trxn_count\tbiomass\n')
    for network in pruned_nets.keys():
        rxn_count = count_bitstring(network)
        out_row = '\t'.join([
            network, str(rxn_count), ','.join(pruned_nets[network])
        ])
        out.write(out_row + '\n')
