# exhaustive_prune.py
# an attempt at finding every possible subnetwork that can produce the given
# biomass components with the given food sources; probably only usable on very
# small networks of all possible reactions

import sys
import string_chem_net as scn
import random
import cobra

# given a network and the food sources it should be capable of using, find all
# reactions that can be removed from said network without removing its ability
# to produce biomass given only those food sources
# network must already have a biomass reaction set as its objective
def find_prunables(network, foods, bm_rxn):
    # list of reactions that could be removed from this network
    could_drop = list()
    for rxn in network.reactions:
        # don't try to remove the biomass reaction
        if rxn.id == bm_rxn.id:
            continue
        # reaction removal happens in-place
        network.remove_reactions([rxn])
        # see if this network can produce biomass
        solution = network.optimize()
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        if solution.status != 'infeasible' and bm_rxn_flux > 10e-10:
            could_drop.append(rxn)
        # regardless of whether this reaction was removeable or not, add it 
        # back to the network so we can try another one
        network.add_reaction(rxn)
    return(could_drop)

# get command-line arguments
try:
    (monos, max_pol, ins, outs) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments: monomers, max polymer length, number of food sources, \
number of biomass precursors.')

# create the reference network and pick a biomass rxn
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
cobra_model.objective = bm_rxn

# pick some food sources and add input reactions to the network
foods = random.sample(cobra_model.metabolites, int(ins))
for food in foods:
    # add an input reaction for this metabolite
    in_rxn = cobra.Reaction(
        '->' + food.id,
        upper_bound = 100.0,
        lower_bound = 0.0
    )
    in_rxn.add_metabolites({food: 1.0})
    cobra_model.add_reaction(in_rxn)

# with randomly-chosen food metabolites, you'll frequently get no feasible
# solutions even on the complete network, especially for small networks, so
# keep choosing new food sources and biomass reactions until one combo
# actually works
could_remove = find_prunables(cobra_model, foods, bm_rxn)
i = 0
while could_remove == list():
    i += 1
    if i % 100 == 0:
        print(f'On the {i}th combination of foods and biomass')
    # remove existing food input reactions
    cobra_model.remove_reactions(cobra_model.boundary)
    foods = random.sample(cobra_model.metabolites, int(ins))
    for food in foods:
        # add an input reaction for this metabolite
        in_rxn = cobra.Reaction(
            '->' + food.id,
            upper_bound = 100.0,
            lower_bound = 0.0
        )
        in_rxn.add_metabolites({food: 1.0})
        cobra_model.add_reaction(in_rxn)
    cobra_model.remove_reactions([bm_rxn])
    bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
    cobra_model.objective = bm_rxn
    could_remove = find_prunables(cobra_model, foods, bm_rxn)

print(
    f'Reselected foods and biomass precursors {i} times before getting a ' +
    'feasible solution on the complete network.'
)
# now we can actually try to prune this network, since we know this combination
# of food sources and biomass precursors has at least one feasible solution on
# this network

# list of all possible pruned networks from this network
end_prunes = list()
# at first we just have one network to find possible prunes for: the complete
# network of all possible reactions
networks = [cobra_model]
print(f'{len(cobra_model.reactions)} reactions in full network')
while len(networks) != 0:
    print(f'Testing {len(networks)} networks')
    # this will hold all networks that have one fewer reaction than the current
    # network (or group of equally-sized networks) that still produce biomass
    new_networks = list()
    for network in networks:
        # find all possible reactions that could be pruned from this network
        # and add the networks missing those reactions to new_networks
        could_drop = find_prunables(network, foods, bm_rxn)
        if could_drop == list():
            # if no reactions could be removed from this network, add it to the
            # list of completely-pruned networks
            end_prunes.append(network)
            # if all the networks are like this, then new_networks will remain
            # an empty list and we'll exit the while loop
        else:
            for rxn in could_drop:
                # removing reactions happens in-place, so we need to make a new
                # network each time we want to prune
                pruned = network.copy()
                # annoyingly, the reaction object from "network" won't exist in
                # "pruned", but it will have the same id as the relevant one in
                # "pruned"
                rxn = pruned.reactions.get_by_id(rxn.id)
                pruned.remove_reactions([rxn])
                new_networks.append(pruned)
    # once you've gone through all the networks of this size, start over on all
    # of the networks with one fewer reaction until there are no valid smaller
    # networks
    networks = new_networks
print(f'There were {len(end_prunes)} ways to prune this network.')
