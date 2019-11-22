# multiple_env_prune.py
# makes one network, picks one biomass reaction, picks some groups of n input
# metabolites, prunes one network using one of those groups a lot of times
# to get several networks that work with that group of metabolites, and then
# filters those networks based on whether or not they still grow when supplied
# with all the other groups of metabolites

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
    (monos, max_pol, ins, in_groups, outs, reps) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments: monomers, max polymer length, number of food sources, \
number of groups of food sources, number of biomass precursors, number of \
times to prune the first network')

# create the reference network and pick a biomass reaction
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
print(f'Biomass reaction: {bm_rxn.id}')
cobra_model.objective = bm_rxn

# generate the specified number of food source groups
food_groups = list()
i = 0
while i < int(in_groups):
    # make sure we don't pick any biomass precursors but don't worry about 
    # picking metabolites that are already in another group
    in_group = random.sample([
        met for met in cobra_model.metabolites if met not in bm_rxn.metabolites
    ], int(ins))
    # make sure this group of inputs works on the full network
    for met in in_group:
        in_rxn = cobra.Reaction(
            '->' + met.id,
            upper_bound = 100.0,
            lower_bound = 0.0
        )
        in_rxn.add_metabolites({met: 1.0})
        cobra_model.add_reaction(in_rxn)
    solution = cobra_model.optimize()
    # sometimes flux through biomass is abyssmally small but solution is still
    # technically feasible, but we don't actually want those solutions
    bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
    if solution.status != 'infeasible' and bm_rxn_flux > 10e-10:
        # only proceed if this was a valid group of metabolites
        i = i + 1
        food_groups.append(in_group)
    # if this wasn't a workable choice of inputs, don't advance the loop
    # counter
    # regardless of whether or not this was a valid choice of metabolites,
    # remove the input reactions that were added so that more can be tested
    cobra_model.remove_reactions(cobra_model.boundary)

# start by pruning the network on one group of input metabolites a bunch of
# times, then worry about the other groups we made
for met in food_groups[0]:
    in_rxn = cobra.Reaction(
        '->' + met.id,
        upper_bound = 100.0,
        lower_bound = 0.0
    )
    in_rxn.add_metabolites({met: 1.0})
    cobra_model.add_reaction(in_rxn)

print('Pruning full network on first set of input metabolites')
# will hold bitstrings of all unique reactions and the count of times each
# one came up
pruned_dict = dict()
# will hold all the unique networks found by random_prune after reps runs
pruned_nets = list()
i = 0
while i < int(reps):
    i += 1
    if i % 100 == 0:
        print(f'Pruned {i} times')
    pruned_net = scn.random_prune(cobra_model, bm_rxn)
    # in order to know whether we've seen this model before, we can't just 
    # compare models, since no two models are ever 'equal', so we'll compare
    # reaction presence bitstrings. 
    # we also want to keep track of how many times we see each model, so we 
    # will make a dict with the bitstrings as keys
    # remove the input reactions from this network so the next step (where we
    # change the input reactions) actually works
    pruned_net.remove_reactions(pruned_net.boundary)
    bitstring = scn.make_bitstring(cobra_model, pruned_net)
    if bitstring not in pruned_dict.keys():
        # make sure all reaction lists are sorted so that all isomorphic
        # networks have the same reaction list
        pruned_dict[bitstring] = 1
        pruned_nets.append(pruned_net)
    else:
        # if we already found this network once, then increment the
        # appropriate counter by 1
        pruned_dict[bitstring] += 1

# now that we have (approximately, if reps was large) all of the networks that
# the first food group worked on, see which of those networks also work with
# the other food groups
# store this information as a dict with network bitstrings as keys and the list
# of food groups as values
usable_foods = dict()
for network in pruned_nets:
    bitstring = scn.make_bitstring(cobra_model, network)
    # we already determined that they all grow when given the first food group
    usable_foods[bitstring] = [met.id for met in food_groups[0]]
for group in food_groups[1:]:
    for network in pruned_nets:
        for met in group:
            in_rxn = cobra.Reaction(
                '->' + met.id,
                upper_bound = 100.0,
                lower_bound = 0.0
            )
            in_rxn.add_metabolites({met: 1.0})
            network.add_reaction(in_rxn)
        # now see if this network can produce biomass
        solution = network.optimize()
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        if solution.status != 'infeasible' and bm_rxn_flux > 10e-10:
            # if the network grew using this food group, add it to the
            # appropriate list in usable_foods
            bitstring = scn.make_bitstring(cobra_model, network)
            usable_foods[bitstring].append([met.id for met in group])

for network in usable_foods.keys():
    print(network)
    print(f'This network showed up {pruned_dict[network]} times')
    print('This network could produce biomass in these environments:' + 
        ','.join([
            '(' + ','.join(foods) + ')' for foods in usable_foods[network]
        ]
    )
