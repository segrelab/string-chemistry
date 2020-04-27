# random_prune_env_tests.py
# does random pruning on a single universal scale network and then tests
# whether or not each of the unique randomly-pruned networks can grow in a 
# specified number of environments
# the environments are separated into environments you expect growth in and
# environments you don't expect growth in to simulate what you might do with
# real metabolic networks if you had growth data for a real organism

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
    (monos, max_pol, ins, yes_groups, no_groups, outs, reps) = sys.argv[1:]
except ValueError:
    sys.exit('Arguments:\nmonomers\nmax polymer length\nnumber of food ' +
        'sources in each environment\nnumber of environments to grow in\n' +
        'number of environments to not grow in\nnumber of biomass ' + 
        'precursors\nnumber of times to prune the first network.'
    )

# create the reference network and pick a biomass reaction
SCN = scn.CreateNetwork(monos, int(max_pol))
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
bm_rxn = scn.choose_bm_mets(int(outs), cobra_model)
print(f'Biomass reaction: {bm_rxn.id}')
cobra_model.objective = bm_rxn

print('Generating random environments the networks should grow in')
# generate the specified number of environments that the network should grow in
yes_envs = list()
i = 0
while i < int(yes_groups):
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
        yes_envs.append(in_group)
    # if this wasn't a workable choice of inputs, don't advance the loop
    # counter
    # regardless of whether or not this was a valid choice of metabolites,
    # remove the input reactions that were added so that more can be tested
    cobra_model.remove_reactions(cobra_model.boundary)

print('Generating random environments the networks should not grow in')
# generate the specified number of environments that the network should not 
# grow in
no_envs = list()
for group in range(int(no_groups)):
    # make sure we don't pick any biomass precursors but don't worry about 
    # picking metabolites that are already in another group
    in_group = random.sample([
        met for met in cobra_model.metabolites if met not in bm_rxn.metabolites
    ], int(ins))
    # this time, we don't need to check if the full network can grow given
    # these metabolites as input, since we are looking for networks that won't
    # grow with these
    no_envs.append(in_group)

# start by pruning the network on one group of input metabolites a bunch of
# times, then worry about the other groups we made
for met in yes_envs[0]:
    in_rxn = cobra.Reaction(
        '->' + met.id,
        upper_bound = 100.0,
        lower_bound = 0.0
    )
    in_rxn.add_metabolites({met: 1.0})
    cobra_model.add_reaction(in_rxn)

print(
    f'Pruning full {len(cobra_model.reactions)}-reaction network on first ' +
    'environment'
)
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

print('Seeing which networks grow in the environments they should grow in.')
# now that we have (approximately, if reps was large) all of the networks that
# the first food group worked on, see which of those networks also work with
# the other food groups it should be able to grow on
# store this information as a dict with network bitstrings as keys and the list
# of food groups as values
usable_foods = dict()
for network in pruned_nets:
    bitstring = scn.make_bitstring(cobra_model, network)
    # we already determined that they all grow when given the first food group
    usable_foods[bitstring] = [[met.id for met in yes_envs[0]]]
# start at the second group
for group in yes_envs[1:]:
    for network in pruned_nets:
        # add input reactions for all metabolites in this group
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
            # drop the boundary reactions before making the bitstirng
            # since those will vary
            network.remove_reactions(network.boundary)
            bitstring = scn.make_bitstring(cobra_model, network)
            usable_foods[bitstring].append([met.id for met in group])
        else:
            # even if this group wasn't viable, remove these input reactions
            # so the next group can be tested separately
            network.remove_reactions(network.boundary)

print('Seeing which networks don\'t grow in the environments they shouldn\'t')
# now see if each network doesn't grow on the food sources we don't want them
# to grow on
# once again store this as a dict with bitstring keys and lists of metabolites
# as values
unusable_foods = dict()
for network in pruned_nets:
    bitstring = scn.make_bitstring(cobra_model, network)
    # this time we don't know if any of the pruned networks won't grow on these
    # so we need to initialize with empty lists
    unusable_foods[bitstring] = list()
for group in no_envs:
    for network in pruned_nets:
        # add input reactions for all metabolites in this group
        for met in group:
            in_rxn = cobra.Reaction(
                '->' + met.id, upper_bound = 100.0, lower_bound = 0.0
            )
            in_rxn.add_metabolites({met: 1.0})
            network.add_reaction(in_rxn)
        # now see if this network can produe biomass
        solution = network.optimize()
        # sometimes the solution is feasible but the flux through the biomass
        # reaction is negligible, so check for that as well
        bm_rxn_flux = solution.fluxes.get(key = bm_rxn.id)
        if solution.status == 'infeasible' or bm_rxn_flux < 10e-10:
            # if this network couldn't produce biomass, add it to the
            # appropriate list
            # drop the boundary reactions before making the bitstring since 
            # those will be variable
            network.remove_reactions(network.boundary)
            bitstring = scn.make_bitstring(cobra_model, network)
            unusable_foods[bitstring].append([met.id for met in group])
        else:
            # if we didn't add this to the list, drop the boundary reactions
            # so that they're not around for the next iteration
            network.remove_reactions(network.boundary)

# print a bunch of info but also write it out to a tsv
with open(
        f'../data/{monos}_{max_pol}_{yes_groups}yes_{no_groups}no_{ins}_' +
        f'{outs}outs.tsv', 'w'
    ) as out:
    out.write('bitstring\trxn_count\toccurrences\tyes_count\tyes_envs\t' +
        'no_count\tno_envs\tbiomass\n')
    for network in usable_foods.keys():
        rxn_count = count_bitstring(network)
        print(
            f'A {rxn_count}-reaction network showed up ' +
            f'{pruned_dict[network]} times and could produce biomass in ' +
            f'these {len(usable_foods[network])} environments:'
        )
        for foods in usable_foods[network]:
            print(','.join(foods))
        print(
            'It could not produce biomass in these ' +
            f'{len(unusable_foods[network])} environments:'
        )
        for foods in unusable_foods[network]:
            print(','.join(foods))
        out_row = '\t'.join([
            network, str(rxn_count), str(pruned_dict[network]),
            str(len(usable_foods[network])), 
            ';'.join([','.join(foods) for foods in usable_foods[network]]),
            str(len(unusable_foods[network])),
            ';'.join([','.join(foods) for foods in unusable_foods[network]]),
            bm_rxn.id
        ])
        out.write(out_row + '\n')
