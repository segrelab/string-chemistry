# ARtificial CHemistry NEtwork Toolbox (ARCHNET)

Tools for creating string chemistry networks and doing flux-balance analysis on them. All scripts have a summary of what they do in the first few lines, and I've given more detailed descriptions of a few key ones below. The package itself is in string\_chem\_net.py, and everything else just uses the package in various ways. Nearly every script expects there to be another directory called `data` in the same directory as `scripts`, and that you're calling the script from that directory (rather than `scripts`). Every script should print the command-line arguments it expects if you run it without any arguments.

## The `ARCHNET` package
Contains a class for generating arbitrary string chemistry networks and several functions for maniuplating them and turning them into [COBRApy models](https://cobrapy.readthedocs.io/en/latest/index.html) for FBA.

### Quick Example Use:

(from same directory as string\_chem\_net.py, i.e. `scripts`)
```
import string_chem_net as scn

monos = 'ab' # characters to use as monomers
max_len = 5  # maximum length of each string chemical
ins = 2      # number of food sources / environmental nutrients
outs = 5     # number of biomass precursors

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

# now we can be confident that pruning will work
pruned_model = scn.min_flux_prune(model, bm_rxn)
```

If you want to work with the stoichiometric matrix directly, you can pass `make_stoich_mat = True` to CreateNetwork; it's set to false by default because it can take a while to create for large networks.

### Classes:

- `CreateNetwork`

    Generates all possible metabolites given the constraints and generates all possible bimolecular reactions involving only those metabolites.

    Arguments:

    - `monos`: number of unique metabolites
    - `max_len`: maximum polymer length
    - `no_mirrors`: should mirror-image metabolites be removed? Default is False
    - `make_stoich`: should a stoichiometric matrix for the model be generated (as a numpy array)? Default is False

    Returns:

    An object with the following attributes:

    - `met_list`: list of all possible string metabolites given `monos` and `max_len`
    - `met_set`: set of all possible string metabolites given `monos` and `max_len`
    - `rxn_list`: list of all possible bimolecular reactions involving only those metabolites in `met_list`
    - `S`: stoichiometric matrix corresponding to the network of reactions in `rxn_list` (only defined if `make_stoich` is True)

### Functions:

- `choose_bm_mets`

    Randomly choose n metabolites (without replacement) from a given network and make an exchange reaction that consumes all of them (i.e. a biomass reaction)
    NOTE: Does not set the biomass reaction as the objective for the model; if you want to do FBA with this biomass reaction as the objective, you need to do that separately

    Arguments:
    - `n`: number of biomass precursors to choose
    - `model`: a COBRApy model to add this biomass reaction to

    Returns:

    - `bm_rxn`: the COBRApy reaction object corresponding to the newly-created biomass reaction, so you can set it as the objective for the model or just have it around for whatever purpose

- `choose_inputs`

    Randomly choose n metabolites (without replacement) from a given network and make exchange reactions that produce them

    Arguments:

    - `n`: number of input metabolites to choose
    - `model`: a COBRApy model to add these exchange reactions to
    - `bm_rxn`: the biomass reaction of the model, if one exists (so you don't get an exchange reaction for a metabolite that's in the biomass reaction). Defaults to an empty reaction (in case your model doesn't already have a biomass reaction defined)

    Returns:
    None; COBRApy models are edited in-place

- `make_bitstring`

    Given two COBRApy models where one is a subset of the other (e.g. the input and output from one of the pruning functions), makes a string of 1s and 0s representing which reactions in the larger model are present in the smaller one. Reactions from a particular COBRApy model are always output in the same order, so you can use this on multiple models pruned from the same initial model and get multiple comparable bitstrings (many of the plotting scripts depend on this)

    Arguments:

    - `full_model`: the larger COBRApy model (must contain all of the reactions in `pruned_model`)
    - `pruned_model`: the smaller (probably pruned) COBRApy model

    Returns:

    A binary vector indicating which reactions in `full_model` are present in `pruned_model`

- `make_cobra_model`

    Makes a COBRApy model so you can do FBA on a string chemistry network

    Arguments:

    - `met_list`: list of metabolites in the network
    - `rxn_list`: list of reactions in the network
    - `allow_export`: True/False (default is False); if True, all metabolites in network can be exported

    Returns:

    A [COBRApy model](https://cobrapy.readthedocs.io/en/latest/index.html)

- `make_edgelist`

    Makes an edgelist to facilitate visualization of a network using, e.g. Cytoscape or graphviz

    Arguments:

    - `rxn_list`: list of reactions in a network
    - `rxns_as_nodes`: should reactions be nodes connected to their associated metabolites? Default is True; if False, all nodes are metabolites, which are connected iff they both participate in the same reaction

    Returns:

    A list of lists where each sublist is two elements: the names of the pair of metabolites or the metabolite-reaction pair connected by that edge

- `min_flux_prune`

    Removes reactions from a COBRApy model by:
    1. Doing FBA (the model has to have a reaction set as the objective and needs some input reactions)
    2. Removing all reactions that have no flux through them in the FBA solution
    3. Removing the reaction with the smallest flux (of those that remain)
    4. Doing FBA again
    5. If there was a feasible solution, go back to step 2. If there was not a feasible solution, replace the reaction that was just removed and return the model with however many reactions it has left

    Note that this function will always give the same answer if given the same inputs, while the other pruning function will not

    Arguments:
    
    - `cobra_model`: the COBRApy model to be pruned

    Returns:

    A COBRApy model; it makes a copy of the input model so that the input model isn't modified

- `random_prune`

    Removes reactions from a COBRApy model by:
    1. Doing FBA (the model has to have a reaction set as the objective and needs some input reactions)
    2. Removing all reactions that have no flux through them in the FBA solution
    3. Randomly removing one of the remaining reactions that is not the biomass reaction
    4. Doing FBA again
    5. If there was a feasible solution, go back to step 2. If there was not a feasible solution, replace the reaction that was just removed and return the model with however many reactions it has left

    Note that this function can return a wide variety of networks as output given the same input, while the other pruning function will always return the same result given the same input

    Arguments:
    
    - `cobra_model`: the COBRApy model to be pruned
    - `bm_rxn`: the biomass reaction of the model (if the biomass reaction got randomly removed the model would have no solution)

    Returns:

    A COBRApy model; it makes a copy of the input model so that the input model isn't modified

- `remove_random_rxns`

    Randomly removes reactions from a network given a removal probability.

    Arguments:

    - `more_rxns`: reaction list for the network to be edited
    - `S`: stoichiometric matrix for the network to be edited
    - `prob`: probability that a reaction should be removed

    Returns:

    - `less_rxns`: reaction list for smaller network
    - `smaller_S`: stoichiometric for smaller network

    At some point, I'll change this to make including S optional, since CreateNetwork no longer creates an S by default.

- `reverse_rxns`

    Makes n randomly-chosen (without replacement) reactions in a COBRApy model reversible (by setting their lower bounds to -100)
    Note that CreateNetwork makes all reactions irreversible by default, hence the existence of this function

    Arguments:

    - `model`: a COBRApy model
    - `n`: number of reactions to make reversible


## Scripts That Use `ARCHNET`

### `exhaustive_prune.py`

Finds every possible way to prune a given network, i.e. every subnetwork that can still grow (i.e. have flux through the biomass reaction). Runs forever / crashes due to memory issues on all but the smallest of networks; was mostly created to see just how fast the number of possible pruned networks exploded beyond the point of computational feasibility.

Arguments:

- `monos`: number of different kinds of monomers to include in the network
- `max_pol`: maximum length of a polymer in the network
- `ins`: number of input reactions to add to the network before pruning (i.e. "food" sources; randomly chosen)
- `outs`: number of biomass precursors (randomly chosen)

### `pruning.py`

Makes a string chemistry network, prunes it once using the minimum-flux pruner, prunes it n times with the random pruner, recording which reactions were left after each round of pruning and how many times each network was returned as a result.

Arguments:

- `monos`: number of different kinds of monomers to include in the network
- `max_pol`: maximum length of a polymer in the network
- `ins`: number of input reactions to add to the network before pruning (i.e. "food" sources; randomly chosen)
- `outs`: number of biomass precursors (randomly choosen)

Returns:

CSV file with the following columns:

- Reaction inclusion bitstrings for all pruned networks (using the full un-pruned network as a reference for all of them so all the bitstrings are directly comparable)
- Number of times the network was returned by the random pruner (1 for the minimum-flux pruned network)
- Number of reactions in the network

## Other Scripts

- `bm_edit_dist.py`
    
    Looks at the output of multiple_env_prune.py, which has a set of food source metabolites and biomass precursors on each line, and calculates the average edit distance between every possible pair of either the food source metabolite lists or the biomass precursor lists. E.g.: line 1 has ab and bb as food source metabolites, while line 2 has aa and bc as food source metabolites. The average edit distance between the two sets of metabolites is (1 [ab <-> aa] + 2 [ab <-> bc] + 2 [bb <-> aa] + 1 [bb <-> bc]) / 4 = 1.5. This could concievably be helpful information when interpreting the results of multiple_env_prune.py, since you might expect networks with more similar input and output metabolites to get pruned into more similar networks, and this would be a way to quantify that.
    I don't think I ever finished debugging this script so I doubt it works yet. 

    Arguments:
        
    None; at the moment it is hard-coded to look at the default-named output file from multiple_env_prune.py. Eventually, it would be good if multiple_env_prune named its output files differently for each run and this script took a command-line argument specifying a target file, but that is not likely to happen soon.

- `network_sizes.py`

    Prints number of metabolites and reactions in a string chemistry network with those parameters; does not actually generate that network.

    Arguments:

    - `n`: number of unique metabolites
    - `l`: maximum polymer length

The other scripts mostly do a lot of pruning by pruning the same network in a wide variety of ways and saving varying pieces of information about the pruned networks. The rest have some examples of how one might visualize various things about string chemistry networks, both pruned and not.
