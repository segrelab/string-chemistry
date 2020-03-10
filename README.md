# string-chemistry

Tools for creating string chemistry networks and doing flux-balance analysis on them.

## Scripts That Make Networks

### string\_chem\_net.py

THe key script; contains a class for generating arbitrary string chemistry networks and several functions for maniuplating them and turning them into COBRApy models for FBA.

#### Classes:

##### CreateNetwork

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

#### Functions:

##### choose\_bm\_mets

Randomly choose n metabolites (without replacement) from a given network and make an exchange reaction that consumes all of them (i.e. a biomass reaction)
NOTE: Does not set the biomass reaction as the objective for the model; if you want to do FBA with this biomass reaction as the objective, you need to do that separately

Arguments:
- `n`: number of biomass precursors to choose
- `model`: a COBRApy model to add this biomass reaction to

Returns:

- `bm_rxn`: the COBRApy reaction object corresponding to the newly-created biomass reaction, so you can set it as the objective for the model or just have it around for whatever purpose

##### choose\_inputs

Randomly choose n metabolites (without replacement) from a given network and make exchange reactions that produce them

Arguments:

- `n`: number of input metabolites to choose
- `model`: a COBRApy model to add these exchange reactions to
- `bm_rxn`: the biomass reaction of the model, if one exists (so you don't get an exchange reaction for a metabolite that's in the biomass reaction). Defaults to an empty reaction (in case your model doesn't already have a biomass reaction defined)

Returns:
None; COBRApy models are edited in-place

##### make\_cobra\_model

Makes a COBRApy model so you can do FBA on a string chemistry network

Arguments:

- `met_list`: list of metabolites in the network
- `rxn_list`: list of reactions in the network

Returns:

A [COBRApy mode](lhttps://cobrapy.readthedocs.io/en/latest/index.html)

##### make\_edgelist

Makes an edgelist to facilitate visualization of a network using, e.g. Cytoscape or graphviz

Arguments:

- `rxn_list`: list of reactions in a network
- `rxns_as_nodes`: should reactions be nodes connected to their associated metabolites? Default is True; if False, all nodes are metabolites, which are connected iff they both participate in the same reaction

Returns:

A list of lists where each sublist is two elements: the names of the pair of metabolites or the metabolite-reaction pair connected by that edge

##### remove\_random\_rxns

Randomly removes reactions from a network given a removal probability.

Arguments:

- `more\_rxns`: reaction list for the network to be edited
- `S`: stoichiometric matrix for the network to be edited
- `prob`: probability that a reaction should be removed

Returns:

- `less_rxns`: reaction list for smaller network
- `smaller_S`: stoichiometric for smaller network

At some point, I'll change this to make including S optional, since CreateNetwork no longer creates an S by default.

##### reverse\_rxns

Makes n randomly-chosen (without replacement) reactions in a COBRApy model reversible (by setting their lower bounds to -100)
Note that CreateNetwork makes all reactions irreversible by default, hence the existence of this function

Arguments:

- `model`: a COBRApy model
- `n`: number of reactions to make reversible

## Plotting Scripts

## Other Scripts

### counting.py

Arguments:

- n: number of unique metabolites
- l: maximum polymer length

Prints number of metabolites and reactions in a string chemistry network with those parameters; does not actually generate that network.
