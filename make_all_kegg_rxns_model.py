# make_all_kegg_rxns_model.py
# using the list of all reactions annotated in KEGG, create a COBRApy model

import cobra
import re
import sys

# start by reading in the file from KEGG; it came from this URL:
# rest.kegg.jp/list/reaction
rxn_strings = list()
with open('data/kegg_reactions.txt', 'r') as f:
    rxn_strings = f.readlines()

# create an empty COBRA model
cobra_model = cobra.Model('full_kegg')

# now loop over each string, splitting them into the reaction ID, name (if 
# there is one), reactants and products, while assembling a COBRA reaction
# object, adding it to the above model when it's done
for rxn_string in rxn_strings:
    # reaction ID is separated from everything else by a tab
    (rxn_id, rest) = rxn_string.rstrip('\n').split('\t')

    # there may or may not be an enzyme/reaction name at the beginning of the
    # rest, separated by a semicolon. there may be multiple names also 
    # separated by semicolons, so we have to just take the last element of the
    # list returned by split
    maybe = rest.split('; ')

    # initialize rxn_name and reaction as empty strings
    rxn_name = str()
    reaction = str()
    if len(maybe) == 1:
        # there was no enzyme/reaction name
        rxn_name = str()
        reaction = maybe[0]
    elif len(maybe) != 1:
        # there was at least one enzyme/reaction name; the last element is the
        # actual reaction, and everything else is name(s)
        rxn_name = '; '.join(maybe[:-1])
        reaction = maybe[-1]
    
    # create the COBRA reaction object using the id and name we just got
    cobra_rxn = cobra.Reaction(
        rxn_id, name = rxn_name,
        # make all reactions reversible; this seems to be the norm
        lower_bound = -100, upper_bound = 100
    )

    # split reaction into reactants and products
    (reactant_string, product_string) = reaction.split(' <=> ')
    
    # separate individual reactants and products; split on ' + ' and not just
    # '+' because some metabolites have +s in them, like H+ or NAD+
    reactants = reactant_string.split(' + ')
    products = product_string.split(' + ')

    # make dictionary of metabolites and stoichiometric coefficients
    met_dict = dict()

    # start with the reactants
    for reactant in reactants:
        # see if there's a number at the beginning of this metabolite's name
        # separated from the rest of the name with a space; this is the
        # stoichiometric coefficient
        stoich = 0
        reac_name = str()
        if re.match('^\d+ ', reactant) is not None:
            # this string contains both the name and stoichiometric coefficient
            # separate the two things and make the coefficient an float
            stoich = float(reactant.split(' ')[0])
            reac_name = re.sub('^\d+ ', '', reactant)
            # metabolite names can't have whitespace characters
            reac_name = re.sub(' ', '_', reac_name)
        elif re.match('^\d+ ', reactant) is None:
            # this string only has the metabolite's name, so the stoichiometric
            # coefficient must be 1
            stoich = 1.0
            # metabolite names can't have whitespace characters
            reac_name = re.sub(' ', '_', reactant)
        # since these are reactants, make the stoichiometric coefficients
        # negative
        # the compartment argument is necessary if we want to optimize this
        # network down the road (and we do)
        met_dict[cobra.Metabolite(reac_name, compartment = 'c')] = -stoich

    # do the same thing for the products
    for product in products:
        # deal with stoichiometric coefficient
        stoich = 0
        prod_name = str()
        if re.match('^\d+ ', product) is not None:
            # separate stoich and reac_name
            stoich = float(product.split(' ')[0])
            prod_name = re.sub('^\d+ ', '', product)
            # deal with whitespace
            prod_name = re.sub(' ', '_', prod_name)
        elif re.match('^\d+ ', product) is None:
            stoich = 1.0
            # deal with whitespace
            prod_name = re.sub(' ', '_', product)
        # don't need to negate coefficients this time
        met_dict[cobra.Metabolite(prod_name, compartment = 'c')] = stoich

    # add this dictionary to the COBRA reaction object
    cobra_rxn.add_metabolites(met_dict)

    # now this COBRA reaction object is complete and can be added to the model
    cobra_model.add_reactions([cobra_rxn])

# print some numbers
print(f'Number of reactions: {len(cobra_model.reactions)}')
print(f'Number of metabolites: {len(cobra_model.metabolites)}')

# write model
cobra.io.save_json_model(cobra_model, 'data/full_kegg_cobra_model.json')
