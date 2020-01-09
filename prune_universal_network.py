# prune_universal_network.py
# prunes the universal BIGG model using the E. coli environment and biomass
# reaction from iJO1336

import cobra
import string_chem_net as scn
import sys

# need to load both the universal model and an E. coli model to get the biomass
# reaction
print('Loading COBRA models')
e_coli_model = cobra.io.load_json_model('data/iJO1366.json')
u_model = cobra.io.load_json_model('data/bigg_universal_model.json')

# remove all exchange reactions from the unviersal model
print('Removing non-E. coli exchange reactions from universal model')
u_model.remove_reactions(u_model.boundary)

# add each remaining reaction in the universal model to the E. coli model
# one at a time
print(
    'Testing remaining reactions to see if adding them to E. coli model ' +
    'breaks the model.'
)
i = 0
bad = 0
for rxn in u_model.reactions:
    i += 1
    if i % 100 == 0:
        print(f'On reaction {i}; {bad} reactions were not added.')
    # skip reactions already in the E. coli model
    if rxn not in e_coli_model.reactions:
        e_coli_model.add_reactions([rxn])
    else:
        next
    try:
        e_coli_model.medium
    except AttributeError:
        bad += 1
        e_coli_model.remove_reactions([rxn])

print(f'Added {i - bad} of {i} possible reactions.')

# optimize the expanded E. coli model
solution = e_coli_model.optimize()
print(e_coli_model.summary())
