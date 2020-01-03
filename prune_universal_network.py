# prune_universal_network.py
# adds a glucose import reaction and the E. coli biomnass reaction to the
# universal biochemical reaction network from KEGG and prunes it

import cobra
import string_chem_net as scn

# need to load both the universal model and an E. coli model to get the biomass
# reaction
e_coli_model = cobra.io.read_sbml_model('data/ecoli_core_model.xml')
u_model = cobra.io.load_json_model('data/full_kegg_cobra_model.json')

# make a glucose-supplying reaction
glc_in_rxn = cobra.Reaction()
glucose_met = u_model.metabolites.get_by_id('D-Glucose')
glc_in_rxn.add_metabolites({glucose_met : 1.0})

# get biomass reaction for E. coli
bm_rxn = e_coli_model.reactions.get_by_id('Biomass_Ecoli_core_w_GAM')

# add these to the universal model and set the objective
u_model.add_reactions([glc_in_rxn, bm_rxn])
u_model.objective = bm_rxn

# see if the network is viable
solution = u_model.optimize()
print(solution.summary())
