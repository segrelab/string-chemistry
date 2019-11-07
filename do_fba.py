# do_fba.py
import string_chem_net as scn
import time

SCN = scn.CreateNetwork('abc', 5)
cobra_model = scn.make_cobra_model(SCN.met_list, SCN.rxn_list)
scn.reverse_rxns(cobra_model, len(cobra_model.reactions))
scn.choose_inputs(cobra_model, 5)
print([rxn.id for rxn in cobra_model.boundary])
bm_rxn = scn.choose_bm_mets(cobra_model, 5)
cobra_model.objective = bm_rxn
solution = cobra_model.optimize()
print(cobra_model.summary())
