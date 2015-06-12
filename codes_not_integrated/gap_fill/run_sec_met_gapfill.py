'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file
from gapfill_network_manipulation import (
    get_mnxr_bigg_in_target_model,
    get_mnxr_unique_to_universal_model,
    get_balanced_rxns_from_mnxr,
    get_manipulated_target_universal_models
)
from gapfill_core import gapfilling_precursor
import pickle

target_model = create_cobra_model_from_sbml_file('./E_coli_iAF1260_rewrite.xml')
#cobra_model = create_cobra_model_from_sbml_file('./sma_target_model_sco_complete.xml')

#ivcoa_c, 3-Methylbutanoyl-CoA: one of the precursors for non-ribosomal peptide
if 'ivcoa_c' in target_model.metabolites:
    print "\n"
    print "ivcoa_c (ivcoa_c) is already included in template_model"
else:
    print "\n"
    print "'Metabolite (ivcoa_c) is not included in template_model"

#Load universal network in sbml format
#universal_model = create_cobra_model_from_sbml_file('./universal_network_fiexed_bigg_mnxref.xml')

# fixing the definition of metabolite id (e.g. ala__L_c --> ala_DASH_L_c)
#universal_model = fix_special_characters_compoundid(universal_model)

#pickle.dump(universal_model, open('./universal_model.p','wb'))
universal_model = pickle.load(open("universal_model.p","rb"))

#From gapfill_network_manipulation.py
bigg_mnxr_dict = pickle.load(open("bigg_mnxr_dict.p","rb"))

mnxr_bigg_target_model_dict = get_mnxr_bigg_in_target_model(target_model, bigg_mnxr_dict)

mnxr_unique_to_universal_model_list = get_mnxr_unique_to_universal_model(mnxr_bigg_target_model_dict, universal_model)

mnxr_rxn_all_dict = pickle.load(open("mnxr_rxn_all_dict.p","rb"))

balanced_unique_mnxr_list = get_balanced_rxns_from_mnxr(mnxr_unique_to_universal_model_list, mnxr_rxn_all_dict)

target_model, universal_model2 = get_manipulated_target_universal_models(balanced_unique_mnxr_list, target_model, universal_model)


############ Generation of metabolic reaction of objective function for desired products
#TO REMOVE
metabolite_id = 'ivcoa'
comp_met_id = metabolite_id+'_c'

if comp_met_id in target_model.metabolites:
    demand_reaction = Reaction("Gapfill_DEMAND_RXN_%s" % metabolite_id)
    demand_reaction.name = "Gapfill_DEMAND_RXN_%s" % metabolite_id

    demand_reaction.lower_bound = 0
    demand_reaction.upper_bound = 999999
    demand_reaction.reversibility = 1

    comp_met_id_c = target_model.metabolites.get_by_id(comp_met_id)

    demand_reaction.add_metabolites({comp_met_id_c: -1})

    comp_met_id_e = Metabolite(metabolite_id+'_e', name = '', compartment='e')
    demand_reaction.add_metabolites({comp_met_id_e: 1})

    target_model.add_reaction(demand_reaction)

    #Creating an exchange reaction
    #Creating reaction ID
    export_reaction = Reaction("Gapfill_Export_RXN_%s" % metabolite_id)

    #Setting bounds
    export_reaction.reversibility = 1 # 1: reversible 0: irreversible
    export_reaction.lower_bound = 0
    export_reaction.upper_bound = 999999

    #Adding a substrate metabolite
    #print target_model.metabolites.get_by_id(str(product_c))
    export_reaction.add_metabolites({comp_met_id_e:-1})

    #Adding the new reaction to the model
    target_model.add_reaction(export_reaction)


#Change objective function from biomass to desired precursor
for reaction in target_model.reactions:
    if reaction.objective_coefficient == 1:
        print "the objective function in original model is %s" % (reaction.id)
        print "the related reaction is %s" % (reaction.reaction)

for reaction in target_model.reactions:
    reaction.objective_coefficient = 0

demand_reaction.objective_coefficient = 1


#Run gap-filling algorithm based on MILP in gurobipy
obj = gapfilling_precursor()

#Load merged model
obj.load_cobra_model(target_model)

obj.fill_gap(demand_reaction.id, target_model, universal_model2)
