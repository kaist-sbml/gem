'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
import copy

#Similar to replace_special_characters_compoundid(biggid) in independentModule.py
def fix_special_characters_compoundid(model):

    for each_metabolite in model.metabolites:
        
        met_id = each_metabolite.id
        
        met_id = met_id.replace('__', '_DASH_')
        met_id = met_id.replace('/', '_FSLASH_')
        met_id = met_id.replace("\\", '_BSLASH_')
        met_id = met_id.replace('(', '_LPAREN_')
        met_id = met_id.replace('[', '_LSQBKT_')
        met_id = met_id.replace(']', '_RSQBKT_')
        met_id = met_id.replace(')', '_RPAREN_')
        met_id = met_id.replace(',', '_COMMA_')
        met_id = met_id.replace('.', '_PERIOD_')
        met_id = met_id.replace("'", '_APOS_')
        met_id = met_id.replace('&', '&amp;')
        met_id = met_id.replace('<', '&lt;')
        met_id = met_id.replace('>', '&gt;')
        met_id = met_id.replace('"', '&quot;')
        
        conv_met_id = met_id
        each_metabolite.id = conv_met_id

    return model


def get_mnxr_bigg_in_target_model(target_model, bigg_mnxr_dict):

    mnxr_bigg_target_model_dict = {}

    for each_reaction in bigg_mnxr_dict.keys():
        if each_reaction in target_model.reactions:
            biggid = each_reaction
            mnxr = bigg_mnxr_dict[each_reaction]
            mnxr_bigg_target_model_dict[mnxr] = biggid

    return mnxr_bigg_target_model_dict


def get_mnxr_unique_to_universal_model(mnxr_bigg_target_model_dict, universal_model):

    mnxr_unique_to_universal_model_list = []

    for mnxr in universal_model.reactions:
        if mnxr not in mnxr_bigg_target_model_dict.keys():
            mnxr_unique_to_universal_model_list.append(str(mnxr))

    return mnxr_unique_to_universal_model_list


def get_balanced_rxns_from_mnxr(mnxr_unique_to_universal_model_list, mnxr_rxn_all_dict):

    balanced_unique_mnxr_list = []

    for mnxr in mnxr_unique_to_universal_model_list:
        if mnxr in mnxr_rxn_all_dict.keys():
            if mnxr_rxn_all_dict[mnxr][2] == 'true':
                balanced_unique_mnxr_list.append(mnxr)

    return balanced_unique_mnxr_list


def get_manipulated_target_universal_models(balanced_unique_mnxr_list, target_model, universal_model):

    target_model2 = copy.deepcopy(target_model)

    #This is an empty model
    universal_model2 = Model("universal model wihout overlapping reactions")

    #Add reactions in universal network to template model
    for mnxr in balanced_unique_mnxr_list:
    #for mnxr in universal_model.reactions:
        mnxr_from_universal_model = universal_model.reactions.get_by_id(mnxr)

        #cobrapy doc: "When copying a reaction, it is necessary to deepcopy the components so the list references are not carried over"
        mnxr_from_universal_model = copy.deepcopy(mnxr_from_universal_model)
        #Integrated model having reactions from both target and universal models
        target_model2.add_reaction(mnxr_from_universal_model)
        #A model having only reactions unique to universal model
        universal_model2.add_reaction(mnxr_from_universal_model)

    return target_model2, universal_model2


def add_transport_exchange_rxn_nonprod_monomer(target_model, nonprod_monomer):

    target_model_temp = copy.deepcopy(target_model)

    #Creat a transport reaction
    #Creat reaction ID
    rxn = Reaction("Transport_"+nonprod_monomer)
    rxn.name = "Transport_"+nonprod_monomer

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible

    #Add a substrate metabolite
    rxn.add_metabolites({target_model_temp.metabolites.get_by_id(nonprod_monomer+'_c'):-1})

    #Add product metabolite(s)
    product_e = Metabolite(nonprod_monomer+"_e", name='', compartment='e')
    rxn.add_metabolites({product_e:1})

    #Add the new reaction to the model
    target_model_temp.add_reaction(rxn)

    #Creat an exchange reaction
    #Creat reaction ID
    rxn = Reaction("Ex_"+nonprod_monomer)
    rxn.name = "Ex_"+nonprod_monomer

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible 0: irreversible

    #Add a substrate metabolite
    rxn.add_metabolites({target_model_temp.metabolites.get_by_id(str(product_e)):-1})

    #Add the new reaction to the model
    target_model_temp.add_reaction(rxn)

    return target_model_temp


def check_producibility_nonprod_monomer(target_model_temp, nonprod_monomer):
    #Change objective function from biomass to desired precursor
    #This is for Sco as a template model
    target_model_temp.reactions.get_by_id('Biomass_SCO').objective_coefficient = 0
    target_model_temp.reactions.get_by_id("Ex_"+nonprod_monomer).objective_coefficient = 1

    #Model reloading and overwrtting are necessary for model stability
    #Without these, model does not produce an accurate prediction
    write_cobra_model_to_sbml_file(target_model_temp, "target_model_temp_%s.xml" %nonprod_monomer)
    target_model_temp = create_cobra_model_from_sbml_file("target_model_temp_%s.xml" %nonprod_monomer)
    target_model_temp.optimize()

    print target_model_temp.reactions.get_by_id("Ex_"+nonprod_monomer)
    print target_model_temp.reactions.get_by_id("Ex_"+nonprod_monomer).reaction
    print "Flux:", target_model_temp.solution.f, "\n"

    #fp1 = open("%s_fba.txt" %nonprod_monomer,"w")
    #for the_reaction, the_value in target_model_temp.solution.x_dict.items():
    #    fp1.write(str(the_reaction)+"\t"+str(the_value)+"\n")
    #fp1.close()

    if target_model_temp.solution.f >= 0.001:
        return None
    else:
        return target_model_temp

def get_unique_nonprod_monomers_list(nonprod_sec_met):
    unique_nonprod_monomers_list = []

    for nonprod_monomers_list in nonprod_sec_met.keys():
        for nonprod_monomer in nonprod_sec_met[nonprod_monomers_list]:
            if nonprod_monomer not in unique_nonprod_monomers_list:
                unique_nonprod_monomers_list.append(nonprod_monomer)

    print unique_nonprod_monomers_list
    return unique_nonprod_monomers_list

