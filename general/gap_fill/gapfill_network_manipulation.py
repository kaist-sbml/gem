'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

from cobra import Model, Reaction, Metabolite
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

    #This is an empty model
    universal_model2 = Model("universal model wihout overlapping reactions")

    #Add reactions in universal network to template model
    for mnxr in balanced_unique_mnxr_list:
    #for mnxr in universal_model.reactions:
        mnxr_from_universal_model = universal_model.reactions.get_by_id(mnxr)

        #cobrapy doc: "When copying a reaction, it is necessary to deepcopy the components so the list references are not carried over"
        mnxr_from_universal_model = copy.deepcopy(mnxr_from_universal_model)
        #Integrated model having reactions from both target and universal models
        target_model.add_reaction(mnxr_from_universal_model)
        #A model having only reactions unique to universal model
        universal_model2.add_reaction(mnxr_from_universal_model)

    return target_model, universal_model2
