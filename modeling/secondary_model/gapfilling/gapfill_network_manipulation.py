'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

import copy
import logging
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from cobra.manipulation.delete import prune_unused_metabolites
from gapfill_core import gapfilling_precursor
from ..primary_model.augPhase_utils import(
    get_exrxnid_flux,
    check_exrxn_flux_direction
)


def get_mnxr_bigg_in_target_model(target_model, options):

    mnxr_bigg_target_model_dict = {}

    for each_reaction in options.bigg_mnxr_dict.keys():
        if each_reaction in target_model.reactions:
            biggid = each_reaction
            mnxr = options.bigg_mnxr_dict[each_reaction]
            mnxr_bigg_target_model_dict[mnxr] = biggid

    options.mnxr_bigg_target_model_dict = mnxr_bigg_target_model_dict


def get_mnxr_unique_to_universal_model(universal_model, options):

    mnxr_unique_to_universal_model_list = []

    for mnxr in universal_model.reactions:
        if mnxr not in options.mnxr_bigg_target_model_dict.keys():
            mnxr_unique_to_universal_model_list.append(str(mnxr))

    options.mnxr_unique_to_universal_model_list = mnxr_unique_to_universal_model_list


def integrate_target_universal_models(target_model, universal_model, options):
    target_model2 = copy.deepcopy(target_model)

    #Add reactions in universal network to template model
    for mnxr in options.mnxr_unique_to_universal_model_list:
        mnxr_from_universal_model = universal_model.reactions.get_by_id(mnxr)

        #cobrapy doc: "When copying a reaction, it is necessary to deepcopy the components so the list references are not carried over"
        mnxr_from_universal_model = copy.deepcopy(mnxr_from_universal_model)
        #Integrated model having reactions from both target and universal models
        target_model2.add_reaction(mnxr_from_universal_model)

    return target_model2


def get_unique_nonprod_monomers_list(options):

    unique_prod_monomers_list = []
    unique_nonprod_monomers_list = []

    for prod_monomers_list in options.prod_sec_met_dict.keys():
        for prod_monomer in options.prod_sec_met_dict[prod_monomers_list]:
            if prod_monomer not in unique_prod_monomers_list:
                unique_prod_monomers_list.append(prod_monomer)

    for nonprod_monomers_list in options.nonprod_sec_met_dict.keys():
        for nonprod_monomer in options.nonprod_sec_met_dict[nonprod_monomers_list]:
            if nonprod_monomer not in unique_nonprod_monomers_list and nonprod_monomer not in unique_prod_monomers_list:
                unique_nonprod_monomers_list.append(nonprod_monomer)

    logging.debug(unique_nonprod_monomers_list)
    return unique_nonprod_monomers_list


def add_transport_exchange_rxn_nonprod_monomer(target_model, nonprod_monomer, options):

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

    #Model reloading and overwrtting are necessary for model stability
    write_cobra_model_to_sbml_file(target_model_temp, '%s/3_temp_models/target_model_temp_%s.xml' %(options.output, nonprod_monomer))
    target_model_temp = create_cobra_model_from_sbml_file('%s/3_temp_models/target_model_temp_%s.xml' %(options.output, nonprod_monomer))

    return target_model_temp


def check_producibility_nonprod_monomer(cobra_model, nonprod_monomer):
    for rxn in cobra_model.reactions:
        rxn.objective_coefficient = 0

    cobra_model.reactions.get_by_id("Ex_"+nonprod_monomer).objective_coefficient = 1
    cobra_model.optimize()

    logging.debug(cobra_model.reactions.get_by_id("Ex_"+nonprod_monomer))
    logging.debug("Flux: %s" %cobra_model.solution.f)

    return cobra_model


def execute_gapfill(target_model_temp, nonprod_monomer, options):

    #Run gap-filling algorithm based on MILP in gurobipy
    obj = gapfilling_precursor()

    #Load merged model
    obj.load_cobra_model(target_model_temp)
    added_reaction = obj.fill_gap("Ex_"+nonprod_monomer, options.mnxr_unique_to_universal_model_list)

    return added_reaction


#Find gap-filling reactions that cause large biomass values and unrealistic fluxes for critial nutrients (e.g., O2, CO2 etc)
#Remove such gap-filling reactions from the list of reactions to be added to the model
def check_gapfill_rxn_biomass_effects(target_model, universal_model, added_reaction, options):

    added_reaction2 = copy.deepcopy(added_reaction)
    target_model_gapFilled = copy.deepcopy(target_model)

    for gapfill_rxn in added_reaction:
        target_model_gapFilled.add_reaction(universal_model.reactions.get_by_id(gapfill_rxn))

        write_cobra_model_to_sbml_file(target_model_gapFilled, "./%s/3_temp_models/target_model_gapFilled.xml" %options.output)
        target_model_gapFilled = create_cobra_model_from_sbml_file("./%s/3_temp_models/target_model_gapFilled.xml" %options.output)

        target_exrxnid_flux_dict = augPhase.get_exrxnid_flux(target_model_gapFilled, options.template_exrxnid_flux_dict)
        exrxn_flux_change_list = augPhase.check_exrxn_flux_direction(options.template_exrxnid_flux_dict, target_exrxnid_flux_dict)

        #Remove gap-filling reactions if they cause wrong flux values for nutrients transport
        #Removal of such reactions does not affect producibility of corresponding secondary metabolites, and even generates more realistic flux values
        if 'F' in exrxn_flux_change_list:
            target_model_gapFilled.remove_reactions(universal_model.reactions.get_by_id(gapfill_rxn))
            write_cobra_model_to_sbml_file(target_model_gapFilled, "./%s/3_temp_models/target_model_gapFilled.xml" %options.output)
            target_model_gapFilled = create_cobra_model_from_sbml_file("./%s/3_temp_models/target_model_gapFilled.xml" %options.output)

            logging.debug("Gap-filling reaction causing wrong fluxes: %s" %gapfill_rxn)
            added_reaction2.remove(gapfill_rxn)

    return added_reaction2


def add_gapfill_rxn_target_model(target_model, universal_model, added_reaction2, options):

    for gapfill_rxn in added_reaction2:
        target_model.add_reaction(universal_model.reactions.get_by_id(gapfill_rxn))

        logging.debug("Reaction added to the target_model: %s" %gapfill_rxn)

    #Cleanup of the final version of the target model
    #This function causes an error if executed after 'copy.deepcopy'ing the model
    prune_unused_metabolites(target_model)
    target_model_complete = copy.deepcopy(target_model)

    return target_model_complete
