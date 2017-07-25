
import copy
import logging
from cobra import Reaction, Metabolite
from gmsm import utils

def get_unique_nonprod_monomers_list(options):

    unique_prod_monomers_list = []
    unique_nonprod_monomers_list = []

    for prod_monomers_list in options.prod_sec_met_dict:
        for prod_monomer in options.prod_sec_met_dict[prod_monomers_list]:
            if prod_monomer not in unique_prod_monomers_list:
                unique_prod_monomers_list.append(prod_monomer)

    for nonprod_monomers_list in options.nonprod_sec_met_dict:
        for nonprod_monomer in options.nonprod_sec_met_dict[nonprod_monomers_list]:
            if nonprod_monomer not in unique_nonprod_monomers_list \
                and nonprod_monomer not in unique_prod_monomers_list:
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
    utils.stabilize_model(target_model_temp, options.outputfolder5, nonprod_monomer)

    return target_model_temp


def check_producibility_nonprod_monomer(cobra_model, nonprod_monomer):
    for rxn in cobra_model.reactions:
        rxn.objective_coefficient = 0

    cobra_model.reactions.get_by_id("Ex_"+nonprod_monomer).objective_coefficient = 1
    flux_dist = cobra_model.optimize()

    logging.debug("%s; Flux value: %f",
                    cobra_model.reactions.get_by_id("Ex_"+nonprod_monomer),
                    flux_dist.objective_value)

    return cobra_model, flux_dist


#Find gap-filling reactions that cause large biomass values
#and unrealistic fluxes for critial nutrients (e.g., O2, CO2 etc)
#Remove such gap-filling reactions from the list of reactions to be added to the model
def check_gapfill_rxn_biomass_effects(target_model, universal_model,
                                     gapfill_rxns, options):

    gapfill_rxns2 = copy.deepcopy(gapfill_rxns)
    target_model_gapFilled = copy.deepcopy(target_model)

    for gapfill_rxn in gapfill_rxns:
        target_model_gapFilled.add_reaction(
                universal_model.reactions.get_by_id(gapfill_rxn))

        utils.stabilize_model(target_model_gapFilled, options.outputfolder5, '')

        target_exrxnid_flux_dict = utils.get_exrxnid_flux(
                target_model_gapFilled, options.template_exrxnid_flux_dict)
        exrxn_flux_change_list = utils.check_exrxn_flux_direction(
                options.template_exrxnid_flux_dict, target_exrxnid_flux_dict)

        #Remove gap-filling reactions
        #if they cause wrong flux values for nutrients transport
        #Removal of such reactions does not affect
        #producibility of corresponding secondary metabolites,
        #and even generates more realistic flux values
        if 'F' in exrxn_flux_change_list:
            target_model_gapFilled.remove_reactions(
                    universal_model.reactions.get_by_id(gapfill_rxn))

            utils.stabilize_model(target_model_gapFilled, options.outputfolder5, '')

            logging.debug("Gap-filling reaction causing wrong fluxes: %s"
                            %str(gapfill_rxn))

            gapfill_rxns2.remove(gapfill_rxn)

    return gapfill_rxns2


def add_gapfill_rxn_target_model(target_model, universal_model, gapfill_rxns2, options):

    for gapfill_rxn in gapfill_rxns2:
        target_model.add_reaction(universal_model.reactions.get_by_id(gapfill_rxn))

        logging.debug("Reaction added to the target_model: %s" %gapfill_rxn)

    target_model_complete = copy.deepcopy(target_model)

    return target_model_complete
