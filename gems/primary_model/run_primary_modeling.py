
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
from augPhase_utils import(
    get_targetGenome_locusTag_ec_nonBBH_dict,
    get_rxnid_info_dict_from_kegg,
    get_mnxr_list_from_modelPrunedGPR,
    get_rxnid_to_add_list,
    get_mnxr_to_add_list,
    get_rxnid_mnxm_coeff_dict,
    add_nonBBH_rxn
)
from prunPhase_utils import (
    labelRxnToRemove,
    pruneModel,
    swap_locusTag_tempModel
)

def run_prunPhase(model, options):
    logging.info("Pruning phase starting..")
    logging.info("Labeling reactions with nonhomologous genes to remove from the template model..")
    labelRxnToRemove(model, options)

    logging.info("Removing reactions with nonhomologous genes from the template model..")
    modelPruned = pruneModel(model, options)

    logging.info("Correcting GPR associations in the template model..")
    modelPrunedGPR = swap_locusTag_tempModel(modelPruned, options)

    return modelPrunedGPR


def run_augPhase(modelPrunedGPR, options):
    logging.info("Augmentation phase starting..")
    logging.info("Creating various dictionary files for the nonBBH gene-associted reactions... (time-consuming)")

    get_targetGenome_locusTag_ec_nonBBH_dict(options)

    #Two nested functions
    #def get_rxnid_from_ECNumber(enzymeEC):
    #def get_rxnInfo_from_rxnid(rxnid):
    get_rxnid_info_dict_from_kegg(options)

    get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, options)

    logging.info("Adding the nonBBH gene-associated reactions... (time-consuming)")
    get_rxnid_to_add_list(options)

    get_mnxr_to_add_list(options)

    get_rxnid_mnxm_coeff_dict(options)

    target_model = add_nonBBH_rxn(modelPrunedGPR, options)

    return target_model
