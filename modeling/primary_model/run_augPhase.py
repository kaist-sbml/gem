
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
rom augPhase_utils import(
    get_targetGenome_locusTag_ec_nonBBH_dict,
    make_all_rxnInfo_fromRefSeq,
    get_mnxr_list_from_modelPrunedGPR,
    check_existing_rxns,
    get_mnxr_using_kegg,
    extract_rxn_mnxm_coeff,
    add_nonBBH_rxn
)


def run_augPhase(modelPrunedGPR, options):
    logging.info("Augmentation phase starting..")
    logging.info("Creating various dictionary files for the nonBBH gene-associted reactions...")

    get_targetGenome_locusTag_ec_nonBBH_dict(options)

    #Two nested functions
    #def get_rxnid_from_ECNumber(enzymeEC):
    #def get_rxnInfo_from_rxnid(rxnid):
    make_all_rxnInfo_fromRefSeq(options)

    get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, options)

    logging.info("Adding the nonBBH gene-associated reactions...")
    check_existing_rxns(options)

    get_mnxr_using_kegg(options)

    extract_rxn_mnxm_coeff(options)

    target_model = add_nonBBH_rxn(modelPrunedGPR, options)

    return target_model

