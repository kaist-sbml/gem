
import logging
from gmsm.primary_model.augPhase_utils import(
    get_targetGenome_locusTag_ec_nonBBH_dict,
    get_rxnid_info_dict_from_kegg,
    get_mnxr_list_from_modelPrunedGPR,
    get_mnxr_to_add_list,
    add_nonBBH_rxn,
    get_rxn_newComp_list_from_model,
    create_rxn_newComp,
    remove_inactive_rxn_newComp
)
from gmsm.primary_model.prunPhase_utils import (
    label_rxn_to_remove,
    prune_model,
    swap_locustag_with_homolog
)

def run_prunPhase(model, io_ns, config_ns, homology_ns, primary_model_ns):
    logging.info("Pruning phase starting..")
    logging.info("Removing reactions with nonhomologous genes from the template model..")
    label_rxn_to_remove(model, io_ns, homology_ns, primary_model_ns)

    modelPruned = prune_model(model, config_ns, primary_model_ns)

    logging.info("Correcting GPR associations in the template model..")
    modelPrunedGPR = swap_locustag_with_homolog(modelPruned, homology_ns)

    return modelPrunedGPR


def run_augPhase(modelPrunedGPR, run_ns, io_ns, config_ns, homology_ns, primary_model_ns):
    logging.info("Augmentation phase starting..")
    logging.info("Creating various dict data for the nonBBH gene-associted reactions..")
    logging.info("(time-consuming)")

    get_targetGenome_locusTag_ec_nonBBH_dict(io_ns, homology_ns, primary_model_ns)

    #Two nested functions
    #def get_rxnid_from_ECNumber(enzymeEC):
    #def get_rxnInfo_from_rxnid(rxnid):
    get_rxnid_info_dict_from_kegg(io_ns, config_ns, primary_model_ns)

    get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, io_ns, primary_model_ns)

    logging.info("Adding the nonBBH gene-associated reactions..")
    logging.info("(time-consuming)")

    get_mnxr_to_add_list(io_ns, primary_model_ns)

    target_model = add_nonBBH_rxn(modelPrunedGPR, io_ns, config_ns, primary_model_ns)

    if run_ns.comp:
        logging.info("Adding reactions with the new compartment..")
        logging.info("(time-consuming)")

        rxn_newComp_list = get_rxn_newComp_list_from_model(target_model, io_ns)

        target_model, added_rxn_newComp_list = \
                create_rxn_newComp(rxn_newComp_list, target_model, io_ns, primary_model_ns)

        target_model = remove_inactive_rxn_newComp(
                added_rxn_newComp_list, target_model, io_ns, primary_model_ns)

    return target_model
