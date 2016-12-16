
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
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

