'''
2014-2015
Hyun Uk Kim, Tilmann Weber, Kyu-Sang Hwang and Jae Yong Ryu
'''

#Wildcard imports should never be used in production code.
#import argparse
import logging
#import multiprocessing
#import os
#import pickle
#import sys
#import time

#from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
#from argparse import Namespace
#from modeling import prunPhase
#from modeling import augPhase
#from modeling import sec_met_rxn_generation
#from modeling.gapfilling import gapfill_network_manipulation
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
    modelPruned = pruneModel(model, options, 'gurobi')

    logging.info("Correcting GPR associations in the template model..")
    modelPrunedGPR = swap_locusTag_tempModel(modelPruned, options)

    return modelPrunedGPR

