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
from augPhase_utils import(
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

    target_model = augPhase.add_nonBBH_rxn(modelPrunedGPR, options)

    return target_model

