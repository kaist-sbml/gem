'''
2014-2015
Hyun Uk Kim, Tilmann Weber, Kyu-Sang Hwang and Jae Yong Ryu
'''

#Wildcard imports should never be used in production code.
#import argparse
import logging
#import multiprocessing
#import os
import pickle
#import sys
#import time

#from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
#from argparse import Namespace
#from modeling import prunPhase
#from modeling import augPhase
#from modeling import sec_met_rxn_generation
#from modeling.gapfilling import gapfill_network_manipulation
from io_utils import (
    get_temp_fasta,
    get_target_gbk,
    get_targetGenomeInfo,
    get_target_fasta
)


def get_genome_files(options):
    logging.info("Reading input genome files..")

    logging.info("Looking for a gbk file of a template model genome..")
    get_temp_fasta(options)

    logging.info("Looking for a gbk file of a target genome..")
    get_target_gbk(options)

    logging.info("Reading genbank file of the target genome..")
    get_targetGenomeInfo(options, 'genbank')

    logging.info("Looking for a fasta file of a target genome..")
    get_target_fasta(options)


#For model augmentation phase in both primary and secondary modeling
def get_pickles_add_rxn(options):
    logging.info("Loading pickle files necessary for the model augmentation phase..")

    bigg_mnxr_dict = pickle.load(open('./modeling/io/data/input2/bigg_mnxr_dict.p','rb'))
    options.bigg_mnxr_dict = bigg_mnxr_dict
    kegg_mnxr_dict = pickle.load(open('./modeling/io/data/input2/kegg_mnxr_dict.p','rb'))
    options.kegg_mnxr_dict = kegg_mnxr_dict
    mnxr_kegg_dict = pickle.load(open('./modeling/io/data/input2/mnxr_kegg_dict.p','rb'))
    options.mnxr_kegg_dict = mnxr_kegg_dict
    mnxr_rxn_dict = pickle.load(open('./modeling/io/data/input2/mnxr_rxn_dict.p','rb'))
    options.mnxr_rxn_dict = mnxr_rxn_dict

    bigg_mnxm_compound_dict = pickle.load(open('./modeling/io/data/input2/bigg_mnxm_compound_dict.p','rb'))
    options.bigg_mnxm_compound_dict = bigg_mnxm_compound_dict
    mnxm_bigg_compound_dict = pickle.load(open('./modeling/io/data/input2/mnxm_bigg_compound_dict.p','rb'))
    options.mnxm_bigg_compound_dict = mnxm_bigg_compound_dict
    kegg_mnxm_compound_dict = pickle.load(open('./modeling/io/data/input2/kegg_mnxm_compound_dict.p','rb'))
    options.kegg_mnxm_compound_dict = kegg_mnxm_compound_dict
    mnxm_kegg_compound_dict = pickle.load( open('./modeling/io/data/input2/mnxm_kegg_compound_dict.p','rb'))
    options.mnxm_kegg_compound_dict = mnxm_kegg_compound_dict

    mnxm_compoundInfo_dict = pickle.load(open('./modeling/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
    options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(options.input1),'rb'))
    options.template_exrxnid_flux_dict = template_exrxnid_flux_dict


#For model pruning phase
#Only model file is not saved in Namespace
def get_pickles_prunPhase(options):
    logging.info("Loading pickle files of the parsed template model and its relevant genbank data..")
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    tempModel_biggRxnid_locusTag_dict = pickle.load(open('%s/tempModel_biggRxnid_locusTag_dict.p' %(options.input1),'rb'))
    options.tempModel_biggRxnid_locusTag_dict = tempModel_biggRxnid_locusTag_dict

    return model

