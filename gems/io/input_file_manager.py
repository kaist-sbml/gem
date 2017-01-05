
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import pickle
from io_utils import (
    get_temp_fasta,
    get_targetGenomeInfo,
    get_target_fasta
)


def get_genome_files(options):
    logging.info("Reading input genome files..")

    logging.info("Reading a genbank file of a target genome..")
    get_targetGenomeInfo(options, 'genbank')

    #Following data are needed only for primary metabolic modeling
    if options.pmr_generation:
        logging.info("Looking for a fasta file of a target genome..")
        get_target_fasta(options)

        logging.info("Looking for a fasta file of template model genes..")
        get_temp_fasta(options)


#For model pruning phase
#Only model file is not saved in Namespace
def get_pickles_prunPhase(options):
    logging.info("Loading pickle files associated with a template model..")
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    tempModel_biggRxnid_locusTag_dict = pickle.load(open(
        '%s/tempModel_biggRxnid_locusTag_dict.p' %(options.input1),'rb'))
    options.tempModel_biggRxnid_locusTag_dict = tempModel_biggRxnid_locusTag_dict

    return model


#For model augmentation phase in both primary and secondary modeling
def get_pickles_augPhase(options):
    logging.info("Loading pickle files necessary for the model augmentation phase..")

    bigg_mnxr_dict = pickle.load(open('./gems/io/data/input2/bigg_mnxr_dict.p','rb'))
    options.bigg_mnxr_dict = bigg_mnxr_dict
    kegg_mnxr_dict = pickle.load(open('./gems/io/data/input2/kegg_mnxr_dict.p','rb'))
    options.kegg_mnxr_dict = kegg_mnxr_dict
    mnxr_kegg_dict = pickle.load(open('./gems/io/data/input2/mnxr_kegg_dict.p','rb'))
    options.mnxr_kegg_dict = mnxr_kegg_dict
    mnxr_rxn_dict = pickle.load(open('./gems/io/data/input2/mnxr_rxn_dict.p','rb'))
    options.mnxr_rxn_dict = mnxr_rxn_dict

    bigg_mnxm_compound_dict = pickle.load(open('./gems/io/data/input2/bigg_mnxm_compound_dict.p','rb'))
    options.bigg_mnxm_compound_dict = bigg_mnxm_compound_dict
    mnxm_bigg_compound_dict = pickle.load(open('./gems/io/data/input2/mnxm_bigg_compound_dict.p','rb'))
    options.mnxm_bigg_compound_dict = mnxm_bigg_compound_dict

    mnxm_compoundInfo_dict = pickle.load(open('./gems/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
    options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(options.input1),'rb'))
    options.template_exrxnid_flux_dict = template_exrxnid_flux_dict

