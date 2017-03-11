
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import glob
import logging
import os
import pickle
from io_utils import (
    get_temp_fasta,
    get_targetGenomeInfo,
    get_target_fasta
)


def make_folder(folder):
    if not os.path.isdir(folder):
        os.makedirs(folder)


def setup_outputfolders(options):
    folders = ['1_EFICAz_results', '2_blastp_results',
            '3_primary_metabolic_model', '4_complete_model', 'tmp_files']

    if '/' in options.outputfolder:
        options.outputfolder = options.outputfolder[:-1]

    if options.eficaz:
        #'1_EFICAz_results'
        options.outputfolder1 = os.path.join(options.outputfolder, folders[0])
        make_folder(options.outputfolder1)
    if options.pmr_generation:
        #'2_blastp_results'
        options.outputfolder2 = os.path.join(options.outputfolder, folders[1])
        make_folder(options.outputfolder2)
        #'3_primary_metabolic_model'
        options.outputfolder3 = os.path.join(options.outputfolder, folders[2])
        make_folder(options.outputfolder3)
    if options.smr_generation:
        #'4_complete_model'
        options.outputfolder4 = os.path.join(options.outputfolder, folders[3])
        make_folder(options.outputfolder4)

    #'tmp_files'
    options.outputfolder5 = os.path.join(options.outputfolder, folders[4])
    make_folder(options.outputfolder5)


def check_input_filetype(options):
    logging.info('Input file: %s', options.input)
    input_ext = os.path.splitext(options.input)[1]

    if input_ext in ('.gbk', '.gb', '.genbank', '.gbff'):
        logging.debug("GenBank file format is supported")
    elif input_ext in ('.fa', '.fasta', '.fna', '.faa', '.fas'):
        logging.warning("FASTA file format is not currently supported")
        if options.eficaz:
            options.eficaz = False
        elif options.pmr_generation:
            options.pmr_generation = False
        elif options.smr_generation:
            options.smr_generation = False


def get_target_gbk(options):
    logging.info("Reading input genome files..")
    logging.info("Reading a genbank file of a target genome..")

    try:
        gbk_file = glob.glob(os.path.join(options.outputfolder1, '*.gbk'))
        seq_record = get_targetGenomeInfo(gbk_file[0], options, 'genbank')
    except:
        seq_record = get_targetGenomeInfo(options.input, options, 'genbank')

    return seq_record


def get_fasta_files(options):
    #Following data are needed only for primary metabolic modeling
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

    mnxm_compoundInfo_dict = pickle.load(open('./gems/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
    options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    mnxr_kegg_dict = pickle.load(open('./gems/io/data/input2/mnxr_kegg_dict.p','rb'))
    options.mnxr_kegg_dict = mnxr_kegg_dict

    mnxref = pickle.load(open('./gems/io/data/input2/MNXref.p','rb'))
    options.mnxref = mnxref

    template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(options.input1),'rb'))
    options.template_exrxnid_flux_dict = template_exrxnid_flux_dict

