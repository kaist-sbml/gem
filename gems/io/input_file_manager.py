
import glob
import logging
import os
import pickle
from Bio import SeqIO
from io_utils import (
    get_temp_fasta,
    get_features_from_gbk,
    get_features_from_fasta,
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
        #'3_primary_metabolic_model'
        options.outputfolder3 = os.path.join(options.outputfolder, folders[2])
        make_folder(options.outputfolder3)
        #'4_complete_model'
        options.outputfolder4 = os.path.join(options.outputfolder, folders[3])
        make_folder(options.outputfolder4)

    #'tmp_files'
    options.outputfolder5 = os.path.join(options.outputfolder, folders[4])
    make_folder(options.outputfolder5)


def check_input_filetype(options):
    logging.info('Input file: %s', options.input)
    input_ext = os.path.splitext(options.input)[1]

    if input_ext in ('.gbk', '.gb', '.genbank', '.gbf', '.gbff'):
        logging.debug("A GenBank file is found for input")
        return 'genbank'
    elif input_ext in ('.fa', '.fasta', '.fna', '.faa', '.fas'):
        logging.debug("A FASTA file is found for input")
        return 'fasta'


def get_target_genome_from_input(filetype, options):

    options.targetGenome_locusTag_aaSeq_dict = {}
    options.targetGenome_locusTag_ec_dict = {}
    options.targetGenome_locusTag_prod_dict = {}
    options.total_cluster = 0

    seq_records = list(SeqIO.parse(options.input, filetype))

    # len(seq_records) == 1: e.g., A complete bacterial genome (1 contig)
    # len(seq_records) > 1: e.g., An incomplete bacterial genome (multiple contigs)
    if len(seq_records) >= 1:
        if len(seq_records) == 1:
            logging.debug("One record is found in genome data")
        elif len(seq_records) > 1:
            logging.debug("Multiple records are found in genome data")

        if filetype == 'genbank':
            for seq_record in seq_records:
                get_features_from_gbk(seq_record, options)

        elif filetype == 'fasta':
            for seq_record in seq_records:
                get_features_from_fasta(seq_record, options)


    #Number of 'locus_tag's obtained above may be different from
    #the number directly obtained from genbank file
    #because some 'locus_tag's exit in the features 'tRNA' and 'rRNA',
    #which are not considered herein.
    #Same above comment for 'product's
    logging.debug(
                "len(options.targetGenome_locusTag_prod_dict.keys):%s"
                %len(options.targetGenome_locusTag_prod_dict.keys()))

    logging.debug(
                "len(options.targetGenome_locusTag_ec_dict.keys): %s"
                %len(options.targetGenome_locusTag_ec_dict.keys()))

    return seq_records


def get_target_genome_from_eficaz(options):

    #'1_EFICAz_results': following argument should not be changed
    gbk_file = glob.glob(os.path.join(options.outputfolder1, '*.gbk'))
    options.input = gbk_file[0]
    seq_record = get_target_genome_from_input('genbank', options)

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

