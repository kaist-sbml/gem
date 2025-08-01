
import glob
import logging
import os
import pickle
import re
from Bio import SeqIO
from gmsm.io.io_utils import (
    get_temp_fasta,
    get_features_from_gbk,
    get_features_from_fasta,
    get_target_fasta
)


def make_folder(folder):
    if not os.path.isdir(folder):
        os.makedirs(folder)


def setup_outputfolders(run_ns, io_ns):
    folders = ['1_EFICAz_results', '2_blastp_results',
            '3_primary_metabolic_model', '4_complete_model',
            'tmp_model_files', 'tmp_data_files']

    # Keep "-o test/test" from creating "test/tes", not "test/test"
    if '/' in run_ns.outputfolder[-1]:
        run_ns.outputfolder = run_ns.outputfolder[:-1]

    if run_ns.eficaz:
        #'1_EFICAz_results'
        io_ns.outputfolder1 = os.path.join(run_ns.outputfolder, folders[0])
        make_folder(io_ns.outputfolder1)
    if run_ns.pmr_generation:
        #'2_blastp_results'
        io_ns.outputfolder2 = os.path.join(run_ns.outputfolder, folders[1])
        make_folder(io_ns.outputfolder2)
        #'3_primary_metabolic_model'
        io_ns.outputfolder3 = os.path.join(run_ns.outputfolder, folders[2])
        make_folder(io_ns.outputfolder3)
    if run_ns.smr_generation:
        #'3_primary_metabolic_model'
        io_ns.outputfolder3 = os.path.join(run_ns.outputfolder, folders[2])
        make_folder(io_ns.outputfolder3)
        #'4_complete_model'
        io_ns.outputfolder4 = os.path.join(run_ns.outputfolder, folders[3])
        make_folder(io_ns.outputfolder4)

    #'tmp_model_files'
    io_ns.outputfolder5 = os.path.join(run_ns.outputfolder, folders[4])
    make_folder(io_ns.outputfolder5)

    #'tmp_data_files'
    io_ns.outputfolder6 = os.path.join(run_ns.outputfolder, folders[5])
    make_folder(io_ns.outputfolder6)


def show_input_options(run_ns):

    logging.debug("input_file: %s", run_ns.input)
    logging.debug("outputfolder: %s", run_ns.outputfolder)
    logging.debug("template_model_organism: %s", run_ns.orgName)
    logging.debug("eficaz: %s", run_ns.eficaz)
    logging.debug("primary_metabolic_modeling: %s", run_ns.pmr_generation)
    logging.debug("secondary_metabolic_modeling: %s", run_ns.smr_generation)
    logging.debug("ec_number_file: %s", run_ns.ec_file)
    logging.debug("compartment_file: %s", run_ns.comp)


def check_input_filetype(run_ns):

    input_ext = os.path.splitext(run_ns.input)[1]

    if input_ext in ('.gbk', '.gb', '.genbank', '.gbf', '.gbff'):
        logging.debug("A GenBank file is found for input")
        return 'genbank'
    elif input_ext in ('.fa', '.fasta', '.fna', '.faa', '.fas'):
        logging.debug("A FASTA file is found for input")
        return 'fasta'


def get_target_genome_from_input(filetype, run_ns, io_ns):
    io_ns.targetGenome_locusTag_aaSeq_dict = {}
    io_ns.targetGenome_locusTag_ec_dict = {}
    io_ns.targetGenome_locusTag_prod_dict = {}
    io_ns.seq_record_BGC_num_lists = []
    io_ns.total_region = 0

    seq_records = list(SeqIO.parse(run_ns.input, filetype))
    # len(seq_records) == 1: e.g., A complete bacterial genome (1 contig)
    # len(seq_records) > 1: e.g., An incomplete bacterial genome (multiple contigs)
    if len(seq_records) >= 1:
        if len(seq_records) == 1:
            logging.debug("One record is found in genome data")
        elif len(seq_records) > 1:
            logging.debug("Multiple records are found in genome data")

        # Ignore existing annotations of EC numbers in an input gbk file as they are from a different source.
        if run_ns.eficaz or run_ns.ec_file:
            logging.info("Ignoring EC annotations from input gbk file")
        else:
            logging.info("Using EC annotations from input gbk file")

        if filetype == 'genbank':
            for seq_record in seq_records:
                get_features_from_gbk(seq_record, run_ns, io_ns)
            if io_ns.total_region == 0:
                logging.debug("This gbk file needs to be processed by antiSMASH or does not contain any BGC")

        elif filetype == 'fasta':
            for seq_record in seq_records:
                get_features_from_fasta(seq_record, io_ns)

    #Number of 'locus_tag's obtained above may be different from
    #the number directly obtained from genbank file
    #because some 'locus_tag's exit in the features 'tRNA' and 'rRNA',
    #which are not considered herein.
    #Same above comment for 'product's
    logging.debug(
                "len(io_ns.targetGenome_locusTag_prod_dict.keys):%s"
                %len(io_ns.targetGenome_locusTag_prod_dict.keys()))

    logging.debug(
                "len(io_ns.targetGenome_locusTag_ec_dict.keys): %s"
                %len(io_ns.targetGenome_locusTag_ec_dict.keys()))


def get_ec_file(run_ns, io_ns):

    logging.info("Reading EC number prediction file..")
    
    with open(run_ns.ec_file, 'r') as f:
        header = next(f)
        for line in f:
            locustag, pred = line.strip().split('\t')
            if not pred.startswith('EC:'):
                continue
            ec_num = pred[len('EC:'):]
            parts = ec_num.split('.')
            # select EC number with 4 digits
            if len(parts) == 4 and all(part.isdigit() for part in parts):
                if not locustag in io_ns.targetGenome_locusTag_ec_dict.keys():
                    io_ns.targetGenome_locusTag_ec_dict[locustag] = []
                io_ns.targetGenome_locusTag_ec_dict[locustag].append(ec_num)

                logging.debug(
                        "Locus tag: %s; EC number with 4 digits: %s", locustag, ec_num)

    logging.debug("len(io_ns.targetGenome_locusTag_ec_dict.keys): %s",
                  len(io_ns.targetGenome_locusTag_ec_dict.keys()))


def get_fasta_files(run_ns, io_ns):
    #Following data are needed only for primary metabolic modeling
    logging.info("Looking for a fasta file of a target genome..")
    get_target_fasta(io_ns)

    logging.info("Looking for a fasta file of template model genes..")
    get_temp_fasta(run_ns, io_ns)


#For model pruning phase
#Only model file is not saved in Namespace
def get_pickles_prunPhase(io_ns):
    logging.info("Loading pickle files associated with a template model..")
    model = pickle.load(open('%s/model.p' %(io_ns.input1),'rb'))
    tempModel_biggRxnid_locusTag_dict = pickle.load(open(
        '%s/tempModel_biggRxnid_locusTag_dict.p' %(io_ns.input1),'rb'))
    io_ns.tempModel_biggRxnid_locusTag_dict = tempModel_biggRxnid_locusTag_dict

    return model


#For model augmentation phase in both primary and secondary modeling
def get_pickles_augPhase(io_ns):
    logging.info("Loading pickle files necessary for the model augmentation phase..")

    bigg_mnxr_dict = pickle.load(open('./gmsm/io/data/input2/bigg_mnxr_dict.p','rb'))
    io_ns.bigg_mnxr_dict = bigg_mnxr_dict

    mnxm_compoundInfo_dict = pickle.load(open('./gmsm/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
    io_ns.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    mnxr_kegg_dict = pickle.load(open('./gmsm/io/data/input2/mnxr_kegg_dict.p','rb'))
    io_ns.mnxr_kegg_dict = mnxr_kegg_dict

    mnxref = pickle.load(open('./gmsm/io/data/input2/MNXref.p','rb'))
    io_ns.mnxref = mnxref

    template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(io_ns.input1),'rb'))
    io_ns.template_exrxnid_flux_dict = template_exrxnid_flux_dict


def get_locustag_comp_dict(run_ns, io_ns):

    logging.info("Reading file on subcellular localizations (compartments)..")

    io_ns.locustag_comp_dict = {}

    try:
        f = open(run_ns.comp,"r")
    except OSError as e:
         logging.error("No file %s (subcellular localizations (compartments)) found",
                        run_ns.comp)
         #continue
    except IOError as e:
         logging.error("No file %s (subcellular localizations (compartments)) found",
                        run_ns.comp)
         #continue

    for line in f.read().splitlines():
        (locustag, comp) = line.split('\t', 1)

        locustag = locustag.strip()
        comp = comp.strip()
        logging.debug("Locus tag: %s; Predicted compartment: %s", locustag, comp)

        if locustag not in io_ns.locustag_comp_dict:
            io_ns.locustag_comp_dict[locustag] = [comp]
        else:
            io_ns.locustag_comp_dict[locustag].append(comp)

    f.close()

    logging.debug("len(io_ns.locustag_comp_dict.keys): %s",
                  len(io_ns.locustag_comp_dict.keys()))

