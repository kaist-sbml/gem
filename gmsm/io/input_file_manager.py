
import glob
import logging
import os
import pickle
import re
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
            '3_primary_metabolic_model', '4_complete_model',
            'tmp_model_files', 'tmp_data_files']

    # Second if statement is to keep "-o ./test" from creating "tes", not "test"
    if '/' in options.outputfolder and './' not in options.outputfolder:
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

    #'tmp_model_files'
    options.outputfolder5 = os.path.join(options.outputfolder, folders[4])
    make_folder(options.outputfolder5)

    #'tmp_data_files'
    options.outputfolder6 = os.path.join(options.outputfolder, folders[5])
    make_folder(options.outputfolder6)


def show_input_options(options):
    logging.debug("input_file: %s", options.input)
    logging.debug("outputfolder: %s", options.outputfolder)
    logging.debug("template_model_organism: %s", options.orgName)
    logging.debug("eficaz: %s", options.eficaz)
    logging.debug("primary_metabolic_modeling: %s", options.pmr_generation)
    logging.debug("secondary_metabolic_modeling: %s", options.smr_generation)
    logging.debug("eficaz_file: %s", options.eficaz_file)
    logging.debug("compartment_file: %s", options.comp)


def check_input_filetype(options):
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


def get_eficaz_file(options):

    logging.info("Reading EFICAz output file..")

#    EC4Info = {}
#    EC3Info = {}

    try:
        f = open(options.eficaz_file,"r")
    except OSError as e:
         logging.error("No EFICAz output file %s found", options.eficaz_file)
         #continue
    except IOError as e:
         logging.error("No EFICAz output file %s found", options.eficaz_file)
         #continue

    for line in f.read().splitlines():
        (locustag, eficazResultString) = line.split(',', 1)
        eficazResultString = eficazResultString.strip()
        if eficazResultString == 'No EFICAz EC assignment':
            continue

        #if eficazResultString.strip().startswith("3EC"):
        #    r = re.match('3EC: (\d+\.\d+\.\d+), (.*)', eficazResultString)
        #    if r:
        #        EC = r.group(1) + ".-"
        #        #ECDesc = r.group(2)
        #        if not options.targetGenome_locusTag_ec_dict.has_key(locustag):
        #            options.targetGenome_locusTag_ec_dict[locustag] = []
        #            #EC3Info[locustag] = []
        #        options.targetGenome_locusTag_ec_dict[locustag].append(EC)
        #        #EC3Info[locustag].append(ECDesc)
        #        continue

        if eficazResultString.strip().startswith("4EC"):
            r = re.match('4EC: (\d+\.\d+\.\d+\.\d+), (.*)', eficazResultString)
            if r:
                EC = r.group(1)
                #ECDesc = r.group(2)
                if not options.targetGenome_locusTag_ec_dict.has_key(locustag):
                    options.targetGenome_locusTag_ec_dict[locustag] = []
                    #EC4Info[locustag] = []
                options.targetGenome_locusTag_ec_dict[locustag].append(EC)
                #EC4Info[locustag].append(ECDesc)
                continue

                logging.debugging(
                        "Locus tag: %s; EC number with 4 digits: %s", locustag, EC)
    f.close()

    logging.debug("len(options.targetGenome_locusTag_ec_dict.keys): %s",
                  len(options.targetGenome_locusTag_ec_dict.keys()))


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

    bigg_mnxr_dict = pickle.load(open('./gmsm/io/data/input2/bigg_mnxr_dict.p','rb'))
    options.bigg_mnxr_dict = bigg_mnxr_dict

    mnxm_compoundInfo_dict = pickle.load(open('./gmsm/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
    options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    mnxr_kegg_dict = pickle.load(open('./gmsm/io/data/input2/mnxr_kegg_dict.p','rb'))
    options.mnxr_kegg_dict = mnxr_kegg_dict

    mnxref = pickle.load(open('./gmsm/io/data/input2/MNXref.p','rb'))
    options.mnxref = mnxref

    template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(options.input1),'rb'))
    options.template_exrxnid_flux_dict = template_exrxnid_flux_dict


def get_locustag_comp_dict(options):

    logging.info("Reading file on subcellular localizations (compartments)..")

    options.locustag_comp_dict = {}

    try:
        f = open(options.comp,"r")
    except OSError as e:
         logging.error("No file %s (subcellular localizations (compartments)) found",
                        options.comp)
         #continue
    except IOError as e:
         logging.error("No file %s (subcellular localizations (compartments)) found",
                        options.comp)
         #continue

    for line in f.read().splitlines():
        (locustag, comp) = line.split('\t', 1)

        locustag = locustag.strip()
        comp = comp.strip()
        logging.debug("Locus tag: %s; Predicted compartment: %s", locustag, comp)

        if locustag not in options.locustag_comp_dict:
            options.locustag_comp_dict[locustag] = [comp]
        else:
            options.locustag_comp_dict[locustag].append(comp)

    f.close()

    logging.debug("len(options.locustag_comp_dict.keys): %s",
                  len(options.locustag_comp_dict.keys()))

