# -*- coding: utf-8 -*-

import argparse
import ast
import cobra
import glob
import os
import logging
import pyparsing
import pickle
import subprocess
import sys
import urllib2
from Bio import Entrez, SeqIO
from input2_manager import ParseMNXref
from os.path import join, abspath, dirname

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import gems

def get_options():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    group = parser.add_argument_group(
        "Provide:",
        "1) BiGG Models ID (e.g., iAF1260) only or "
        "2) Genbank accession number (e.g., NC_003888.3) with a relevant SBML file or"
        "3) Genome file (GenBank or FASTA) with a relevant SBML file."
        "For the options 2) and 3), place files in 'scripts/input1_data/[organism-specific folder'].")

    group.add_argument('-m', '--model',
            dest='model',
            help = "Specify BiGG ID of a metabolic model to prepare as a template model")
    group.add_argument('-a', '--acc_number',
            dest='acc_number',
            help = "Specify an organism's Genbank accession number\n"
            "Also provide a relevant SBML file")
    group.add_argument('-g', '--genome',
            dest='genome',
            help = "Provide an organism's genome file (Genbank or FASTA)\n"
            "Also provide a relevant SBML file")
    group.add_argument('-f', '--folder',
            dest='folder',
            help = "Specify an output folder name with the KEGG organism code\n")

    options = parser.parse_args()
    logging.debug(options)

    return options

def get_output_dirs(options):

    if options.folder:
        input1_dir = join(os.pardir, 'gems', 'io', 'data', 'input1', options.folder)
        if not os.path.isdir(input1_dir):
            os.makedirs(input1_dir)

        input1_tmp_dir = join(dirname(abspath(__file__)), 'input1_data', options.folder)
        if not os.path.isdir(input1_tmp_dir):
            os.makedirs(input1_tmp_dir)

    return input1_dir, input1_tmp_dir


def download_model_from_biggDB(input1_tmp_dir, options):
    model_file = ''.join([options.model, '.xml'])
    url = ''.join(['http://bigg.ucsd.edu/static/models/', model_file])
    logging.debug('URL for downloading a model from the BiGG Models:')
    logging.debug(url)

    model = urllib2.urlopen(url).read()

    with open(join(input1_tmp_dir, model_file), 'wb') as f:
        f.write(model)

    model = cobra.io.read_sbml_model(join(input1_tmp_dir, model_file))
    model = gems.utils.stabilize_model(model, input1_tmp_dir, options.model)

    if len(model.reactions) > 1:
        logging.debug('%s downloaded successfully', options.model)
    else:
        logging.debug('%s NOT downloaded successfully', options.model)

    return model


def get_model_details(options):
    model_info_dict = {}

    url = ''.join(['http://bigg.ucsd.edu/api/v2/models/', options.model])
    logging.debug('URL for accessing model details the BiGG Models:')
    logging.debug(url)

    model_info = urllib2.urlopen(url).read()

    # 'null' causes "ValueError: malformed string"
    if 'null' in model_info:
        model_info = model_info.replace('null', '""')

    model_info_dict = ast.literal_eval(model_info)

    if not model_info_dict['organism']:
        model_info_dict['organism'] = raw_input('Organism name?')

        if not model_info_dict['organism']:
            logging.error("Organism name ('model_info_dict['organism']') is not provided")
            sys.exit(1)

    if not model_info_dict['genome_name']:
        model_info_dict['genome_name'] = raw_input('Genome name?')

        if not model_info_dict['genome_name']:
            logging.error("Genome name ('model_info_dict['genome_name']') is not provided")
            sys.exit(1)

    logging.debug('%s details:', options.model)
    logging.debug('model_bigg_id: %s', model_info_dict['model_bigg_id'])
    logging.debug('organism: %s', model_info_dict['organism'])
    logging.debug('genome_name: %s', model_info_dict['genome_name'])
    logging.debug('gene_count: %s', model_info_dict['gene_count'])

    return model_info_dict


def prepare_nonstd_model(input1_tmp_dir, options):
    mnx_parser = ParseMNXref()
    bigg_old_new_dict = mnx_parser.fix_legacy_id_using_BiGGModels()

    sbml_list = glob.glob(join(input1_tmp_dir, '*.xml'))
    logging.debug('Model found: %s', sbml_list)
    # This considers 'fix_legacy_id'
    model = cobra.io.read_legacy_sbml(join(input1_tmp_dir, sbml_list[0]))

    for i in range(len(model.metabolites)):
        metab = model.metabolites[i]
        if metab.id in bigg_old_new_dict:
            metab.id = bigg_old_new_dict[metab.id]

    for j in range(len(model.reactions)):
        rxn = model.reactions[j]

        # NOTE:
        #See '\BiGG\170405\nar-02327-data-e-2015-File017_All modifications_KHU_v2'
        #This is for the 'sco' template model.
        if rxn.id in bigg_old_new_dict and \
                rxn.id != 'FACOAL80' and \
                rxn.id != 'FE3abc':

            # Otherwise a duplicate reaction can be inserted: e.g., ME1 -> ME2
            if bigg_old_new_dict[rxn.id] not in model.reactions:
                logging.debug('Reaction: %s -> %s ', rxn.id, bigg_old_new_dict[rxn.id])
                rxn.id = bigg_old_new_dict[rxn.id]

        if rxn.id == 'THRPS':
            rxn.id = 'LTHRK'
            logging.debug('Reaction: %s -> %s ', 'THRPS', rxn.id)

    model = gems.utils.stabilize_model(model, input1_tmp_dir, '')

    model_info_dict = {}
    model_info_dict['genome_name'] = options.acc_number

    return model, model_info_dict


def download_gbk_from_ncbi(input1_tmp_dir, model_info_dict):
    gbk_file = ''.join([model_info_dict['genome_name'], '.gb'])
    Entrez.email = "ehukim@kaist.ac.kr"

    handle = Entrez.efetch(db='nucleotide',
            id=model_info_dict['genome_name'], rettype='gbwithparts', retmode='text')

    seq_record = handle.read()

    with open(join(input1_tmp_dir, gbk_file), 'wb') as f:
        f.write(seq_record)

    return gbk_file


def get_tempGenome_locusTag_aaSeq_dict(input1_tmp_dir, gbk_file):

    tempGenome_locusTag_aaSeq_dict = {}

    seq_record = SeqIO.read(join(input1_tmp_dir, gbk_file), 'genbank')

    for feature in seq_record.features:
        if feature.type == 'CDS':
            locusTag = feature.qualifiers['locus_tag'][0]

            if feature.qualifiers.get('translation'):
                translation = feature.qualifiers.get('translation')[0]
                tempGenome_locusTag_aaSeq_dict[locusTag] = translation

    return tempGenome_locusTag_aaSeq_dict


def get_tempModel_exrxnid_flux_dict(model):
    tempModel_exrxnid_flux_dict = {}

    model.optimize()

    # NOTE: 'f' and 'x_dict' are deprecated properties in cobra>=0.6.1.
    # TODO: This function should be upgraded upon use of cobra>=0.6.1.
    tempModel_exrxnid_flux_dict['EX_pi_e'] = model.solution.x_dict['EX_pi_e']
    tempModel_exrxnid_flux_dict['EX_co2_e'] = model.solution.x_dict['EX_co2_e']
    tempModel_exrxnid_flux_dict['EX_glc__D_e'] = model.solution.x_dict['EX_glc__D_e']
    tempModel_exrxnid_flux_dict['EX_nh4_e'] = model.solution.x_dict['EX_nh4_e']
    tempModel_exrxnid_flux_dict['EX_h2o_e'] = model.solution.x_dict['EX_h2o_e']
    tempModel_exrxnid_flux_dict['EX_h_e'] = model.solution.x_dict['EX_h_e']
    tempModel_exrxnid_flux_dict['EX_o2_e'] = model.solution.x_dict['EX_o2_e']
    tempModel_exrxnid_flux_dict[str(model.objective.keys()[0])] = model.solution.f

    return tempModel_exrxnid_flux_dict

# NOTE: RJY combined two lines of 'pyparsing.oneOf' into one
#    booleanop = pyparsing.oneOf('AND and OR or')
#    expr = pyparsing.infixNotation(gpr_regex,
#                                [
#                                (booleanop, 2, pyparsing.opAssoc.LEFT)
#                                ])
def get_gpr_fromString_toList(gpr):
    # Some locus tags contain underscores: Pseudomonas putida KT2440
    gpr_regex = pyparsing.Word(pyparsing.alphanums + '_')
    and_booleanop = pyparsing.oneOf('AND and')
    or_booleanop = pyparsing.oneOf('OR or')
    expr = pyparsing.infixNotation(gpr_regex,
                                [
                                (and_booleanop, 2, pyparsing.opAssoc.LEFT),
                                (or_booleanop, 2, pyparsing.opAssoc.LEFT)
                                ])
    gpr_list = expr.parseString(gpr)[0].asList()

    return gpr_list


def get_tempModel_biggRxnid_locusTag_dict(model):
    tempModel_biggRxnid_locusTag_dict = {}

    for j in range(len(model.reactions)):
        rxn = model.reactions[j]
        logging.debug('%s; %s', rxn.id, rxn.gene_reaction_rule)

        if rxn.gene_reaction_rule and \
                ('and' in rxn.gene_reaction_rule or 'AND' in rxn.gene_reaction_rule \
                or 'or' in rxn.gene_reaction_rule or 'OR' in rxn.gene_reaction_rule):
            gene_list = get_gpr_fromString_toList(rxn.gene_reaction_rule)
            tempModel_biggRxnid_locusTag_dict[rxn.id] = gene_list
            logging.debug('%s; %s', rxn.id, rxn.gene_reaction_rule)
        elif rxn.gene_reaction_rule:
            tempModel_biggRxnid_locusTag_dict[rxn.id] = [rxn.gene_reaction_rule]
            logging.debug('%s; %s', rxn.id, rxn.gene_reaction_rule)
        else:
            logging.debug('%s; %s', rxn.id, rxn.gene_reaction_rule)
            continue

    return tempModel_biggRxnid_locusTag_dict


def get_tempModel_locusTag_aaSeq_dict(model, tempGenome_locusTag_aaSeq_dict, options):
    tempModel_locusTag_aaSeq_dict = {}

    for g in range(len(model.genes)):
        gene = model.genes[g]

        if gene.id in tempGenome_locusTag_aaSeq_dict:
            tempModel_locusTag_aaSeq_dict[gene] = tempGenome_locusTag_aaSeq_dict[gene.id]
        else:
            logging.warning('Sequence of following gene NOT available: %s', gene.id)

    if options.model:
        logging.debug('Number of genes in %s (total %s genes): %s',
                options.model, len(model.genes), len(tempModel_locusTag_aaSeq_dict))
    elif options.acc_number:
        logging.debug('Number of genes in the model for %s (total %s genes): %s',
                options.acc_number, len(model.genes), len(tempModel_locusTag_aaSeq_dict))

    return tempModel_locusTag_aaSeq_dict


def generate_output_files(
        input1_dir,
        input1_tmp_dir,
        model,
        tempGenome_locusTag_aaSeq_dict,
        tempModel_biggRxnid_locusTag_dict,
        tempModel_locusTag_aaSeq_dict):

    # Text and FASTA files in tmp folder
    with open(join(input1_tmp_dir, 'tempModel_exrxnid_flux_dict.txt'), 'w') as f:
        for k, v in tempModel_exrxnid_flux_dict.iteritems():
            print >>f, '%s\t%s' %(k, v)

    with open(join(input1_tmp_dir, 'tempGenome_locusTag_aaSeq_dict.txt'), 'w') as f:
        for k, v in tempGenome_locusTag_aaSeq_dict.iteritems():
            print >>f, '%s\t%s' %(k, v)

    with open(join(input1_tmp_dir, 'tempModel_biggRxnid_locusTag_dict.txt'), 'w') as f:
        for k, v in tempModel_biggRxnid_locusTag_dict.iteritems():
            print >>f, '%s\t%s' %(k, v)

    with open(join(input1_dir, 'tempModel_locusTag_aaSeq.fa'), 'w') as f:
        for k, v in tempModel_locusTag_aaSeq_dict.iteritems():
            print >>f, '>%s\n%s' %(k, v)

    # Pickles in `input1` data folder
    with open(join(input1_dir, 'model.p'), 'wb') as f:
        pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input1_dir, 'tempModel_exrxnid_flux_dict.p'), 'wb') as f:
        pickle.dump(tempModel_exrxnid_flux_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input1_dir, 'tempModel_biggRxnid_locusTag_dict.p'), 'wb') as f:
        pickle.dump(tempModel_biggRxnid_locusTag_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


def make_blastDB(input1_dir):
    db_dir = join(input1_dir, 'tempBlastDB')
    query_fasta = join(input1_dir, 'tempModel_locusTag_aaSeq.fa')

    try:
        DBprogramName = gems.utils.locate_executable('makeblastdb')
        subprocess.call(
                [DBprogramName,'-in',
                query_fasta,'-out',
                db_dir,'-dbtype',
                'prot'])
    except:
        logging.warning("Failed to locate file: 'makeblastdb'")


if __name__ == '__main__':
    import time
    import warnings

    warnings.filterwarnings("ignore")

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

    options = get_options()

    input1_dir, input1_tmp_dir = get_output_dirs(options)

    if options.model:
        model = download_model_from_biggDB(input1_tmp_dir, options)
        model_info_dict = get_model_details(options)
    elif options.acc_number:
        model, model_info_dict = prepare_nonstd_model(input1_tmp_dir, options)

    gbk_file = download_gbk_from_ncbi(input1_tmp_dir, model_info_dict)
    tempGenome_locusTag_aaSeq_dict = \
            get_tempGenome_locusTag_aaSeq_dict(input1_tmp_dir, gbk_file)
    tempModel_exrxnid_flux_dict = get_tempModel_exrxnid_flux_dict(model)
    tempModel_biggRxnid_locusTag_dict = get_tempModel_biggRxnid_locusTag_dict(model)
    tempModel_locusTag_aaSeq_dict = \
        get_tempModel_locusTag_aaSeq_dict(model, tempGenome_locusTag_aaSeq_dict, options)

    generate_output_files(
            input1_dir,
            input1_tmp_dir,
            model,
            tempGenome_locusTag_aaSeq_dict,
            tempModel_biggRxnid_locusTag_dict,
            tempModel_locusTag_aaSeq_dict
            )
    make_blastDB(input1_dir)

    logging.info("Make sure to update template model options in 'run_gems.py'!")
    logging.info("Input files have been created in '/gems/io/data/input1'")
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
