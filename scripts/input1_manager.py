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
import zipfile
from Bio import Entrez, SeqIO
from cobra.util.solver import linear_reaction_coefficients
from input2_manager import ParseMNXref
from os.path import join, abspath, dirname

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import gmsm
import gmsm.io.io_utils as io_utils
from gmsm.config import load_config

def get_options():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    group = parser.add_argument_group(
        "Provide:",
        "1) BiGG Models ID (e.g., iAF1260) only or "
        "2) Genbank accession number (e.g., NC_003888.3) with a relevant SBML file or"
        "3) Genome file (GenBank or FASTA) with a relevant SBML file."
        "For the options 2) and 3), place files in 'scripts/input1_data/[organism-specific folder'].")

    group.add_argument('-m', '--model',
            dest='bigg',
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
        input1_dir = join(os.pardir, 'gmsm', 'io', 'data', 'input1', options.folder)
        if not os.path.isdir(input1_dir):
            os.makedirs(input1_dir)

        input1_tmp_dir = join(dirname(abspath(__file__)), 'input1_data', options.folder)
        if not os.path.isdir(input1_tmp_dir):
            os.makedirs(input1_tmp_dir)

    return input1_dir, input1_tmp_dir


def download_model_from_biggDB(input1_tmp_dir, options):
    model_file = ''.join([options.bigg, '.xml'])
    url = ''.join(['http://bigg.ucsd.edu/static/models/', model_file])
    logging.debug('URL for downloading a model from the BiGG Models:')
    logging.debug(url)

    model = urllib2.urlopen(url).read()

    with open(join(input1_tmp_dir, model_file), 'wb') as f:
        f.write(model)

    model = cobra.io.read_sbml_model(join(input1_tmp_dir, model_file))
    model = gmsm.utils.stabilize_model(
            model, input1_tmp_dir, '%s_2' %options.bigg, diff_name=True)

    if len(model.reactions) > 1:
        logging.debug('%s downloaded successfully', options.bigg)
    else:
        logging.debug('%s NOT downloaded successfully', options.bigg)

    return model


def get_model_details(options):
    model_info_dict = {}

    url = ''.join(['http://bigg.ucsd.edu/api/v2/models/', options.bigg])
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

    logging.debug('%s details:', options.bigg)
    logging.debug('model_bigg_id: %s', model_info_dict['model_bigg_id'])
    logging.debug('organism: %s', model_info_dict['organism'])
    logging.debug('genome_name: %s', model_info_dict['genome_name'])
    logging.debug('gene_count: %s', model_info_dict['gene_count'])

    return model_info_dict


# NOTE: Issue #333 was broken again in cobra >= 0.6.2:
#https://github.com/opencobra/cobrapy/issues/333
#Following fixation should be manually incorporated before running this code:
#https://github.com/opencobra/cobrapy/commit/ac2f2e8bd31c982e87d5f385f7a8eea35e7ba811
def get_nonstd_model(input1_tmp_dir, options):
    sbml_list = glob.glob(join(input1_tmp_dir, '*.xml'))
    logging.debug('Model found: %s', sbml_list)

    # NOTE: cobra 0.6.2 causes an error in reloading the model that is \
    #initially loaded with 'read_legacy_sbml', which considers 'fix_legacy_id'
    #model = cobra.io.read_legacy_sbml(join(input1_tmp_dir, sbml_list[0]))
    #model = gmsm.utils.stabilize_model(model, input1_tmp_dir, '')

    # NOTE: Reading 'iMK1208Edited4.xml' via 'read_sbml_model' was slightly
    #more accurate than 'read_legacy_sbml'
    #Metabolites with consistent BiGG IDs: 1045 vs 1037 (out of 1364)
    #Reactions with consistent BiGG IDs: 1231 vs 1232 (out of 1768)
    model = cobra.io.read_sbml_model(join(input1_tmp_dir, sbml_list[0]))

    return model


def get_bigg_old_new_dict_for_nonstd_model():
    mnx_parser = ParseMNXref()
    bigg_old_new_dict = mnx_parser.fix_legacy_id_using_BiGGModels()
    return bigg_old_new_dict


def fix_nonstd_model(bigg_old_new_dict, input1_tmp_dir, model, options):
    tempModel_exrxnid_flux_dict = get_tempModel_exrxnid_flux_dict(
                                                                model, bigg_old_new_dict)

    for j in range(len(model.reactions)):
        rxn = model.reactions[j]
        old_id = rxn.id

        # NOTE:
        #See '\BiGG\170405\nar-02327-data-e-2015-File017_All modifications_KHU_v2'
        #This is for the 'sco' template model.
        if rxn.id in bigg_old_new_dict and \
                rxn.id != 'FACOAL80' and \
                rxn.id != 'FE3abc':

            new_id = bigg_old_new_dict[old_id]

            # Otherwise a duplicate reaction can be inserted: e.g., ME1 -> ME2
            if new_id not in model.reactions:

                logging.debug('Reaction: %s -> %s', old_id, new_id)
                rxn.id = new_id

                model = gmsm.utils.stabilize_model(model, input1_tmp_dir, '')
                result = check_model_fluxes(
                            model, tempModel_exrxnid_flux_dict, bigg_old_new_dict, options)
                if result == 'fluxAffected':
                    rxn = model.reactions.get_by_id(new_id)
                    rxn.id = old_id
                    model = gmsm.utils.stabilize_model(model, input1_tmp_dir, '')
                    logging.debug('Reaction: %s -> %s canceled', old_id, new_id)

        if rxn.id == 'THRPS':
            logging.debug('Reaction: THRPS -> LTHRK')
            rxn.id = 'LTHRK'

            model = gmsm.utils.stabilize_model(model, input1_tmp_dir, '')
            result = check_model_fluxes(
                        model, tempModel_exrxnid_flux_dict, bigg_old_new_dict, options)
            if result == 'fluxAffected':
                rxn = model.reactions.get_by_id('LTHRK')
                rxn.id = old_id
                model = gmsm.utils.stabilize_model(model, input1_tmp_dir, '')
                logging.debug('Reaction: THRPS -> LTHRK canceled')

    tempModel_exrxnid_flux_dict = get_tempModel_exrxnid_flux_dict(
                                                                model, bigg_old_new_dict)

    for i in range(len(model.metabolites)):
        metab = model.metabolites[i]
        if metab.id in bigg_old_new_dict:
            old_id = metab.id

            if metab.compartment == 'c': #cytosol
                new_id = '_'.join([bigg_old_new_dict[old_id], 'c'])

            elif metab.compartment == 'e': #extra-organism
                new_id = '_'.join([bigg_old_new_dict[old_id], 'e'])

            elif metab.compartment == 'f': #flagellum
                new_id = '_'.join([bigg_old_new_dict[old_id], 'f'])

            elif metab.compartment == 'g': #golgi apparatus
                new_id = '_'.join([bigg_old_new_dict[old_id], 'g'])

            elif metab.compartment == 'h': #chloroplast
                new_id = '_'.join([bigg_old_new_dict[old_id], 'h'])

            elif metab.compartment == 'm': #mitochondria
                new_id = '_'.join([bigg_old_new_dict[old_id], 'm'])

            elif metab.compartment == 'n': #nucleus
                new_id = '_'.join([bigg_old_new_dict[old_id], 'n'])

            elif metab.compartment == 's': #eyespot
                new_id = '_'.join([bigg_old_new_dict[old_id], 's'])

            elif metab.compartment == 'u': #thylakoid lumen
                new_id = '_'.join([bigg_old_new_dict[old_id], 'u'])

            elif metab.compartment == 'x': #glyoxysome
                new_id = '_'.join([bigg_old_new_dict[old_id], 'x'])

            logging.debug('Metabolite: %s -> %s', old_id, new_id)
            metab.id = new_id

            model = gmsm.utils.stabilize_model(model, input1_tmp_dir, '')
            result = check_model_fluxes(
                        model, tempModel_exrxnid_flux_dict, bigg_old_new_dict, options)
            if result == 'fluxAffected':
                metab = model.metabolites.get_by_id(new_id)
                metab.id = old_id
                model = gmsm.utils.stabilize_model(model, input1_tmp_dir, '')
                logging.debug('Metabolite: %s -> %s canceled', old_id, new_id)

    model_info_dict = {}

    if options.acc_number:
        model_info_dict['genome_name'] = options.acc_number

    return model, model_info_dict


def check_model_fluxes(model, tempModel_exrxnid_flux_dict, bigg_old_new_dict, options):
    model.optimize()

    for rxnid in tempModel_exrxnid_flux_dict:
        if rxnid in model.reactions:
            rxn = model.reactions.get_by_id(rxnid)
        if rxnid not in model.reactions:
            rxn = model.reactions.get_by_id(bigg_old_new_dict[rxnid])

        tempModel_exrxnid_flux = tempModel_exrxnid_flux_dict[rxnid]
        nonzero = float(options.cobrapy.non_zero_flux_cutoff)

        if rxn.flux <= (tempModel_exrxnid_flux + nonzero) or \
                rxn.flux >= (tempModel_exrxnid_flux - nonzero):
            return ''
        else:
            logging.debug("Flux affected for %s: %s vs %s",
                          rxnid,
                          rxn.flux,
                          tempModel_exrxnid_flux_dict[rxnid])
            return 'fluxAffected'


def download_gbk_from_ncbi(input1_tmp_dir, model_info_dict):
    gbk_file = ''.join([model_info_dict['genome_name'], '.gb'])
    Entrez.email = "ehukim@kaist.ac.kr"

    handle = Entrez.efetch(db='nucleotide',
            id=model_info_dict['genome_name'], rettype='gbwithparts', retmode='text')

    seq_record = handle.read()

    with open(join(input1_tmp_dir, gbk_file), 'wb') as f:
        f.write(seq_record)

    return gbk_file


def get_tempGenome_locusTag_aaSeq_dict(input1_tmp_dir, options, **kwargs):

    tempGenome_locusTag_aaSeq_dict = {}
    options.targetGenome_locusTag_aaSeq_dict = {}
    options.targetGenome_locusTag_prod_dict = {}
    options.targetGenome_locusTag_ec_dict = {}
    options.total_cluster = 0

    if 'gbk_file' in kwargs:
        gbk_file = kwargs['gbk_file']
        filetype = 'genbank'
        seq_records = list(SeqIO.parse(join(input1_tmp_dir, gbk_file), filetype))

    elif 'gbk_file' not in kwargs and options.genome:
        input_ext = os.path.splitext(options.genome)[1]

        if input_ext in ('.gbk', '.gb', '.genbank', '.gbf', '.gbff'):
            filetype = 'genbank'
        elif input_ext in ('.fa', '.fasta', '.fna', '.faa', '.fas'):
            filetype = 'fasta'

        seq_records = list(SeqIO.parse(join(input1_tmp_dir, options.genome), filetype))

    if filetype == 'genbank':
        for seq_record in seq_records:
            io_utils.get_features_from_gbk(seq_record, options)

    elif filetype == 'fasta':
        for seq_record in seq_records:
            io_utils.get_features_from_fasta(seq_record, options)

    tempGenome_locusTag_aaSeq_dict = options.targetGenome_locusTag_aaSeq_dict

    return tempGenome_locusTag_aaSeq_dict


def get_tempModel_exrxnid_flux_dict(model):
    tempModel_exrxnid_flux_dict = {}

    flux_dist = model.optimize()

    if 'EX_pi_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_pi_e'] = float(flux_dist.fluxes.EX_pi_e)
    elif 'EX_pi(e)' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_pi(e)'] = float(flus_dist.fluxes.EX_pi(e))
    elif 'EX_pi_LPAREN_e_RPAREN_' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_pi_LPAREN_e_RPAREN_'] = \
                float(flux_dist.fluxes.EX_pi_LPAREN_e_RPAREN_)
    else:
        logging.error("'EX_pi_e' not available in the model")

    if 'EX_co2_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_co2_e'] = float(flux_dist.fluxes.EX_co2_e)
    elif 'EX_co2(e)' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_co2(e)'] = float(flux_dist.fluxes.EX_co2(e))
    elif 'EX_co2_LPAREN_e_RPAREN_' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_co2_LPAREN_e_RPAREN_'] = \
                float(flux_dist.fluxes.EX_co2_LPAREN_e_RPAREN_)
    else:
        logging.error("'EX_co2_e' not available in the model")

    if 'EX_glc__D_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_glc__D_e'] = float(flux_dist.fluxes.EX_glc__D_e)
    elif 'EX_glc(e)' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_glc(e)'] = float(flux_dist.fluxes.EX_glc(e))
    elif 'EX_glc_D_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_glc_D_e'] = float(flux_dist.fluxes.EX_glc_D_e)
    elif 'EX_glc_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_glc_e'] = float(flux_dist.fluxes.EX_glc_e)
    elif 'EX_glc_LPAREN_e_RPAREN_' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_glc_LPAREN_e_RPAREN_'] = \
                float(flux_dist.fluxes.EX_glc_LPAREN_e_RPAREN_)
    else:
        logging.error("'EX_glc__D_e' not available in the model")

    if 'EX_nh4_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_nh4_e'] = float(flux_dist.fluxes.EX_nh4_e)
    elif 'EX_nh4(e)' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_nh4(e)'] = float(flux_dist.fluxes.EX_nh4(e))
    elif 'EX_nh4_LPAREN_e_RPAREN_' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_nh4_LPAREN_e_RPAREN_'] = \
                float(flux_dist.fluxes.EX_nh4_LPAREN_e_RPAREN_)
    else:
        logging.error("'EX_nh4_e' not available in the model")

    if 'EX_h2o_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_h2o_e'] = float(flux_dist.fluxes.EX_h2o_e)
    elif 'EX_h2o(e)' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_h2o(e)'] = float(flux_dist.fluxes.EX_h2o(e))
    elif 'EX_h2o_LPAREN_e_RPAREN_' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_h2o_LPAREN_e_RPAREN_'] = \
                float(flux_dist.fluxes.EX_h2o_LPAREN_e_RPAREN_)
    else:
        logging.error("'EX_h2o_e' not available in the model")

    if 'EX_h_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_h_e'] = float(flux_dist.fluxes.EX_h_e)
    elif 'EX_h(e)' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_h(e)'] = float(flux_dist.fluxes.EX_h(e))
    elif 'EX_h_LPAREN_e_RPAREN_' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_h_LPAREN_e_RPAREN_'] = \
                float(flux_dist.fluxes.EX_h_LPAREN_e_RPAREN_)
    else:
        logging.error("'EX_h_e' not available in the model")

    if 'EX_o2_e' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_o2_e'] = float(flux_dist.fluxes.EX_o2_e)
    elif 'EX_o2_LPAREN_e_RPAREN_' in model.reactions:
        tempModel_exrxnid_flux_dict['EX_o2_LPAREN_e_RPAREN_'] = \
                float(flux_dist.fluxes.EX_o2_LPAREN_e_RPAREN_)
    else:
        logging.error("'EX_o2_e' not available in the model")

    if linear_reaction_coefficients(model).keys()[0].id:
        tempModel_exrxnid_flux_dict[linear_reaction_coefficients(model).keys()[0].id] = \
                float(flux_dist.objective_value)
    else:
        logging.error("Objective function should be designated in the model")

    for rxn in tempModel_exrxnid_flux_dict:
        logging.debug("Reaction: %s; Flux: %s", rxn, tempModel_exrxnid_flux_dict[rxn])

    return tempModel_exrxnid_flux_dict


# NOTE: RJY combined two lines of 'pyparsing.oneOf' into one
#    booleanop = pyparsing.oneOf('AND and OR or')
#    expr = pyparsing.infixNotation(gpr_regex,
#                                [
#                                (booleanop, 2, pyparsing.opAssoc.LEFT)
#                                ])
def get_gpr_fromString_toList(gpr):
    # Some locus tags contain underscores: Pseudomonas putida KT2440
    # Some locus tags contain periods ("."): Chlamydomonas reinhardtii
    gpr_regex = pyparsing.Word(pyparsing.alphanums + '_' + '.')
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

    if options.bigg:
        logging.debug('Number of genes in %s (total %s genes): %s',
                options.bigg, len(model.genes), len(tempModel_locusTag_aaSeq_dict))
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
        tempModel_locusTag_aaSeq_dict,
        options):

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
        DBprogramName = gmsm.utils.locate_executable('makeblastdb')
        subprocess.call(
                [DBprogramName,'-in',
                query_fasta,'-out',
                db_dir,'-dbtype',
                'prot'])
    except:
        logging.warning("Failed to locate file: 'makeblastdb'")


def create_zip_file(input1_tmp_dir, options):
    input1_tmp_dir_list = glob.glob(join(input1_tmp_dir, '*.*'))
    zip_file = '%s_input1_data.zip' %options.folder

    for output in input1_tmp_dir_list:
        if zip_file in output:
            input1_tmp_dir_list.remove(output)

    zip = zipfile.ZipFile(join(input1_tmp_dir, '%s_input1_data.zip' %options.folder),
                            'w',
                            zipfile.ZIP_DEFLATED)

    for output in input1_tmp_dir_list:
        zip.write(output, os.path.basename(output))

    zip.close()
    return input1_tmp_dir_list


def remove_input1_tmp_dir_files(input1_tmp_dir, input1_tmp_dir_list, options):

    for output in input1_tmp_dir_list:
        os.remove(join(input1_tmp_dir, output))


if __name__ == '__main__':
    import time
    import warnings

    sys.path.insert(0, abspath(join(dirname(__file__), '..')))
    from gmsm.config import load_config

    warnings.filterwarnings("ignore")

    start = time.time()

    # Logging setup
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

    options = get_options()
    load_config(options)
    input1_dir, input1_tmp_dir = get_output_dirs(options)

    # Logfile setup
    logger = logging.getLogger('')
    fomatter = logging.Formatter('[%(levelname)s|%(filename)s:%(lineno)s] > %(message)s')
    fh = logging.FileHandler(join(input1_tmp_dir, '%s.log' %options.folder), mode = 'w')
    fh.setFormatter(fomatter)
    logger.addHandler(fh)

    if options.bigg:
        model = download_model_from_biggDB(input1_tmp_dir, options)
        model_info_dict = get_model_details(options)
    elif options.acc_number or options.genome:
        model = get_nonstd_model(input1_tmp_dir, options)
        bigg_old_new_dict = get_bigg_old_new_dict_for_nonstd_model()
        model, model_info_dict = fix_nonstd_model(
                bigg_old_new_dict,
                input1_tmp_dir,
                model,
                options)

    if options.bigg or options.acc_number:
        gbk_file = download_gbk_from_ncbi(input1_tmp_dir, model_info_dict)
        tempGenome_locusTag_aaSeq_dict = \
            get_tempGenome_locusTag_aaSeq_dict(
                    input1_tmp_dir, options, gbk_file = gbk_file)
    elif options.genome:
        tempGenome_locusTag_aaSeq_dict = \
                get_tempGenome_locusTag_aaSeq_dict(input1_tmp_dir, options)

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
            tempModel_locusTag_aaSeq_dict,
            options
            )
    make_blastDB(input1_dir)

    logging.info("Make sure to update template model options in 'run_gmsm.py' and 'output_file_manager.py!")
    logging.info("Input files have been created in '/gmsm/io/data/input1/%s'",
                    options.folder)
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))

    input1_tmp_dir_list = create_zip_file(input1_tmp_dir, options)
    remove_input1_tmp_dir_files(input1_tmp_dir, input1_tmp_dir_list, options)
