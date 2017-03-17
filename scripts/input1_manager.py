# -*- coding: utf-8 -*-

import argparse
import ast
import cobra
import os
import pickle
import subprocess
import sys
import urllib2
from Bio import Entrez, SeqIO
from os.path import join, abspath, dirname

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import gems

input1_dir = join(os.pardir, 'gems', 'io', 'data', 'input1')
input1_tmp_dir = join(dirname(abspath(__file__)), 'input1_data')

def get_bigg_model_id():
    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--model',
            dest='model',
            help = "Specify BiGG ID of a metabolic model to prepare as a template model")

    options = parser.parse_args()
    logging.debug(options)

    return options


def download_model_from_biggDB(options):
    model_file = ''.join([options.model, '.xml'])
    url = ''.join(['http://bigg.ucsd.edu/static/models/', model_file])
    logging.debug('URL for downloading a model from the BiGG Models:')
    logging.debug(url)

    model = urllib2.urlopen(url).read()

    with open(join(input1_tmp_dir, model_file), 'wb') as f:
        f.write(model)

    model = cobra.io.read_sbml_model(join(input1_tmp_dir, model_file))
    cobra.io.write_sbml_model(model,
            join(input1_tmp_dir, model_file), use_fbc_package=False)
    model = cobra.io.read_sbml_model(join(input1_tmp_dir, model_file))

    if len(model.reactions) > 1:
        logging.debug('%s downloaded successfully', options.model)
    else:
        logging.debug('%s NOT downloaded successfully', options.model)

    return model


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


def get_model_details(options):
    model_info_dict = {}

    url = ''.join(['http://bigg.ucsd.edu/api/v2/models/', options.model])
    logging.debug('URL for accessing model details the BiGG Models:')
    logging.debug(url)

    model_info = urllib2.urlopen(url).read()

    model_info_dict = ast.literal_eval(model_info)
    logging.debug('%s details:', options.model)
    logging.debug('model_bigg_id: %s', model_info_dict['model_bigg_id'])
    logging.debug('organism: %s', model_info_dict['organism'])
    logging.debug('genome_name: %s', model_info_dict['genome_name'])
    logging.debug('gene_count: %s', model_info_dict['gene_count'])

    return model_info_dict


def download_gbk_from_ncbi(model_info_dict):
    gbk_file = ''.join([model_info_dict['genome_name'], '.gb'])
    Entrez.email = "ehukim@kaist.ac.kr"

    handle = Entrez.efetch(db='nucleotide',
            id=model_info_dict['genome_name'], rettype='gbwithparts', retmode='text')

    seq_record = handle.read()

    with open(join(input1_tmp_dir, gbk_file), 'wb') as f:
        f.write(seq_record)

    return gbk_file


def get_tempGenome_locusTag_aaSeq_dict(gbk_file):

    tempGenome_locusTag_aaSeq_dict = {}

    seq_record = SeqIO.read(join(input1_tmp_dir, gbk_file), 'genbank')

    for feature in seq_record.features:
        if feature.type == 'CDS':
            locusTag = feature.qualifiers['locus_tag'][0]

            if feature.qualifiers.get('translation'):
                translation = feature.qualifiers.get('translation')[0]
                tempGenome_locusTag_aaSeq_dict[locusTag] = translation

    return tempGenome_locusTag_aaSeq_dict


def get_gpr_fromString_toList(line):
    calcNewList = []
    line = line.strip()
    calcList = line.split('or')
    for c in calcList:
        c = c.replace('(','')
        c = c.replace(')','')
        c = c.replace(' ','')
        c = c.strip()
        if 'and' in c:
            newlist = c.split('and')
            newlist = list(set(newlist))
            newlist.sort()
            calcNewList.append(newlist)
        else:
            geneid=c.strip()
            if geneid not in calcNewList:
                calcNewList.append(geneid)

    return calcNewList


def get_tempModel_biggRxnid_locusTag_dict(model):
    tempModel_biggRxnid_locusTag_dict = {}

    for j in range(len(model.reactions)):
        rxn = model.reactions[j]
        gene_list = get_gpr_fromString_toList(rxn.gene_reaction_rule)
        tempModel_biggRxnid_locusTag_dict[rxn.id] = gene_list

    return tempModel_biggRxnid_locusTag_dict


def get_tempModel_locusTag_aaSeq_dict(model, tempGenome_locusTag_aaSeq_dict, options):
    tempModel_locusTag_aaSeq_dict = {}

    for g in range(len(model.genes)):
        gene = model.genes[g]

        if gene.id in tempGenome_locusTag_aaSeq_dict:
            tempModel_locusTag_aaSeq_dict[gene] = tempGenome_locusTag_aaSeq_dict[gene.id]
        else:
            logging.warning('Sequence of following gene NOT available: %s', gene.id)

    logging.debug('Number of genes in %s (total %s genes): %s',
            options.model, len(model.genes), len(tempModel_locusTag_aaSeq_dict))

    return tempModel_locusTag_aaSeq_dict


def generate_output_files(model,
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

    with open(join(input1_tmp_dir, 'tempModel_locusTag_aaSeq_dict.fa'), 'w') as f:
        for k, v in tempModel_locusTag_aaSeq_dict.iteritems():
            print >>f, '>%s\n%s' %(k, v)

    # Pickles in `input1` data folder
    with open(join(input1_dir, 'model.p'), 'wb') as f:
        pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input1_dir, 'tempModel_exrxnid_flux_dict.p'), 'wb') as f:
        pickle.dump(tempModel_exrxnid_flux_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input1_dir, 'tempModel_biggRxnid_locusTag_dict.p'), 'wb') as f:
        pickle.dump(tempModel_biggRxnid_locusTag_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input1_dir, 'tempModel_locusTag_aaSeq_dict.p'), 'wb') as f:
        pickle.dump(tempModel_locusTag_aaSeq_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


def make_blastDB():
    db_dir = join(input1_dir, 'tempBlastDB')
    query_fasta = join(input1_tmp_dir, 'tempModel_locusTag_aaSeq_dict.fa')

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
    import logging
    import time

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

    options = get_bigg_model_id()
    model = download_model_from_biggDB(options)
    tempModel_exrxnid_flux_dict = get_tempModel_exrxnid_flux_dict(model)
    model_info_dict = get_model_details(options)
    gbk_file = download_gbk_from_ncbi(model_info_dict)
    tempGenome_locusTag_aaSeq_dict = get_tempGenome_locusTag_aaSeq_dict(gbk_file)
    tempModel_biggRxnid_locusTag_dict = get_tempModel_biggRxnid_locusTag_dict(model)
    tempModel_locusTag_aaSeq_dict = \
        get_tempModel_locusTag_aaSeq_dict(model, tempGenome_locusTag_aaSeq_dict, options)
    generate_output_files(
            model,
            tempGenome_locusTag_aaSeq_dict,
            tempModel_biggRxnid_locusTag_dict,
            tempModel_locusTag_aaSeq_dict
            )
    make_blastDB()

    logging.info("Make sure to update template model options in 'run_gems.py'!")
    logging.info("Input files have been created in '/gems/io/data/input1'")
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
