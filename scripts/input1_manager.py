# -*- coding: utf-8 -*-

import argparse
import ast
import cobra
import os
import pickle
import re
import subprocess
import sys
import urllib2
from Bio import Entrez, SeqIO
from os.path import join, abspath, dirname

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

    if len(model.reactions) > 1:
        logging.debug('%s downloaded successfully', options.model)
    else:
        logging.debug('%s NOT downloaded successfully', options.model)


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


#Looks for .xml and .gb(k) files in the pre-defined folder
#def get_tempInfo(orgName):
#    for root, dirs, files in os.walk('./input1/%s/' %(orgName)):
#        for file in files:
#            if file.endswith('.gb') or file.endswith('.gbk'):
#		tempGenome = os.path.join(root, file)
#            if file.endswith('.xml'):
#		tempModel = os.path.join(root, file)
#	if tempGenome and tempModel:
#            return root, tempGenome, tempModel
#	else:
#            sys.exit(1)

#Files contents in SBML formatted into Dictionary {key:value} type.
#Created by Jae Yong Ryu, on 2014. 1. 13.
#Updated by Hyun Uk Kim, Feb and May 2014
def readCobraModel(root, sbmlFile):
    sbml_level = 2
    cobra_model = create_cobra_model_from_sbml_file(sbmlFile, print_time=True)
    cobra_reaction_dic={}
    cobra_metabolite_dic={}

#Double dictionary:
    for each_reaction in cobra_model.reactions:
        cobra_reaction_dic[each_reaction.id]={'id':each_reaction.id, 'name':each_reaction.name,'reactants':each_reaction.get_reactants(), 'products':each_reaction.get_products(), 'reactants_coff':each_reaction.get_coefficients(each_reaction.get_reactants()), 'products_coff':each_reaction.get_coefficients(each_reaction.get_products()), 'gpr':each_reaction.gene_reaction_rule, 'lowerbound':each_reaction.lower_bound, 'upperbound':each_reaction.upper_bound, 'subsystem':each_reaction.subsystem, 'reversibility':each_reaction.reversibility}

    for each_metabolite in cobra_model.metabolites:
        abbrmetabolite = each_metabolite.id[:-2]+'[%s]'%(each_metabolite[-1])
        if len(str(each_metabolite.formula).strip()) == 0:
            formula = ''
        else:
            formula = each_metabolite.formula
        cobra_metabolite_dic[each_metabolite.id]={'id':each_metabolite.id, 'name':each_metabolite.name, 'abbrname':abbrmetabolite, 'reaction':each_metabolite.get_reaction(),'compartment':each_metabolite.compartment, 'formula':formula}

    pickle.dump(cobra_model, open('%s/model.p' %(root),'wb'))

    return cobra_model, cobra_reaction_dic


#This file reads the FASTA sequences and saves in Dictionary type.
#Regarding SeqIO code used here, check JaeYong's email on Jan 18, 2014.
#Updated Nov, 014
def readSeq_templateGenome(gbkFile, FileType):
    fp = open('./forChecking/tempGenome_locusTag_aaSeq.fa','w')

    tempGenome_locusTag_aaSeq_dict = {}

#Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)

    for feature in record.features:
        if feature.type == 'CDS':

#Note that the numbers of CDS and "translation" do not match.
#There are occasions that CDS does not have "translation".
            if 'translation' in feature.qualifiers:

#Retrieving "locus_tag (i.e., ORF name)" for each CDS
                feature_locus_tag = str(feature.qualifiers['locus_tag']).split("']")
                feature_locus_tag = str(feature_locus_tag[0]).split("['")
#                fp.write(str(feature_locus_tag[1]))
#                fp.write("\t")

#Retrieving "translation (i.e., amino acid sequences)" for each CDS
                feature_translation = str(feature.qualifiers.get('translation')).split("']")
                feature_translation = str(feature_translation[0]).split("['")

#Saves in Dictionary
                tempGenome_locusTag_aaSeq_dict[feature_locus_tag[1]] = feature_translation[1]
                print >>fp, '>%s\n%s' % (str(feature_locus_tag[1]), str(feature_translation[1]))

    fp.close()
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


#This file reads the model and amino acid sequences, and parses them.
#Created by Hyun Uk Kim, Jan 2014
def parse_templateModel_gpr(outputFile1, root, cobra_reaction_dic, tempGenome_locusTag_aaSeq_dict):
    fp1 = open(outputFile1, 'w')
    tempModel_biggRxnid_locusTag_dict = {}
    tempModel_locusTag_aaSeq_dict = {} #For reactions with genes and AA seq

    #Starts model parsing
    for Reaction_name in cobra_reaction_dic.keys():
        gpr = cobra_reaction_dic[Reaction_name]['gpr']

        #Parsing gene associations into a list
        Genes_list = get_gpr_fromString_toList(gpr)

        #Saves reaction names and their intact GPR associations
        #Input for "makeBooleanFormat" in "Reconstruction_GPRCalculation_Edited.py"
        tempModel_biggRxnid_locusTag_dict[Reaction_name] = Genes_list

        for Genes_each_list in Genes_list:
            #Checks if the element itself is List.
            #'if type(Genes_each_list) == list' also works.
            #Genes connected with 'AND' come in the List type.
            if isinstance(Genes_each_list, list):
                for Genes_each2_list in Genes_each_list:
                    if Genes_each2_list in tempGenome_locusTag_aaSeq_dict.keys():
                        #Dictionary - ORF:AA seq
                        #Stores only reactions with genes and their amino acid seq in Dictionary
                        if tempGenome_locusTag_aaSeq_dict[Genes_each2_list]:
                            tempModel_locusTag_aaSeq_dict[Genes_each2_list] = tempGenome_locusTag_aaSeq_dict[Genes_each2_list]

            #Single genes for a reaction, or genes connected with 'OR'
            else:
                if Genes_each_list in tempGenome_locusTag_aaSeq_dict.keys():

                    #Dictionary - ORF:AA seq
                    #Stores only reactions with genes and their amino acid seq in Dictionary
                    if tempGenome_locusTag_aaSeq_dict[Genes_each_list]:
                        tempModel_locusTag_aaSeq_dict[Genes_each_list] = tempGenome_locusTag_aaSeq_dict[Genes_each_list]

    pickle.dump(tempModel_biggRxnid_locusTag_dict, open('%s/tempModel_biggRxnid_locusTag_dict.p' %(root),'wb'))

    for key in tempModel_biggRxnid_locusTag_dict.keys():
	print >>fp1, '%s\t%s' %(key, tempModel_biggRxnid_locusTag_dict[key])
    fp1.close()
    return tempModel_biggRxnid_locusTag_dict, tempModel_locusTag_aaSeq_dict


#This function was removed, but restored
#"tempModel_locusTag_aaSeq_dict" contains a unique key, which is locusTag
#Creating a fasta file directly using the above function requires additional filters
def makeFasta_locusTag_aaSeq(root, tempModel_locusTag_aaSeq_dict):
    fp1 = open('%s/tempModel_locusTag_aaSeq.fa' %(root),'w')

    for locusTag in tempModel_locusTag_aaSeq_dict.keys():
        aaSeq = tempModel_locusTag_aaSeq_dict[locusTag]
        print >>fp1, '>%s\n%s' % (str(locusTag), str(aaSeq))

    fp1.close()


#making database files using fasta files
def make_blastDB(root, query_fasta):
    db_dir = '%s/tempBlastDB' %(root)
    DBprogramName = './blastpfiles/makeblastdb.exe'
    print DBprogramName, query_fasta, db_dir
    subprocess.call([DBprogramName,'-in',query_fasta,'-out',db_dir,'-dbtype','prot'])
    return


#Reaction IDs and their flux values for major Exchange reactions
def pickle_template_exchange_rxnid_flux_dict():
    fp1 = open('./input1/eco/tempModel_exrxnid_flux_dict.txt',"r")
    template_exrxnid_flux_dict = {}
    text = fp1.readline()

    while text:
	text = text.split('\t')
	text[0] = text[0].strip()
	text[1] = text[1].strip()
	template_exrxnid_flux_dict[text[0]] = text[1]

	print text[0],template_exrxnid_flux_dict[text[0]]
	text = fp1.readline()

    print template_exrxnid_flux_dict
    pickle.dump(template_exrxnid_flux_dict, open('./input1/eco/tempModel_exrxnid_flux_dict.p','wb'))
    fp1.close()

if __name__ == '__main__':
    import logging
    import time

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

    options = get_bigg_model_id()
    download_model_from_biggDB(options)
    model_info_dict = get_model_details(options)
    download_gbk_from_ncbi(model_info_dict)

    logging.info("Make sure to update template model options in 'run_gems.py'!")
    logging.info("Input files have been created in '/gems/io/data/input1'")
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
