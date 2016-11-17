'''
2014-2015
Hyun Uk Kim, Tilmann Weber, Kyu-Sang Hwang and Jae Yong Ryu
'''

#Wildcard imports should never be used in production code.
import argparse
import logging
import multiprocessing
import os
import pickle
import sys
import time

from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
from argparse import Namespace

from prunPhase import (
    get_temp_fasta,
    get_target_gbk,
    get_targetGenomeInfo,
    get_target_fasta,
    make_blastDB,
    run_blastp,
    parseBlaspResults,
    makeBestHits_dict,
    getBBH,
    get_nonBBH,
    labelRxnToRemove,
    pruneModel,
    swap_locusTag_tempModel
)
from augPhase import (
    get_targetGenome_locusTag_ec_nonBBH_dict,
    make_all_rxnInfo_fromRefSeq,
    get_mnxr_list_from_modelPrunedGPR,
    check_existing_rxns,
    get_mnxr_using_kegg,
    extract_rxn_mnxm_coeff,
    add_nonBBH_rxn
)


def main():
    start = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--model', dest='orgName', default='sco', choices=['eco','sco'], help="Specify a template model for the target modeling")
    parser.add_argument('-s', '--smr', dest='smr_generation', default=False, choices=[True,False], help="Specify whether to run secondary metabolic modeling")
    parser.add_argument('-i', '--input', dest='input', default='input', help="Specify input directory")
    parser.add_argument('-o', '--output', dest='output', default='output', help="Specify output directory")
    parser.add_argument('-e', '--ec', dest='eficaz', action='store_true', default=False, help="Run EC number prediction using EFICAz")
    parser.add_argument('-c', '--cpu', dest='cpus', default=multiprocessing.cpu_count(), type=int, help="How many CPUs to use in parallel. (default: %(default)s)")
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help="Print debugging information to stderr")

    options = parser.parse_args()

    if options.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    #Create output folders
    folders = ['0_EFICAz_results', '1_blastp_results', '2_primary_metabolic_model', '3_temp_models', '4_complete_model']

    for folder in folders:
        if not os.path.isdir(options.output+'/'+folder):
            os.makedirs(options.output+'/'+folder)

    options.outputfoldername = options.output+'/'+folders[0]

    get_genome_files(options)

    get_homolgs(options)

    model = get_pickles_prunPhase(options)

    modelPrunedGPR = run_prunPhase(model, options)

    if options.targetGenome_locusTag_ec_dict:
        get_pickles_augPhase(options)
        run_augPhase(modelPrunedGPR, options)
    else:
        logging.debug("No EC_number found in the submitted gbk file")

    #Secondary metabolic modeling
    if options.smr_generation:
        run_smr_generation

    generate_outputs

    logging.debug(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


def get_genome_files(options):
    logging.debug("Reading input genome filesg..")

    logging.debug("Looking for a gbk file of a template model genome..")
    get_temp_fasta(options)

    logging.debug("Looking for a gbk file of a target genome..")
    get_target_gbk(options)

    logging.debug("Reading genbank file of the target genome..")
    get_targetGenomeInfo(options, 'genbank')

    logging.debug("Looking for a fasta file of a target genome..")
    get_target_fasta(options)


def get_homolgs(options):
    logging.debug("Searching bidirectional homolgs..")
    logging.debug("Generating a DB for the genes from the target genome..")
    make_blastDB(options)

    logging.debug("Running BLASTP #1: genes in the target genome against genes in the template model..")
    run_blastp(target_fasta='./%s/1_blastp_results/targetGenome_locusTag_aaSeq.fa' %options.output, blastp_result='./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, db_dir = '%s/tempBlastDB' %(options.input1), evalue=1e-30)

    logging.debug("Running BLASTP #2: genes in the template model against genes in the target genome..")
    run_blastp(target_fasta=options.temp_fasta, blastp_result='./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, db_dir = './%s/1_blastp_results/targetBlastDB' %options.output, evalue=1e-30)

    logging.debug("Parsing the results of BLASTP #1..")
    blastpResults_dict1 = parseBlaspResults('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, './%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)
 
    logging.debug("Parsing the results of BLASTP #2..")
    blastpResults_dict2 = parseBlaspResults('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, './%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

    logging.debug("Selecting the best hits for BLASTP #1..")
    bestHits_dict1 = makeBestHits_dict('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)

    logging.debug("Selecting the best hits for BLASTP #2..")
    bestHits_dict2 = makeBestHits_dict('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

    logging.debug("Selecting the bidirectional best hits..")
    getBBH(bestHits_dict1, bestHits_dict2, options)

    logging.debug("Selecting genes that are not bidirectional best hits..")
    get_nonBBH(options)


#For model pruning phase
#Only model file is not saved in Namespace
def get_pickles_prunPhase(options):
    logging.debug("Loading pickle files of the parsed template model and its relevant genbank data..")
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    tempModel_biggRxnid_locusTag_dict = pickle.load(open('%s/tempModel_biggRxnid_locusTag_dict.p' %(options.input1),'rb'))
    options.tempModel_biggRxnid_locusTag_dict = tempModel_biggRxnid_locusTag_dict

    return model


def run_prunPhase(model, options):
    logging.debug("Pruning phase starting..")
    logging.debug("Labeling reactions with nonhomologous genes to remove from the template model..")
    rxnToRemove_dict = labelRxnToRemove(model, options)

    logging.debug("Removing reactions with nonhomologous genes from the template model..")
    modelPruned, rxnToRemoveEssn_dict, rxnRemoved_dict, rxnRetained_dict = pruneModel(model, options, 'gurobi')

    logging.debug("Correcting GPR associations in the template model..")
    modelPrunedGPR = swap_locusTag_tempModel(modelPruned, options)

    return modelPrunedGPR


#For model augmentation  phase
def get_pickles_augPhase(options):
    logging.debug("Loading pickle files necessary for the model augmentation phase..")

    bigg_mnxr_dict = pickle.load(open('./input2/bigg_mnxr_dict.p','rb'))
    options.bigg_mnxr_dict = bigg_mnxr_dict
    kegg_mnxr_dict = pickle.load(open('./input2/kegg_mnxr_dict.p','rb'))
    options.kegg_mnxr_dict = kegg_mnxr_dict
    mnxr_kegg_dict = pickle.load(open('./input2/mnxr_kegg_dict.p','rb'))
    options.mnxr_kegg_dict = mnxr_kegg_dict
    mnxr_rxn_dict = pickle.load(open('./input2/mnxr_rxn_dict.p','rb'))
    options.mnxr_rxn_dict = mnxr_rxn_dict

    bigg_mnxm_compound_dict = pickle.load(open('./input2/bigg_mnxm_compound_dict.p','rb'))
    options.bigg_mnxm_compound_dict = bigg_mnxm_compound_dict
    mnxm_bigg_compound_dict = pickle.load(open('./input2/mnxm_bigg_compound_dict.p','rb'))
    options.mnxm_bigg_compound_dict = mnxm_bigg_compound_dict
    kegg_mnxm_compound_dict = pickle.load(open('./input2/kegg_mnxm_compound_dict.p','rb'))
    options.kegg_mnxm_compound_dict = kegg_mnxm_compound_dict
    mnxm_kegg_compound_dict = pickle.load( open('./input2/mnxm_kegg_compound_dict.p','rb'))
    options.mnxm_kegg_compound_dict = mnxm_kegg_compound_dict

    mnxm_compoundInfo_dict = pickle.load(open('./input2/mnxm_compoundInfo_dict.p','rb'))
    options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(options.input1),'rb'))
    options.template_exrxnid_flux_dict = template_exrxnid_flux_dict


def run_augPhase(modelPrunedGPR, options):
    logging.debug("Augmentation phase starting..")
    logging.debug("Creating various dictionary files for the nonBBH gene-associted reactions...")

    get_targetGenome_locusTag_ec_nonBBH_dict(options)

    #Two nested functions
    #def get_rxnid_from_ECNumber(enzymeEC):
    #def get_rxnInfo_from_rxnid(rxnid):
    make_all_rxnInfo_fromRefSeq(options)

    get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, options) 

    logging.debug("Adding the nonBBH gene-associated reactions...")
    check_existing_rxns(options)

    get_mnxr_using_kegg(options)

    extract_rxn_mnxm_coeff(options)

    target_model = add_nonBBH_rxn(modelPrunedGPR, options)


def generate_outputs():
    #Output files
    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    #Cobrapy IO module seems to have an error for adding new reactions
    write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(options.output, options.output, options.orgName))
    target_model = create_cobra_model_from_sbml_file('./%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(options.output, options.output, options.orgName))
    write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(options.output, options.output, options.orgName))

    #Output on screen
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    print "Number of genes:", len(model.genes), "/", len(modelPruned.genes), "/", len(target_model.genes)
    print "Number of reactions:", len(model.reactions), "/", len(modelPruned.reactions), "/", len(target_model.reactions)
    print "Number of metabolites:",  len(model.metabolites), "/", len(modelPruned.metabolites), "/", len(target_model.metabolites)

    fp1 = open('./%s/2_primary_metabolic_model/%s_target_model_reactions.txt' %(options.output, options.output), "w")
    fp2 = open('./%s/2_primary_metabolic_model/%s_target_model_metabolites.txt' %(options.output, options.output), "w")
    fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Lower bound"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
    fp2.write("Metabolite ID"+"\t"+"Metabolite name"+"\t"+"Formula"+"\t"+"Compartment"+"\n")

    print "\n"
    for j in range(len(target_model.reactions)):
        rxn = target_model.reactions[j]
        print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.lower_bound, rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)

    for i in range(len(target_model.metabolites)):
        metab = target_model.metabolites[i]
        print >>fp2, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula, metab.compartment)

    fp1.close()
    fp2.close()


if __name__ == '__main__':
    main()
