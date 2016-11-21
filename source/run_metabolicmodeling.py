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
from sec_met_rxn_generation import (
    get_cluster_location,
    get_cluster_info_from_seq_record,
    get_cluster_product,
    get_cluster_domain,
    get_cluster_monomers,
    get_cluster_module,
    get_currency_metabolites,
    get_total_currency_metab_coeff,
    get_all_metab_coeff,
    add_sec_met_rxn,
    check_producibility_sec_met,
    get_monomers_nonprod_sec_met,
    get_monomers_prod_sec_met
)

def main():
    start = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--model', dest='orgName', default='sco', choices=['eco','sco'], help="Specify a template model for the target modeling")
    parser.add_argument('--disable-modeling', dest='pmr_generation', default=True, action='store_false', help='Disable primary metabolic modeling')
    parser.add_argument('-s', '--smr', dest='smr_generation', default=False, action=('store_true'), help="Specify whether to run secondary metabolic modeling")
    parser.add_argument('-i', '--input', dest='input', default='input', help="Specify input directory")
    parser.add_argument('-o', '--output', dest='output', default='output', help="Specify output directory")
    parser.add_argument('-e', '--ec', dest='eficaz', action='store_true', default=False, help="Run EC number prediction using EFICAz")
    parser.add_argument('-c', '--cpu', dest='cpus', default=multiprocessing.cpu_count(), type=int, help="How many CPUs to use in parallel. (default: %(default)s)")
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help="Print debugging information to stderr")

    options = parser.parse_args()

    print options.pmr_generation
    if options.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    #Create output folders
    folders = ['0_EFICAz_results', '1_blastp_results', '2_primary_metabolic_model', '3_temp_models', '4_complete_model']

    for folder in folders:
        if not os.path.isdir(options.output+'/'+folder):
            os.makedirs(options.output+'/'+folder)

    options.outputfoldername = options.output+'/'+folders[0]

    get_genome_files(options)

    get_pickles_add_rxn(options)

    if options.pmr_generation:
        get_homolgs(options)

        model = get_pickles_prunPhase(options)

        modelPrunedGPR = run_prunPhase(model, options)

        if options.targetGenome_locusTag_ec_dict:

            target_model = run_augPhase(modelPrunedGPR, options)
        else:
            logging.debug("No EC_number found in the submitted gbk file")

        generate_outputs(model, modelPrunedGPR, target_model, options)

    #Secondary metabolic modeling
    if options.smr_generation:
        for model_file in os.listdir(options.output+'/'+'2_primary_metabolic_model'):
            if model_file.endswith('.xml'):
                target_model = create_cobra_model_from_sbml_file(options.output+'/'+'2_primary_metabolic_model/'+model_file)

                logging.debug("Generating secondary metabolite biosynthesizing reactions..")
                cluster_nr = 1
                logging.debug("Total number of clusters: %s" %options.total_cluster)
                while cluster_nr <= options.total_cluster:
                    logging.debug("Cluster number: %s" %cluster_nr)
                    sec_met_rxn_generation(cluster_nr, target_model, options)
                    cluster_nr += 1
        if not model_file:
            logging.debug("COBRA-compliant SBML file needed")

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


#For model augmentation phase in both primary and secondary modeling
def get_pickles_add_rxn(options):
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
    labelRxnToRemove(model, options)

    logging.debug("Removing reactions with nonhomologous genes from the template model..")
    modelPruned = pruneModel(model, options, 'gurobi')

    logging.debug("Correcting GPR associations in the template model..")
    modelPrunedGPR = swap_locusTag_tempModel(modelPruned, options)

    return modelPrunedGPR


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

    return target_model


def generate_outputs(model, modelPrunedGPR, target_model, options):
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
    logging.debug("Number of genes: %s; %s; %s" %(len(model.genes), len(modelPrunedGPR.genes), len(target_model.genes)))
    logging.debug("Number of reactions: %s; %s; %s" %(len(model.reactions), len(modelPrunedGPR.reactions), len(target_model.reactions)))
    logging.debug("Number of metabolites: %s; %s; %s" %(len(model.metabolites), len(modelPrunedGPR.metabolites), len(target_model.metabolites)))

    fp1 = open('./%s/2_primary_metabolic_model/%s_target_model_reactions.txt' %(options.output, options.output), "w")
    fp2 = open('./%s/2_primary_metabolic_model/%s_target_model_metabolites.txt' %(options.output, options.output), "w")
    fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Lower bound"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
    fp2.write("Metabolite ID"+"\t"+"Metabolite name"+"\t"+"Formula"+"\t"+"Compartment"+"\n")

    for j in range(len(target_model.reactions)):
        rxn = target_model.reactions[j]
        print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.lower_bound, rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)

    for i in range(len(target_model.metabolites)):
        metab = target_model.metabolites[i]
        print >>fp2, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula, metab.compartment)

    fp1.close()
    fp2.close()


def sec_met_rxn_generation(cluster_nr, target_model, options):

    prod_sec_met_dict = {}
    nonprod_sec_met_dict = {}

    get_cluster_location(cluster_nr, options)

    get_cluster_info_from_seq_record(options)

    get_cluster_product(cluster_nr, options)

    if 't1pks' in options.product or 'nrps' in options.product:
        get_cluster_domain(options)

        get_cluster_monomers(options)

        get_cluster_module(options)

        get_currency_metabolites(options)

        get_total_currency_metab_coeff(options)

        get_all_metab_coeff(options)

        target_model = add_sec_met_rxn(target_model, options)

        target_model = check_producibility_sec_met(target_model, options)

        if target_model.solution.f < 0.0001:
            nonprod_sec_met_metab_list = get_monomers_nonprod_sec_met(options)
            nonprod_sec_met_dict[options.product] = nonprod_sec_met_metab_list
        else:
            prod_sec_met_metab_list = get_monomers_prod_sec_met(options)
            prod_sec_met_dict[options.product] = prod_sec_met_metab_list

    else:
        logging.debug("Not type I polyketide synthase ('t1pks'), nonribosomal synthetase ('nrps') or their hybird")

if __name__ == '__main__':
    main()
