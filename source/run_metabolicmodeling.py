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

start = time.time()

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--model', dest='orgName', default='sco', choices=['eco','sco'], help="Specify a template model for the target modeling")
parser.add_argument('-s', '--smr', dest='smr_generation', default=False, choices=[True,False], help="Specify whether to run secondary metabolic modeling")
parser.add_argument('-o', '--output', dest='output', default='output', help="Specify output directory")
parser.add_argument('-e', '--ec', dest='eficaz', default=False, choices=['eficaz'], help="Run EC number prediction using EFICAz")
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
print options.outputfoldername

#List of input (static) files as pickles
###################################################################
#For model pruning phase
root, temp_fasta = get_temp_fasta(options.orgName)
print root
print temp_fasta

model = pickle.load(open('%s/model.p' %(root),'rb'))
tempModel_biggRxnid_locusTag_dict = pickle.load(open('%s/tempModel_biggRxnid_locusTag_dict.p' %(root),'rb'))

#For model augmentation  phase
print "loading pickle files of the parsed template model and its relevant genbank data.."
bigg_mnxr_dict = pickle.load(open('./input2/bigg_mnxr_dict.p','rb'))
kegg_mnxr_dict = pickle.load(open('./input2/kegg_mnxr_dict.p','rb'))
mnxr_kegg_dict = pickle.load(open('./input2/mnxr_kegg_dict.p','rb'))

mnxr_rxn_dict = pickle.load(open('./input2/mnxr_rxn_dict.p','rb'))

bigg_mnxm_compound_dict = pickle.load(open('./input2/bigg_mnxm_compound_dict.p','rb'))
mnxm_bigg_compound_dict = pickle.load(open('./input2/mnxm_bigg_compound_dict.p','rb'))
kegg_mnxm_compound_dict = pickle.load(open('./input2/kegg_mnxm_compound_dict.p','rb'))
mnxm_kegg_compound_dict = pickle.load( open('./input2/mnxm_kegg_compound_dict.p','rb'))

mnxm_compoundInfo_dict = pickle.load(open('./input2/mnxm_compoundInfo_dict.p','rb'))

template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(root),'rb'))
###################################################################
                           
print "\n", "pruning phase starting..", "\n"
###################################################################
print "looking for a gbk file of a target genome.."
target_gbk = get_target_gbk(options.output)

print "reading genbank file of the target genome.."    
targetGenome_locusTag_ec_dict, targetGenome_locusTag_prod_dict = get_targetGenomeInfo(options.output, target_gbk, "genbank", options)

print "\n", "looking for a fasta file of a target genome..", "\n"
target_fasta = get_target_fasta(options.output)

print "generating a DB for the genes from the target genome.."
make_blastDB(options.output, query_fasta=target_fasta)

print "\n", "running BLASTP #1: genes in the target genome against genes in the template model.."
run_blastp(target_fasta='./%s/1_blastp_results/targetGenome_locusTag_aaSeq.fa' %options.output, blastp_result='./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, db_dir = '%s/tempBlastDB' %(root), evalue=1e-30)

print "\n", "running BLASTP #2: genes in the template model against genes in the target genome.."
run_blastp(target_fasta='%s/tempModel_locusTag_aaSeq.fa' %(root), blastp_result='./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, db_dir = './%s/1_blastp_results/targetBlastDB' %options.output, evalue=1e-30)

print "parsing the results of BLASTP #1.."
blastpResults_dict1 = parseBlaspResults('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, './%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)
 
print "parsing the results of BLASTP #2.."
blastpResults_dict2 = parseBlaspResults('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, './%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

print "selecting the best hits for BLASTP #1.."
bestHits_dict1 = makeBestHits_dict('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)

print "selecting the best hits for BLASTP #2.."
bestHits_dict2 = makeBestHits_dict('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

print "selecting the bidirectional best hits.."
targetBBH_list, temp_target_BBH_dict = getBBH(bestHits_dict1, bestHits_dict2)


print "selecting genes that are not bidirectional best hits.."
nonBBH_list = get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list)
###################################################################

###################################################################
print "labeling reactions with nonhomologous genes to remove from the template model.."
rxnToRemove_dict = labelRxnToRemove(model, temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict)

print "removing reactions with nonhomologous genes from the template model.."
modelPruned, rxnToRemoveEssn_dict, rxnRemoved_dict, rxnRetained_dict = pruneModel(model, rxnToRemove_dict, 'gurobi')

print "\n", "correcting GPR associations in the template model.."
modelPrunedGPR = swap_locusTag_tempModel(modelPruned, temp_target_BBH_dict)
###################################################################


print "\n", "augmentation phase starting..", "\n"
###################################################################
print "creating various dictionary files for the nonBBH gene-associted reactions..."

#NOT USED AT THE MOMENT
#Four nested functions
#def get_species_locusTag(ncbi_geneid):
#def get_ECNumberList_from_locusTag(species_locusTag):
#def get_rxnid_from_ECNumber(enzymeEC):
#def get_rxnInfo_from_rxnid(rxnid):

targetGenome_locusTag_ec_nonBBH_dict = get_targetGenome_locusTag_ec_nonBBH_dict(targetGenome_locusTag_ec_dict, nonBBH_list)

#Two nested functions
#def get_rxnid_from_ECNumber(enzymeEC):
#def get_rxnInfo_from_rxnid(rxnid):
rxnid_info_dict, rxnid_locusTag_dict = make_all_rxnInfo_fromRefSeq(targetGenome_locusTag_ec_nonBBH_dict)

modelPrunedGPR_mnxr_list = get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, bigg_mnxr_dict) 
###################################################################

###################################################################
print "adding the nonBBH gene-associated reactions..."
rxnid_to_add_list = check_existing_rxns(kegg_mnxr_dict, modelPrunedGPR_mnxr_list, rxnid_info_dict)

mnxr_to_add_list = get_mnxr_using_kegg(rxnid_to_add_list, kegg_mnxr_dict)

rxnid_mnxm_coeff_dict = extract_rxn_mnxm_coeff(mnxr_to_add_list, mnxr_rxn_dict, mnxm_bigg_compound_dict, mnxm_kegg_compound_dict, mnxr_kegg_dict)

target_model = add_nonBBH_rxn(modelPrunedGPR, rxnid_info_dict, rxnid_mnxm_coeff_dict, rxnid_locusTag_dict, bigg_mnxm_compound_dict, kegg_mnxm_compound_dict, mnxm_compoundInfo_dict, targetGenome_locusTag_prod_dict, template_exrxnid_flux_dict, options.output)
###################################################################

#Output files
#Model reloading and overwrtting are necessary for model consistency:
#e.g., metabolite IDs with correct compartment suffices & accurate model stats
#This can also mask the effects of model error (e.g., undeclared metabolite ID)
#Cobrapy IO module seems to have an error for adding new reactions
write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(options.output, options.output, options.orgName))
target_model = create_cobra_model_from_sbml_file('./%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(options.output, options.output, options.orgName))
write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(options.output, options.output, options.orgName))

#Output on screen
model = pickle.load(open('%s/model.p' %(root),'rb'))
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


###################################################################
#Secondary metabolic modeling
if options.smr_generation:
    import run_smr_generation

print "Elapsed time:", time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
