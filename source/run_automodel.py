'''
2014-2015
Hyun Uk Kim, Tilmann Weber, Kyu-Sang Hwang and Jae Yong Ryu
'''

#Wildcard imports should never be used in production code.
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
    check_existing_rxns,
    get_mnxr_using_kegg,
    extract_rxn_mnxm_coeff,
    add_nonBBH_rxn
)
from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
import os
import pickle
import sys
import time

start = time.time()

dirname = sys.argv[1]

#Create output folders
folders = ['1_blastp_results', '2_primary_metabolic_model', '3_temp_models', '4_complete_model']

for folder in folders:
   if not os.path.isdir('./%s/'%dirname+folder):
       os.makedirs('./%s/'%dirname+folder)


#List of input (static) files as pickles
###################################################################
#For model pruning phase
#Choose "eco" or "sco"
orgName = 'sco'
root, temp_fasta = get_temp_fasta(orgName)
print root
print temp_fasta

model = pickle.load(open('%s/model.p' %(root),'rb'))
tempModel_biggRxnid_locusTag_dict = pickle.load(open('%s/tempModel_biggRxnid_locusTag_dict.p' %(root),'rb'))

#For model augmentation  phase
print "loading pickle files of the parsed template model and its relevant genbank data.."
kegg_mnxr_dict = pickle.load(open('./input2/kegg_mnxr_dict.p','rb'))
mnxr_kegg_dict = pickle.load(open('./input2/mnxr_kegg_dict.p','rb'))

mnxr_rxn_dict = pickle.load(open('./input2/mnxr_rxn_dict.p','rb'))

bigg_mnxm_compound_dict = pickle.load(open('./input2/bigg_mnxm_compound_dict.p','rb'))
mnxm_bigg_compound_dict = pickle.load(open('./input2/mnxm_bigg_compound_dict.p','rb'))
kegg_mnxm_compound_dict = pickle.load(open('./input2/kegg_mnxm_compound_dict.p','rb'))
mnxm_kegg_compound_dict = pickle.load( open('./input2/mnxm_kegg_compound_dict.p','rb'))

mnxm_compoundInfo_dict = pickle.load(open('./input2/mnxm_compoundInfo_dict.p','rb'))

templateModel_bigg_mnxr_dict = pickle.load(open('./input2/templateModel_bigg_mnxr_dict.p','rb'))
###################################################################
                           
print "\n", "pruning phase starting..", "\n"
###################################################################
print "looking for a gbk file of a target genome.."
target_gbk = get_target_gbk(dirname)

print "reading genbank file of the target genome.."    
targetGenome_locusTag_ec_dict, targetGenome_locusTag_prod_dict = get_targetGenomeInfo(dirname, target_gbk, "genbank")

print "\n", "looking for a fasta file of a target genome..", "\n"
target_fasta = get_target_fasta(dirname)

print "generating a DB for the genes from the target genome.."
make_blastDB(dirname, query_fasta=target_fasta)

print "\n", "running BLASTP #1: genes in the target genome against genes in the template model.."
run_blastp(target_fasta='./%s/1_blastp_results/targetGenome_locusTag_aaSeq.fa' %dirname, blastp_result='./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %dirname, db_dir = '%s/tempBlastDB' %(root), evalue=1e-30)

print "\n", "running BLASTP #2: genes in the template model against genes in the target genome.."
run_blastp(target_fasta='%s/tempModel_locusTag_aaSeq.fa' %(root), blastp_result='./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %dirname, db_dir = './%s/1_blastp_results/targetBlastDB' %dirname, evalue=1e-30)

print "parsing the results of BLASTP #1.."
blastpResults_dict1 = parseBlaspResults('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %dirname, './%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %dirname)
 
print "parsing the results of BLASTP #2.."
blastpResults_dict2 = parseBlaspResults('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %dirname, './%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %dirname)

print "selecting the best hits for BLASTP #1.."
bestHits_dict1 = makeBestHits_dict('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %dirname)

print "selecting the best hits for BLASTP #2.."
bestHits_dict2 = makeBestHits_dict('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %dirname)

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
###################################################################

###################################################################
print "adding the nonBBH gene-associated reactions..."
rxnid_to_add_list = check_existing_rxns(kegg_mnxr_dict, templateModel_bigg_mnxr_dict, rxnid_info_dict)

mnxr_to_add_list = get_mnxr_using_kegg(rxnid_to_add_list, kegg_mnxr_dict)

rxnid_mnxm_coeff_dict = extract_rxn_mnxm_coeff(mnxr_to_add_list, mnxr_rxn_dict, mnxm_bigg_compound_dict, mnxm_kegg_compound_dict, mnxr_kegg_dict)

target_model = add_nonBBH_rxn(modelPrunedGPR, rxnid_info_dict, rxnid_mnxm_coeff_dict, rxnid_locusTag_dict, bigg_mnxm_compound_dict, kegg_mnxm_compound_dict, mnxm_compoundInfo_dict, targetGenome_locusTag_prod_dict)
###################################################################

#Output files
#Model reloading and overwrtting are necessary for model consistency:
#e.g., metabolite IDs with correct compartment suffices & accurate model stats
#This can also mask the effects of model error (e.g., undeclared metabolite ID)
#Cobrapy IO module seems to have an error for adding new reactions
write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(dirname, dirname, orgName))
target_model = create_cobra_model_from_sbml_file('./%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(dirname, dirname, orgName))
write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/%s_target_model_%s.xml' %(dirname, dirname, orgName))

#Output on screen
model = pickle.load(open('%s/model.p' %(root),'rb'))
print "Number of genes:", len(model.genes), "/", len(modelPruned.genes), "/", len(target_model.genes)
print "Number of reactions:", len(model.reactions), "/", len(modelPruned.reactions), "/", len(target_model.reactions)
print "Number of metabolites:",  len(model.metabolites), "/", len(modelPruned.metabolites), "/", len(target_model.metabolites)

fp1 = open('./%s/2_primary_metabolic_model/%s_target_model_reactions.txt' %(dirname, dirname), "w")
fp2 = open('./%s/2_primary_metabolic_model/%s_target_model_metabolites.txt' %(dirname, dirname), "w")
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


###################################################################
#Secondary metabolic modeling
import run_smr_generation

print "Elapsed time:", time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
