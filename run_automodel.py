'''
2014
Hyun Uk Kim
'''

from prunPhase import *
from augPhase import *
from cobra.io.sbml import write_cobra_model_to_sbml_file
import copy
import pickle
import urllib2
import time

start = time.time()

#Names of output txt files and their pickesl from functions are given herein.
###################################################################
#from "independentModule_"

#Choose "eco" or "sco"
orgName = 'sco'
root, tempGenome = get_tempInfo(orgName)
print root
print tempGenome

#print "loading pickle files of the parsed template model and its relevant genbank data.."
model = pickle.load(open('%s/model.p' %(root),'rb'))
tempModel_biggRxnid_locusTag_dict = pickle.load(open('%s/tempModel_biggRxnid_locusTag_dict.p' %(root),'rb'))
###################################################################
                           
###################################################################
#print "looking for fasta and genbank files for a target genome.."
fastaFile, gbkFile = get_targetGenome()
#print "complete"+'\n'

#print "reading genbank file of the target genome.."    
targetGenome_locusTag_ec_dict = readSeq(gbkFile, "genbank")
#print "complete"+'\n'

print "\n"
print "pruning phase starting.."
print "\n"

#print "generating a DB for the genes from the target genome.."
make_blastDB(query_fasta=fastaFile)
#print "complete"+'\n'

#print "running BLASTP #1: genes in the target genome against genes in the template model.."
run_blastp(target_fasta='./temp1/targetGenome_locusTag_aaSeq.fa', blastp_result='./temp1/blastp_targetGenome_against_tempGenome.txt', db_dir = '%s/tempBlastDB' %(root), evalue=1e-30)
#print "complete"+'\n'    

#print "running BLASTP #2: genes in the template model against genes in the target genome.."
run_blastp(target_fasta='%s/tempModel_locusTag_aaSeq.fa' %(root), blastp_result='./temp1/blastp_tempGenome_against_targetGenome.txt', db_dir = './temp1/targetBlastDB', evalue=1e-30)
#print "complete"+'\n'

#print "parsing the results of BLASTP #1.."
blastpResults_dict1 = parseBlaspResults('./temp1/blastp_targetGenome_against_tempGenome.txt', './temp1/blastp_targetGenome_against_tempGenome_parsed.txt')
#print "complete"+'\n'
 
#print "parsing the results of BLASTP #2.."
blastpResults_dict2 = parseBlaspResults('./temp1/blastp_tempGenome_against_targetGenome.txt', './temp1/blastp_tempGenome_against_targetGenome_parsed.txt')
#print "complete"+'\n'

#print "selecting the best hits for BLASTP #1.."
bestHits_dict1 = makeBestHits_dict(inputFile='./temp1/blastp_targetGenome_against_tempGenome_parsed.txt')
#print "complete"+'\n'

#print "selecting the best hits for BLASTP #2.."
bestHits_dict2 = makeBestHits_dict(inputFile='./temp1/blastp_tempGenome_against_targetGenome_parsed.txt')
#print "complete"+'\n'

#print "selecting the bidirectional best hits.."
targetBBH_list, temp_target_BBH_dict = getBBH(bestHits_dict1, bestHits_dict2, './temp1/blastp_BBH_targetLocusTag_list.txt', './temp1/blastp_BBH_temp_targetLocusTag_dict.txt')
#print "complete"+'\n'

#print "selecting genes that are not bidirectional best hits.."
nonBBH_list = get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list, './temp1/blastp_nonBBH.txt')
#print "complete"+'\n'

###################################################################

###################################################################
#print "labeling reactions with nonhomologous genes to remove from the template model.."
rxnToRemove_dict = labelRxnToRemove(model, temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict, './temp1/rxnToRemove.txt')
#print "complete"+'\n'

#print "removing reactions with nonhomologous genes from the template model.."
modelPruned, rxnToRemoveEssn_dict, rxnRemoved_dict, rxnRetained_dict = pruneModel(model, rxnToRemove_dict, 'gurobi', './temp1/rxnToRemove_essen.txt', './temp1/rxnRemoved.txt', './temp1/rxnRetained.txt')
#print "complete"+'\n'

pickle.dump(modelPruned, open('./temp1/modelPruned.p', 'wb'))
modelPruned = pickle.load(open('./temp1/modelPruned.p', 'rb'))

#print "correcting GPR associations in the template model.."
modelPrunedGPR = swap_locusTag_tempModel(modelPruned, temp_target_BBH_dict, './temp1/rxnGPRcorrected.txt')
pickle.dump(modelPrunedGPR, open('./temp1/modelPrunedGPR.p', 'wb'))
#print "complete"+'\n'
###################################################################


print "augmentation phase starting.."
###################################################################
#from "independentModule_"
#print "loading pickle files of the parsed template model and its relevant genbank data.."
modelPrunedGPR = pickle.load(open('./temp1/modelPrunedGPR.p', 'rb'))

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

#print "creacting dictionary files for the noBBH genes..."
###################################################################
locusTag_geneID_dict, geneID_locusTag_dict = make_locusTag_geneID_nonBBH(gbkFile, "genbank","./temp1/blastp_nonBBH.txt", "./temp2/locusTag_geneID_dict.txt")

pickle.dump(locusTag_geneID_dict, open('./temp2/locusTag_geneID_dict.p', 'wb'))
pickle.dump(geneID_locusTag_dict, open('./temp2/geneID_locusTag_dict.p', 'wb'))
locusTag_geneID_dict = pickle.load(open('./temp2/locusTag_geneID_dict.p', 'rb'))
geneID_locusTag_dict = pickle.load(open('./temp2/geneID_locusTag_dict.p', 'rb'))
###################################################################

#print "creating various dictionary files for the nonBBH gene-associted reactions..."
###################################################################
#Four nested functions
#def get_species_locusTag(ncbi_geneid):
#def get_ECNumberList_from_locusTag(species_locusTag):
#def get_rxnid_from_ECNumber(enzymeEC):
#def get_rxnInfo_from_rxnid(rxnid):
rxnid_info_dict, rxnid_geneid_dict, rxnid_locusTag_dict = make_all_rxnInfo(locusTag_geneID_dict, "./temp2/locusTag_list_withoutReaction.txt", "./temp2/rxnid_info_dict.txt", "./temp2/rxnid_geneid_dict.txt", "./temp2/rxnid_locusTag_dict.txt")

pickle.dump(rxnid_info_dict, open('./temp2/rxnid_info_dict.p', 'wb'))
pickle.dump(rxnid_geneid_dict, open('./temp2/rxnid_geneid_dict.p', 'wb'))
pickle.dump(rxnid_locusTag_dict, open('./temp2/rxnid_locusTag_dict.p', 'wb'))
rxnid_info_dict = pickle.load(open('./temp2/rxnid_info_dict.p', 'rb'))
rxnid_geneid_dict = pickle.load(open('./temp2/rxnid_geneid_dict.p', 'rb'))
rxnid_locusTag_dict = pickle.load(open('./temp2/rxnid_locusTag_dict.p', 'rb'))
###################################################################

#print "adding the nonBBH gene-associated reactions..."
rxnid_to_add_list = check_existing_rxns(kegg_mnxr_dict, templateModel_bigg_mnxr_dict, rxnid_info_dict, "./temp2/rxnid_overlapping.txt", "./temp2/rxnid_to_add.txt")

mnxr_to_add_list = get_mnxr_using_kegg(rxnid_to_add_list, kegg_mnxr_dict, "./temp2/mnxr_to_add.txt")

rxnid_mnxm_coeff_dict = extract_rxn_mnxm_coeff(mnxr_to_add_list, mnxr_rxn_dict, mnxm_bigg_compound_dict, mnxm_kegg_compound_dict, mnxr_kegg_dict, './temp2/rxnid_mnxr_reaction.txt', './temp2/discarded_rxns.txt')


#One nested function
#get_compoundInfo(compoundID)
target_model = add_nonBBH_rxn(modelPrunedGPR, rxnid_info_dict, rxnid_mnxm_coeff_dict, rxnid_locusTag_dict, bigg_mnxm_compound_dict, kegg_mnxm_compound_dict, mnxm_compoundInfo_dict, "./temp2/nonBBH_final_reactions.txt")

pickle.dump(target_model, open('./temp2/target_model_%s.p' %(orgName), 'wb'))

model = pickle.load(open('%s/model.p' %(root),'rb'))
print "Stats of the two models: template model vs. pruned intermediate model vs . target_model"
print "Number of genes:", len(model.genes), "/", len(modelPruned.genes), "/", len(target_model.genes)
print "Number of reactions:", len(model.reactions), "/", len(modelPruned.reactions), "/", len(target_model.reactions)
print "Number of metabolites:",  len(model.metabolites), "/", len(modelPruned.metabolites), "/", len(target_model.metabolites)

write_cobra_model_to_sbml_file(target_model, './temp2/target_model_%s.xml' %(orgName))

print "elpased time:", time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
