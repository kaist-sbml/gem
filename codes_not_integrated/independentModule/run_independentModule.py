'''
Created by Hyun Uk Kim, Nov 2014
2015 Hyun Uk Kim

This file executes functions handling:
1) parsing of the template model when the template model is first loaded,
2) creating various dictionary files derived from MNXref

Current template model: "TemplateModel_MNXref.xml"
'''
from cobra.io.sbml import create_cobra_model_from_sbml_file
from independentModule import *
import pickle
import subprocess

#root, tempGenome, tempModel = get_tempInfo(orgName='sco')
#print root
#print tempGenome
#print tempModel

################################################################################
#Loading and parsing of the template model
#cobra_model, cobra_reaction_dic = readCobraModel(root, tempModel)
#cobra_model = pickle.load(open("Pickle_template_model_original.p","rb"))
#cobra_reaction_dic = pickle.load(open('./forChecking/cobra_reaction_dic.p','rb'))
#cobra_metabolite_dic = pickle.load(open("Pickle_cobra_metabolite_dic.p","rb"))
    
#Reading a genbank file of the genome for the template model
#tempGenome_locusTag_aaSeq_dict = readSeq_templateGenome(tempGenome , 'genbank')

#Parsing GPR of the template model
#tempModel_biggRxnid_locusTag_dict, tempModel_locusTag_aaSeq_dict, tempModel_biggRxnid_wo_gene_list, tempModel_biggRxnidwithGene_woSeq_list = parse_templateModel_gpr('%s/tempModel_biggRxnid_locusTag_dict.txt' %(root), root, cobra_reaction_dic, tempGenome_locusTag_aaSeq_dict)

#makeFasta_locusTag_aaSeq(root, tempModel_locusTag_aaSeq_dict)

#Generating a DB for the genes from the template model
#make_blastDB(root, '%s/tempModel_locusTag_aaSeq.fa' %(root))

#tempModel_locusTag_aaSeq_dict = pickle.load(open("./forChecking/tempModel_locusTag_aaSeq_dict.p","rb"))
################################################################################

################################################################################
#pickle_input_mnxr_rxnid('./forChecking/allDB_mnxr_dict.txt', './forChecking/kegg_mnxr_dict.txt', './forChecking/mnxr_kegg_dict.txt', './forChecking/bigg_mnxr_dict.txt')
#allDB_mnxr_dict = pickle.load(open('./input2/allDB_mnxr_dict.p','rb'))
#kegg_mnxr_dict = pickle.load(open("Pickle_kegg_mnxr_dict.p","rb"))
#mnxr_kegg_dict = pickle.load(open("Pickle_mnxr_kegg_dict.p","rb"))

#pickle_input_mnxr_rxn_info()
#mnxr_rxn_dict = pickle.load(open("Pickle_mnxr_rxn_dict.p","rb"))

#pickle_input_bigg_kegg_mnx_compoundID()
#mnxm_bigg_compound_dict = pickle.load(open("Pickle_mnxm_bigg_compound_dict.p","rb"))
#kegg_mnxm_compound_dict = pickle.load(open("Pickle_kegg_mnxm_compound_dict.p","rb"))
#mnxm_kegg_compound_dict = pickle.load( open("Pickle_mnxm_kegg_compound_dict.p","rb"))

#pickling_Input_MNXM_compoundInfo()
#mnxm_compoundInfo_dict = pickle.load(open("Pickle_mnxm_compoundInfo_dict.p","rb"))

#cobra_model = pickle.load(open('./input2/template_model_pruned_GPRCorrected.p', 'rb'))
#templateModel_bigg_mnxr_dict = pickle_templateModel_bigg_mnxr(cobra_model, allDB_mnxr_dict)

pickle_universal_model_for_gapfill()
################################################################################

