'''
by Hyun Uk Kim, Jan - Nov 2014
2015 Hyun Uk Kim
'''

from Bio import SeqIO
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
import os
import pickle
import re
import subprocess
import sys
import urllib2


#Looks for .xml and .gb(k) files in the pre-defined folder
def get_tempInfo(orgName):
    for root, dirs, files in os.walk('./input1/%s/' %(orgName)):
        for file in files:
            if file.endswith('.gb') or file.endswith('.gbk'):
		tempGenome = os.path.join(root, file)
            if file.endswith('.xml'):
		tempModel = os.path.join(root, file)
	if tempGenome and tempModel:
	    return root, tempGenome, tempModel
	else:
	    sys.exit(1)

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
    pickle.dump(cobra_reaction_dic, open('./forChecking/cobra_reaction_dic.p','wb'))

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
    tempModel_biggRxnid_wo_gene_list = []
    tempModel_biggRxnidwithGene_woSeq_list = [] #For reactions with genes, but no seq info

    #Starts model parsing
    for Reaction_name in cobra_reaction_dic.keys():
        gpr = cobra_reaction_dic[Reaction_name]['gpr']

        #Parsing gene associations into a list
        Genes_list = get_gpr_fromString_toList(gpr)

        #Saves reaction names and their intact GPR associations
        #Input for "makeBooleanFormat" in "Reconstruction_GPRCalculation_Edited.py"
        tempModel_biggRxnid_locusTag_dict[Reaction_name] = Genes_list

        for Genes_each_list in Genes_list:
            #Stores reactions without genes
            if not Genes_each_list:
                tempModel_biggRxnid_wo_gene_list.append(Reaction_name)
 
            #Checks if the element itself is List.
            #'if type(Genes_each_list) == list' also works.
            #Genes connected with 'AND' come in the List type.
            elif isinstance(Genes_each_list, list):
                for Genes_each2_list in Genes_each_list:
                    if Genes_each2_list in tempGenome_locusTag_aaSeq_dict.keys():
                        #Dictionary - ORF:AA seq
                        #Stores only reactions with genes and their amino acid seq in Dictionary
                        if tempGenome_locusTag_aaSeq_dict[Genes_each2_list]:
                            tempModel_locusTag_aaSeq_dict[Genes_each2_list] = tempGenome_locusTag_aaSeq_dict[Genes_each2_list]

                    #Some reactions do not have seq info despite their presence of genes (e.g., b2092)
                    else:
                        tempModel_biggRxnidwithGene_woSeq_list.append(Reaction_name)

            #Single genes for a reaction, or genes connected with 'OR'
	    else:
                if Genes_each_list in tempGenome_locusTag_aaSeq_dict.keys():

                    #Dictionary - ORF:AA seq
                    #Stores only reactions with genes and their amino acid seq in Dictionary
                    if tempGenome_locusTag_aaSeq_dict[Genes_each_list]:
                        tempModel_locusTag_aaSeq_dict[Genes_each_list] = tempGenome_locusTag_aaSeq_dict[Genes_each_list]

                #Some reactions in iAF1260 do not have aa seq despite their presence of genes (e.g., b2092)
                else:
                    tempModel_biggRxnidwithGene_woSeq_list.append(Reaction_name)

    pickle.dump(tempModel_biggRxnid_locusTag_dict, open('%s/tempModel_biggRxnid_locusTag_dict.p' %(root),'wb'))
    pickle.dump(tempModel_locusTag_aaSeq_dict, open('./forChecking/tempModel_locusTag_aaSeq_dict.p','wb'))
    pickle.dump(tempModel_biggRxnid_wo_gene_list, open('./forChecking/tempModel_biggRxnid_wo_gene_list.p','wb'))
    pickle.dump(tempModel_biggRxnidwithGene_woSeq_list, open('./forChecking/tempModel_biggRxnidwithGene_woSeq_list.p','wb'))

    for key in tempModel_biggRxnid_locusTag_dict.keys():
	print >>fp1, '%s\t%s' %(key, tempModel_biggRxnid_locusTag_dict[key])
    fp1.close()
    return tempModel_biggRxnid_locusTag_dict, tempModel_locusTag_aaSeq_dict, tempModel_biggRxnid_wo_gene_list, tempModel_biggRxnidwithGene_woSeq_list


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


#Get reaction information using its ID from KEGG
#Input: KEGG rxnid in string form (e.g., R00362)
#Output: Reaction information for 'Name', 'Definition', and 'Equation' as dictionary form {'NAME': 'citrate oxaloacetate-lyase', 'DEFINITION': Citrate <=> Acetate + Oxaloacetate, 'EQUATION': C00158 <=> C00033 + C00036}
def get_rxnInfo_from_rxnid(rxnid):  
    url = "http://rest.kegg.jp/get/rn:%s"%(rxnid)
    reaction_info_text = urllib2.urlopen(url).read()
    split_text = reaction_info_text.strip().split('\n')
    NAME = ''
    DEFINITION = ''
    EQUATION = ''
    ENZYME = ''
    PATHWAY = ''
    for line in split_text:
        sptlist = line.split()
        if sptlist[0].strip() == 'NAME':
            NAME = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'DEFINITION':
            DEFINITION = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'EQUATION':
            EQUATION = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'ENZYME':
            ENZYME = ' '.join(sptlist[1:])

#Considers only reactions mapped in pathways
        if sptlist[0].strip() == 'PATHWAY':        
            PATHWAY = ' '.join(sptlist[1:])
	    return {'NAME':NAME, 'DEFINITION':DEFINITION, 'EQUATION':EQUATION, 'ENZYME':ENZYME, 'PATHWAY':PATHWAY}


def pickle_input_mnxr_rxnid(outputFile1, outputFile2, outputFile3, outputFile4):
    fp1 = open('Input_MNXR_reactionID_v1_1.tsv',"r")
#    fp2 = open(outputFile1, "w")
#    fp3 = open(outputFile2, "w")
#    fp4 = open(outputFile3, "w")
    fp5 = open(outputFile4, "w")

    allDB_mnxr_dict = {}
    kegg_mnxr_dict = {}
    mnxr_kegg_dict = {}
    bigg_mnxr_dict = {}

    rxnid = fp1.readline()
    while rxnid:
	rxnid = rxnid.split('\t')
	rxnid[0] = rxnid[0].strip()
	rxnid[1] = rxnid[1].strip()
	rxnid[2] = rxnid[2].strip()

        #Output = {reaction IDs in all the DB except for KEGG:MNXR}	
#	if rxnid[0] != 'kegg':
#	    allDB_mnxr_dict[rxnid[1]] = rxnid[2]

        #Output = {KEGG reaction ID:MNXR}
#	elif rxnid[0] == 'kegg':
	    #print rxnid[0], rxnid[1]
            #Some reaction ID do not exist at KEGG, causing HTTPError
#	    try:
                #Some different reaction IDs assigned to the same MNXR
                #Only reactions IDs with "PATHWAY" are considered by implementing this function
#	        if get_rxnInfo_from_rxnid(rxnid[1]):
#	            kegg_mnxr_dict[rxnid[1]] = rxnid[2]
#	            mnxr_kegg_dict[rxnid[2]] = rxnid[1]
#	    except:
#	        print "Not available @ KEGG:", rxnid[1]

        #Output = {bigg reactionID:MNXR}
        if rxnid[0] == 'bigg':
            bigg_mnxr_dict[rxnid[1]] = rxnid[2]

	print rxnid[1], rxnid[2]
	rxnid = fp1.readline()

#    pickle.dump(allDB_mnxr_dict, open('./forChecking/allDB_mnxr_dict.p','wb'))
#    pickle.dump(kegg_mnxr_dict, open('./input2/kegg_mnxr_dict.p','wb'))
#    pickle.dump(mnxr_kegg_dict, open('./input2/mnxr_kegg_dict.p','wb'))
    pickle.dump(bigg_mnxr_dict, open('./input2/bigg_mnxr_dict.p','wb'))
    
#    for key in allDB_mnxr_dict.keys():
#	print >>fp2, '%s\t%s' %(key, allDB_mnxr_dict[key]) 
#    for key in allDB_mnxr_dict.keys():
#	print >>fp3, '%s\t%s' %(key, kegg_mnxr_dict[key]) 
#    for key in allDB_mnxr_dict.keys():
#	print >>fp4, '%s\t%s' %(key, mnxr_kegg_dict[key]) 
    for key in bigg_mnxr_dict.keys():
	print >>fp5, '%s\t%s' %(key, bigg_mnxr_dict[key]) 

    fp1.close()
#    fp2.close()
#    fp3.close()
#    fp4.close()
    fp5.close()


def pickling_Input_MNXreaction():
    fp1 = open('Input_MNXreaction.tsv',"r")
    mnxr_rxn_dict = {}
    rxn = fp1.readline()

    while rxn:
	rxn = rxn.split('\t')
	rxn[0] = rxn[0].strip()
	rxn[1] = rxn[1].strip()

	print rxn[0], rxn[1]
	mnxr_rxn_dict[rxn[0]] = rxn[1]
	rxn = fp1.readline()
    pickle.dump(mnxr_rxn_dict, open('./input2/mnxr_rxn_dict.p','wb'))
    fp1.close()


#Based on fix_legacy_id(id, use_hyphens=False, fix_compartments=False) of cobrapy:
def replace_special_characters_compoundid(biggid):
    biggid = biggid.replace('-', '_DASH_')
    biggid = biggid.replace('/', '_FSLASH_')
    biggid = biggid.replace("\\",'_BSLASH_')
    biggid = biggid.replace('(', '_LPAREN_')
    biggid = biggid.replace('[', '_LSQBKT_')
    biggid = biggid.replace(']', '_RSQBKT_')
    biggid = biggid.replace(')', '_RPAREN_')
    biggid = biggid.replace(',', '_COMMA_')
    biggid = biggid.replace('.', '_PERIOD_')
    biggid = biggid.replace("'", '_APOS_')
    biggid = biggid.replace('&', '&amp;')
    biggid = biggid.replace('<', '&lt;')
    biggid = biggid.replace('>', '&gt;')
    biggid = biggid.replace('"', '&quot;')

    return biggid

def pickle_input_bigg_kegg_mnx_compoundID():
    fp1 = open('Input_BiGG_KEGG_MNX_compoundID_v1_3.tsv',"r")
    bigg_mnxm_compound_dict = {}
    mnxm_bigg_compound_dict = {}
    kegg_mnxm_compound_dict = {}
    mnxm_kegg_compound_dict = {}
    text = fp1.readline()

    while text:
  	text = text.split('\t')
	text[1] = text[1].strip()
	text[2] = text[2].strip()
	if text[0] == 'bigg':
            biggid = replace_special_characters_compoundid(text[1])
	    bigg_mnxm_compound_dict[biggid] = text[2]
	    mnxm_bigg_compound_dict[text[2]] = biggid
	elif text[0] == 'kegg':
	    kegg_mnxm_compound_dict[text[1]] = text[2]

            #Following conditions give priority to compoundID starting with 'C' than others
	    if 'D' not in text[1] and 'E' not in text[1] and 'G' not in text[1]:
	        mnxm_kegg_compound_dict[text[2]] = text[1]
	    else:
		if text[2] not in mnxm_kegg_compound_dict.keys():
	            mnxm_kegg_compound_dict[text[2]] = text[1]

	print text[0], text[1], text[2]
	text = fp1.readline()
    pickle.dump(bigg_mnxm_compound_dict, open('./input2/bigg_mnxm_compound_dict.p','wb'))
    pickle.dump(mnxm_bigg_compound_dict, open('./input2/mnxm_bigg_compound_dict.p','wb'))
    pickle.dump(kegg_mnxm_compound_dict, open('./input2/kegg_mnxm_compound_dict.p','wb'))
    pickle.dump(mnxm_kegg_compound_dict, open('./input2/mnxm_kegg_compound_dict.p','wb'))
    fp1.close()


#mnxm_compoundInfo_dict = {'MNXM128019': ['Methyl trans-p-methoxycinnamate', 'C11H12O3']}
def pickling_Input_MNXM_compoundInfo():
    fp1 = open('Input_MNXM_compoundInfo_v1_0.tsv',"r")
    mnxm_compoundInfo_dict = {}
    text = fp1.readline()

    while text:
	text = text.split('\t')
	text[0] = text[0].strip()
	text[1] = text[1].strip()
	text[2] = text[2].strip()
	mnxm_compoundInfo_dict[text[0]] = [(text[1])]
	mnxm_compoundInfo_dict[text[0]].append((text[2]))

	print text[0], text[1], text[2]
	text = fp1.readline()
    print mnxm_compoundInfo_dict
    pickle.dump(mnxm_compoundInfo_dict, open('./input2/mnxm_compoundInfo_dict.p','wb'))
    fp1.close()


#e.g., "EX_ca2_LPAREN_e_RPAREN_" => "EX_ca2(e)"
def fix_legacy_id(id, use_hyphens=False, fix_compartments=False):
    id = id.replace('_DASH_', '__')
    id = id.replace('_FSLASH_', '/')
    id = id.replace('_BSLASH_', "\\")
    id = id.replace('_LPAREN_', '(')
    id = id.replace('_LSQBKT_', '[')
    id = id.replace('_RSQBKT_', ']')
    id = id.replace('_RPAREN_', ')')
    id = id.replace('_COMMA_', ',')
    id = id.replace('_PERIOD_', '.')
    id = id.replace('_APOS_', "'")
    id = id.replace('&amp;', '&')
    id = id.replace('&lt;', '<')
    id = id.replace('&gt;', '>')
    id = id.replace('&quot;', '"')
    if use_hyphens:
        id = id.replace('__', '-')
    else:
        id = id.replace("-", "__")
    if fix_compartments:
        if len(id) > 2:
            if (id[-3] == "(" and id[-1] == ")") or \
               (id[-3] == "[" and id[-1] == "]"):
                id = id[:-3] + "_" + id[-2]
    return id


def get_templateModel_bigg_mnxr(cobra_model, allDB_mnxr_dict):
    fp1 = open('./input2/templateModel_bigg_mnxr_dict.txt','w')
    templateModel_bigg_mnxr_dict = {}

#Counts the number of metabolite IDs
    index_last = len(cobra_model.reactions)

#This is necessary becasue the index starts from "0", NOT "1"
    index_last = index_last - 1
    index = 0

    while index <= index_last:
        reaction = cobra_model.reactions[index].id
        biggid_edited = fix_legacy_id(reaction, use_hyphens=True, fix_compartments=False)

#Correcting prefixes: "D_" or "L_" => "D-" or "L-"
        if biggid_edited.find("D_") == 0 or biggid_edited.find("L_") == 0:
            sptlist = biggid_edited.split("_")
            biggid_edited = "-".join(sptlist)
        if biggid_edited in allDB_mnxr_dict.keys():

#Stores reaction ID with their MNX IDs having compartment tags in a dictionary form
            templateModel_bigg_mnxr_dict[biggid_edited] = allDB_mnxr_dict[biggid_edited]
	    fp1.write(str(biggid_edited)+'\t'+str(templateModel_bigg_mnxr_dict[biggid_edited])+'\n')
        index+=1

    pickle.dump(templateModel_bigg_mnxr_dict, open('./input2/templateModel_bigg_mnxr_dict.p',"wb"))
    fp1.close()
    return templateModel_bigg_mnxr_dic
