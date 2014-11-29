'''
2013, 2014
by Hyun Uk Kim, Jae Yong Ryu and Kyu-Sang Hwang
'''

from Bio import SeqIO
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis import single_deletion
from cobra import Reaction
from cobra.manipulation.delete import prune_unused_metabolites
import copy
import os
import pickle
import subprocess
import sys

#Looks for .xml and .gb(k) files in the pre-defined folder
def get_tempInfo(orgName):
    for root, dirs, files in os.walk('./input1/%s/' %(orgName)):
        for file in files:
            if file.endswith('.fa'):
                tempGenome = os.path.join(root, file)
            	return root, tempGenome
        else:
            sys.exit(1)

def readSeq(gbkFile, FileType):
    fp = open('./temp1/targetGenome_locusTag_aaSeq.fa','w')
    fp2 = open('./temp1/targetGenome_locusTag_ec_dict.txt','w')
    targetGenome_locusTag_aaSeq_dict = {}
    targetGenome_locusTag_ec_dict = {}

#Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)
    for feature in record.features:
        if feature.type == 'CDS':

#Retrieving "locus_tag (i.e., ORF name)" for each CDS
            locusTag = feature.qualifiers['locus_tag']

#Note that the numbers of CDS and "translation" do not match.
#There are occasions that CDS does not have "translation".
            if 'translation' in feature.qualifiers:

#Retrieving "translation (i.e., amino acid sequences)" for each CDS
                translation = feature.qualifiers.get('translation')
                targetGenome_locusTag_aaSeq_dict[locusTag[0]] = translation[0]
                print >>fp, '>%s\n%s' % (str(locusTag[0]), str(translation[0]))

#Some locus_tag's have multiple EC_number's
	    for ec in feature.qualifiers:
		if ec == 'EC_number':
	            ecnum = feature.qualifiers.get('EC_number')
		    targetGenome_locusTag_ec_dict[locusTag[0]] = ecnum
		    print locusTag[0], targetGenome_locusTag_ec_dict[locusTag[0]]

    for key in targetGenome_locusTag_ec_dict.keys():
	print >>fp2, '%s\t%s' %(key, targetGenome_locusTag_ec_dict[key])
    fp.close()
    fp2.close()
    
    return targetGenome_locusTag_ec_dict


#Looks for .fa and .gbk  files in the pre-defined folder
def get_targetGenome():
    for root, dirs, files in os.walk('./temp1/'):
        for file in files:
            if file.endswith('.fa'):
                fastaFile = os.path.join(root, file)
	    if file.endswith('.gbk') or file.endswith('.gb'):
                gbkFile = os.path.join(root, file)
	if fastaFile and gbkFile:
	    return fastaFile, gbkFile
	else:
	    sys.exit(1)

#making database files using fasta files
def make_blastDB(query_fasta):  
    db_dir = './temp1/targetBlastDB'
    DBprogramName = './blastpfiles/makeblastdb.exe'
    subprocess.call([DBprogramName,'-in',query_fasta,'-out',db_dir,'-dbtype','prot'])  

#Checks if DB is properly created; otherwise shutdown
    if os.path.isfile('./temp1/targetBlastDB.psq') == False:
	print "error in make_blastDB: blast DB not created"
	sys.exit(1)


#Output: b0002,ASPK|b0002,0.0,100.00,820
#"1e-30" is set as a threshold for bidirectional best hits
def run_blastp(target_fasta = '', blastp_result = '', db_dir = '', evalue = 1e-30):
    BLASTPprogramName = './blastpfiles/blastp.exe'  
    subprocess.call([BLASTPprogramName,'-query',target_fasta,'-out',blastp_result,'-db',db_dir,'-evalue', str(evalue),'-outfmt',"10 qseqid sseqid evalue score length pident"])


#Input: Results file from "run_blastp"
#Output: '\t' inserted between each element of the input
def parseBlaspResults(inputFile, outputFile):
    blastpResults_dict = {}
    fp = open(inputFile,'r')
    fp2 = open(outputFile,'w')
    itemnum=0
    for line in fp:                                               
        key = itemnum
        sptList = line.strip().split(',')                        
        qseqid = sptList[0].strip()
        sseqid = sptList[1].strip()
        evalue = sptList[2].strip()
        score = float(sptList[3].strip())
        length = int(sptList[4].strip())
        pident = float(sptList[5].strip())

        blastpResults_dict[key] = {"query_locusTag": qseqid, "db_locusTag": sseqid, "evalue": evalue, "score": score, "length": length, "identity": pident}
        print >>fp2, '%s\t%s\t%s\t%f\t%d\t%f' % (qseqid, sseqid, evalue, score, length, pident)
        itemnum += 1
    fp.close()
    fp2.close()
    return blastpResults_dict


#searching the best hit of a particular gene to a target genome
#Input: Results file from "parseBlaspResults"
#Output: query_locusTag '\t' db_locusTag
def makeBestHits_dict(inputFile):
    bestHits_dict = {}
    fp1 = open(inputFile,'r')
    
    for line in fp1:
        sptList = line.strip().split('\t')
        query_locusTag = sptList[0].strip()
        db_locusTag = sptList[1].strip()
        
        if query_locusTag not in bestHits_dict.keys():
#Value is in List to enable "append" below
            bestHits_dict[query_locusTag] = ([db_locusTag])

#This additional condition is necessary because blastp strangely produces a redundant set of gene pairs: e.g., SCO5892	['SAV_7184', 'SAV_419', 'SAV_419', 'SAV_419', 'SAV_419', ...]
        elif '%s' %(db_locusTag) not in bestHits_dict[query_locusTag]:
            bestHits_dict[query_locusTag].append((db_locusTag))
        
    fp1.close()
    return bestHits_dict


#Finding bidirectional best hits
#Input: two dict data from "selectBestHits" (e.g.,bestHits_dict) 
def getBBH(dic1,dic2, outputFile1, outputFile2):
    targetBBH_list = []
    temp_target_BBH_dict = {}
    fp1 = open(outputFile1,'w')
    fp2 = open(outputFile2,'w')

    for target_locusTag in dic1.keys():
        for temp_locusTag in dic1[target_locusTag]:
            if temp_locusTag in dic2.keys():
                for target_locusTag2 in dic2[temp_locusTag]: 
#The BBH case
                    if target_locusTag == target_locusTag2:
			if target_locusTag not in targetBBH_list:
			    targetBBH_list.append(target_locusTag)

#Some genes in template model have more than one homologous gene in a target genome
		        if temp_locusTag not in temp_target_BBH_dict.keys():
			    temp_target_BBH_dict[temp_locusTag] = ([target_locusTag])
			else:
			    temp_target_BBH_dict[temp_locusTag].append((target_locusTag))
#Output as .txt 
    for eachBBH in targetBBH_list:
	print >>fp1, '%s' %(eachBBH)

    for each_BBH in temp_target_BBH_dict.keys():
	print >>fp2, '%s\t%s' %(each_BBH, temp_target_BBH_dict[each_BBH])
           
    fp1.close()
    fp2.close()
    return targetBBH_list, temp_target_BBH_dict


#A set of locusTag not included in BBH_list were considered nonBBH_list.
#Their respective reactions, if available, are added to the model in augPhase.
def get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list, outputFile1):
    nonBBH_list = []
    fp1 = open(outputFile1,'w')

    for locusTag in targetGenome_locusTag_ec_dict.keys():
	if locusTag not in targetBBH_list:
	    nonBBH_list.append(locusTag)

    nonBBH_list = sorted(set(nonBBH_list))

    for each_nonBBH in nonBBH_list:
	print >>fp1, '%s' %(each_nonBBH)
    fp1.close()
    return nonBBH_list


def calcBoolean(booleanList):
    booleanList2 = copy.deepcopy(booleanList)
    finalList = []
    if len(booleanList) == 0:
        return False
    
    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):
                booleanList2[i][j] = 1
            value=1        
            for j in range(len(booleanList[i])):
                value = value * booleanList2[i][j]
            finalList.append(value)               
        else:            
            booleanList2[i] = 1
            finalList.append(booleanList2[i])   
    
    value=0        
    for i in range(len(finalList)):
        value = value + finalList[i]
    
    if value == 0:
        return False
    else:
        return True

                         
def makeBooleanFormat(temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict): 
    booleanList = tempModel_biggRxnid_locusTag_dict
    valueList = copy.deepcopy(booleanList)     
    booleanFormatList = copy.deepcopy(booleanList)
 
    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):
                geneid = booleanList[i][j]
                if geneid in temp_target_BBH_dict.keys():
		    value = 0
                else:
                    value = 'na'
                    booleanList[i][j]='na'  
                valueList[i][j] = value
        else:
            geneid = booleanList[i]
            if geneid in temp_target_BBH_dict.keys():
		value = 0
            else:
                value = 'na'
                booleanList[i]='na'
            valueList[i] = value 
    
    for i in range(len(booleanList)):
        if type(booleanList[i]) == list:
            while 'na' in booleanList[i]:
                booleanList[i].pop( booleanList[i].index('na') )
            if len(booleanList[i]) == 1:
                booleanList[i]=booleanList[i][0] 
    
    while 'na' in booleanList or [] in booleanList:
        if 'na' in booleanList:
            booleanList.pop( booleanList.index('na') )
        if [] in booleanList:
            booleanList.pop( booleanList.index([]) )
            
    newbooleanList = []
    for i in range(len(booleanList)):
        if type(booleanList[i]) != list:
            newbooleanList.append(booleanList[i])
    newbooleanList = list(set(newbooleanList))

    for i in range(len(booleanList)):
        if type(booleanList[i]) == list:
            tmpList = booleanList[i]
            tmpList = list(set(tmpList))
            newbooleanList.append(tmpList)
    
    booleanList = newbooleanList
    valueList2 = copy.deepcopy(booleanList)
    booleanFormatList = copy.deepcopy(booleanList)
    
    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):                
                geneid = booleanList[i][j]
		value = 0
                valueList2[i][j] = value
        else:            
            geneid = booleanList[i]
	    value = 0
            valueList2[i] = value     
    return valueList2


def labelRxnToRemove(model, temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict, outputFile1):
    fp1 = open(outputFile1, 'w')
    rxnToRemove_dict = {}

    for biggRxnid in tempModel_biggRxnid_locusTag_dict.keys():
	rxn = model.reactions.get_by_id(biggRxnid)
#Prevent removal of transport reactions from the template model
	if 'Transport' not in rxn.name and 'transport' not in rxn.name and 'Exchange' not in rxn.name and 'exchange' not in rxn.name:
            booleanList = makeBooleanFormat(temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict[biggRxnid])
            rxnToRemove_dict[biggRxnid] = calcBoolean(booleanList)

    for rxn in rxnToRemove_dict.keys():
	print >>fp1, '%s\t%s' %(rxn, rxnToRemove_dict[rxn])
    fp1.close()
    return rxnToRemove_dict


def pruneModel(model, rxnToRemove_dict, solver_org, outputFile1, outputFile2, outputFile3):
    fp3 = open(outputFile1, 'w')
    fp4 = open(outputFile2, 'w')
    fp5 = open(outputFile3, 'w')

    rxnToRemoveEssn_dict = {}
    rxnRemoved_dict = {}
    rxnRetained_dict = {}
    
    for rxnid in rxnToRemove_dict.keys():

#Single reaction deletion is performed only for reactions labelled as "False"
        if rxnToRemove_dict[rxnid] == False: 
            growth_rate_dict, solution_status_dict, problem_dict = single_deletion(model, list([rxnid]), element_type='reaction', solver=solver_arg) 

#Checks optimality first.
            if str(solution_status_dict.values()[0]) == 'optimal':

#Full list of reactions and predicted growth rates upon their deletions
                rxnToRemoveEssn_dict[rxnid] = float(growth_rate_dict.values()[0])
                print >>fp3, '%s\t%s' %(rxnid, rxnToRemoveEssn_dict[rxnid])
                
#Checks growth rate upon reaction deletion
                if float(growth_rate_dict.values()[0]) >= 0.01:
                    model.remove_reactions(rxnid)
#List of reactions removed from the template model
                    rxnRemoved_dict[rxnid] = float(growth_rate_dict.values()[0])
                    print "Removed reaction:", rxnid, growth_rate_dict.values()[0], len(model.reactions), len(model.metabolites)
		    print >>fp4, '%s\t%s\t%s' %(rxnid, rxnRemoved_dict[rxnid], len(model.reactions))
                else:
#List of reactions retained in the template model
                    rxnRetained_dict[rxnid] = float(growth_rate_dict.values()[0])
                    print "Retained reaction:", rxnid, growth_rate_dict.values()[0], len(model.reactions), len(model.metabolites)
		    print >>fp5, '%s\t%s\t%s' %(rxnid, rxnRetained_dict[rxnid], len(model.reactions))

#Removing metabolites that are not used in the reduced model
    prune_unused_metabolites(model)               
    modelPruned = copy.deepcopy(model)    

    fp3.close()
    fp4.close()
    fp5.close()
    return modelPruned, rxnToRemoveEssn_dict, rxnRemoved_dict, rxnRetained_dict


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

def swap_locusTag_tempModel(modelPruned, temp_target_BBH_dict, outputFile1):
    fp1 = open(outputFile1, 'w')

#Retrieves reactions associated with each homologous gene in template model
    for BBHrxn in modelPruned.reactions:
	booleanList = []
#Retrieves all the genes associated with a reaction having the homologous gene and transforms String to List
	booleanList = get_gpr_fromString_toList(BBHrxn.gene_reaction_rule)

        modified_booleanList = []
	for tempLocusTag in booleanList:

#Checks if the element itself is List.
#If the element is not List, then gene in template model is directly replaced with genes in target genome
	    if type(tempLocusTag) != list:
		if temp_target_BBH_dict.has_key(tempLocusTag) == True:
		    booleanList.pop(booleanList.index(tempLocusTag))
		    for targetLocusTag in temp_target_BBH_dict[tempLocusTag]:
#			booleanList.append(targetLocusTag)
                        modified_booleanList.append(targetLocusTag)
                else:
                    modified_booleanList.append( tempLocusTag )

#This is the case the element is List
	    else:
                temp_gpr_list = []
		for eachLocusTag in tempLocusTag:
		    if temp_target_BBH_dict.has_key(eachLocusTag) == True:
			for targetLocusTag in temp_target_BBH_dict[eachLocusTag]:
			    temp_gpr_list.append(targetLocusTag)
                    else:
                        temp_gpr_list.append(eachLocusTag)
#This case was not generated, but just in case
                if len(temp_gpr_list)==1:
		    print temp_gpr_list
                    modified_booleanList.append(temp_gpr_list[0])
                elif len(temp_gpr_list) > 1:
                    modified_booleanList.append( temp_gpr_list )


#Converts GPR in List to String:
        booleanList = copy.deepcopy( modified_booleanList )
	stringlist = []
	for i in range(len(booleanList)):
	    booleanstring = ''
	    if type(booleanList[i])==list:
		booleanstring=' and '.join(booleanList[i])
		booleanstring='('+booleanstring+')'
		stringlist.append(booleanstring)
	    else:
		booleanstring=booleanList[i]
		stringlist.append(booleanstring)
		
	booleanstring2 = ''
	if len(stringlist) > 0:
	    booleanstring2 = '('+' or '.join(stringlist)+')'
	BBHrxn.add_gene_reaction_rule(booleanstring2)
	print >>fp1, '%s\t%s' %(BBHrxn, BBHrxn.gene_reaction_rule)

    fp1.close()       
    modelPrunedGPR = copy.deepcopy(modelPruned)            
    return modelPrunedGPR

