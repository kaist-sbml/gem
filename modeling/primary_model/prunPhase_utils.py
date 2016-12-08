'''
2014
Hyun Uk Kim, Tilmann Weber, Jae Yong Ryu and Kyu-Sang Hwang
'''

import copy
import logging
#import os
#import pickle
#import subprocess
#import sys
#from Bio import SeqIO
from cobra.flux_analysis import single_deletion
from cobra.manipulation.delete import prune_unused_metabolites
#from eficaz.__init__ import getECs


def calcBoolean(booleanList):
    booleanList2 = copy.deepcopy(booleanList)
    finalList = []
    threshold = 1

    if len(booleanList) == 0:
        return False

    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):

                #Threshold created to differentiate BBH and nonBBH genes.
                if float(booleanList[i][j]) >= threshold:
                    booleanList[i][j] = 0
                if float(booleanList[i][j]) <= threshold:
                    booleanList2[i][j] = 1
            value=1
            for j in range(len(booleanList[i])):
                value = value * booleanList2[i][j]
            finalList.append(value)
        else:
            #Threshold created to differentiate BBH and nonBBH genes.
            if float(booleanList[i]) >= threshold:
                booleanList[i] = 0
            if float(booleanList[i]) <= threshold:
                booleanList2[i] = 1
            finalList.append(booleanList2[i])

    value=0
    for i in range(len(finalList)):
        value = value + finalList[i]

    if value == 0:
        return False
    else:
        return True


#Output: e.g., [[2, '0'], ['0', 2]]
#Now considers nonBBH genes without removing them in the Boolean list.
#For ( A and B), if one of them is nonBBH, its rxn should be False (subject to removal).
def makeBooleanFormat(temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict):
    booleanList = tempModel_biggRxnid_locusTag_dict
    valueList = copy.deepcopy(booleanList)

    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):
                geneid = booleanList[i][j]
                if geneid in temp_target_BBH_dict.keys():
                    value = 0
                #Now considers nonBBH genes without removing them in the Boolean list
                else:
                    value = 2 #For nonBBH genes
                    booleanList[i][j]='nonBBH'
                valueList[i][j] = value
        else:
            geneid = booleanList[i]
            if geneid in temp_target_BBH_dict.keys():
		value = 0
            else:
                value = 'na'
                booleanList[i]='na'
            valueList[i] = value

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

    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):
                geneid = booleanList[i][j]
		if 'nonBBH' not in geneid:
                    value = 0
		else:
                    value = 2
                valueList2[i][j] = value
        else:
            geneid = booleanList[i]
            value = 0
            valueList2[i] = value
    return valueList2


def labelRxnToRemove(model, options):
    rxnToRemove_dict = {}

    for biggRxnid in options.tempModel_biggRxnid_locusTag_dict.keys():
	rxn = model.reactions.get_by_id(biggRxnid)
        #Prevent removal of transport reactions from the template model
	if 'Transport' not in rxn.name and 'transport' not in rxn.name and 'Exchange' not in rxn.name and 'exchange' not in rxn.name:
            booleanList = makeBooleanFormat(options.temp_target_BBH_dict, options.tempModel_biggRxnid_locusTag_dict[biggRxnid])
            rxnToRemove_dict[biggRxnid] = calcBoolean(booleanList)

    options.rxnToRemove_dict = rxnToRemove_dict


def pruneModel(model, options, solver_arg):
    rxnToRemoveEssn_dict = {}
    rxnRemoved_dict = {}
    rxnRetained_dict = {}

    for rxnid in options.rxnToRemove_dict.keys():

        #Single reaction deletion is performed only for reactions labelled as "False"
        if options.rxnToRemove_dict[rxnid] == False:
            growth_rate_dict, solution_status_dict, problem_dict = single_deletion(model, list([rxnid]), element_type='reaction', solver=solver_arg)

            #Checks optimality first.
            if str(solution_status_dict.values()[0]) == 'optimal':

                #Full list of reactions and predicted growth rates upon their deletions
                rxnToRemoveEssn_dict[rxnid] = float(growth_rate_dict.values()[0])

                #Checks growth rate upon reaction deletion
                if float(growth_rate_dict.values()[0]) >= 0.01:
                    model.remove_reactions(rxnid)
                    #List of reactions removed from the template model
                    rxnRemoved_dict[rxnid] = float(growth_rate_dict.values()[0])
                    logging.debug("Removed reaction: %s; %s; %s; %s" %(rxnid, growth_rate_dict.values()[0], len(model.reactions), len(model.metabolites)))
                else:
                    #List of reactions retained in the template model
                    rxnRetained_dict[rxnid] = float(growth_rate_dict.values()[0])
                    logging.debug("Retained reaction: %s; %s; %s; %s" %(rxnid, growth_rate_dict.values()[0], len(model.reactions), len(model.metabolites)))

    #Removing metabolites that are not used in the reduced model
    prune_unused_metabolites(model)
    modelPruned = copy.deepcopy(model)

    #rxnToRemoveEssn_dict, rxnRemoved_dict and rxnRetained_dict:
    #Not used in the downstream of this pipeline
    return modelPruned


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


def swap_locusTag_tempModel(modelPruned, options):

    #Retrieves reactions associated with each homologous gene in template model
    for BBHrxn in modelPruned.reactions:
	booleanList = []
        #Retrieves all the genes associated with a reaction having the homologous gene
        #and transforms String to List
	booleanList = get_gpr_fromString_toList(BBHrxn.gene_reaction_rule)

        modified_booleanList = []
	for tempLocusTag in booleanList:

            #Checks if the element itself is List.
            #If the element is not List, then gene in template model is
            #directly replaced with genes in target genome
            if type(tempLocusTag) != list:
		if tempLocusTag in options.temp_target_BBH_dict:
                    booleanList.pop(booleanList.index(tempLocusTag))
                    for targetLocusTag in options.temp_target_BBH_dict[tempLocusTag]:
                        modified_booleanList.append(targetLocusTag)
                else:
                    modified_booleanList.append( tempLocusTag )

            #This is the case the element is List
            else:
                temp_gpr_list = []
		for eachLocusTag in tempLocusTag:
                    if eachLocusTag in options.temp_target_BBH_dict:
			for targetLocusTag in options.temp_target_BBH_dict[eachLocusTag]:
                            temp_gpr_list.append(targetLocusTag)
                    else:
                        temp_gpr_list.append(eachLocusTag)
                #This case was not generated, but just in case
                if len(temp_gpr_list)==1:
                    logging.debug(temp_gpr_list)
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

    modelPrunedGPR = copy.deepcopy(modelPruned)
    return modelPrunedGPR

