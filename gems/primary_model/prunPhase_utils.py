
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import copy
import logging
from cobra.flux_analysis import single_reaction_deletion


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
	if 'Transport' not in rxn.name and 'transport' not in rxn.name \
            and 'Exchange' not in rxn.name and 'exchange' not in rxn.name:
            booleanList = makeBooleanFormat(
                    options.temp_target_BBH_dict,
                    options.tempModel_biggRxnid_locusTag_dict[biggRxnid])
            rxnToRemove_dict[biggRxnid] = calcBoolean(booleanList)

    options.rxnToRemove_dict = rxnToRemove_dict


def pruneModel(model, options):
    rxnToRemoveEssn_dict = {}
    rxnRemoved_dict = {}
    rxnRetained_dict = {}

    for rxnid in options.rxnToRemove_dict.keys():

        #Single reaction deletion is performed only for reactions labelled as "False"
        if options.rxnToRemove_dict[rxnid] == False:
            #Solver argument causes an error in cobrapy 0.5.8
            growth_rate_dict, solution_status_dict = single_reaction_deletion(
                    model, reaction_list=list([rxnid]), method='fba')

            #Check optimality first.
            if str(solution_status_dict.values()[0]) == 'optimal':

                #Full list of reactions and predicted growth rates upon their deletions
                rxnToRemoveEssn_dict[rxnid] = float(growth_rate_dict.values()[0])

                #Check growth rate upon reaction deletion
                if float(growth_rate_dict.values()[0]) >= float(options.cobrapy.non_zero_flux_cutoff):
                    model.remove_reactions(rxnid)
                    #List of reactions removed from the template model
                    rxnRemoved_dict[rxnid] = float(growth_rate_dict.values()[0])
                    logging.debug("Removed reaction: %s; %s; %s; %s"
                            %(rxnid, growth_rate_dict.values()[0],
                            len(model.reactions), len(model.metabolites)))
                else:
                    #List of reactions retained in the template model
                    rxnRetained_dict[rxnid] = float(growth_rate_dict.values()[0])
                    logging.debug("Retained reaction: %s; %s; %s; %s"
                            %(rxnid, growth_rate_dict.values()[0],
                            len(model.reactions), len(model.metabolites)))

    modelPruned = copy.deepcopy(model)

    #rxnToRemoveEssn_dict, rxnRemoved_dict and rxnRetained_dict:
    #Not used in the downstream of this pipeline
    return modelPruned


# Retain original GPR associations in the template model.
# Thus, parenthesis structures and Boolean relationships are reteined.
def get_gpr_fromString_toList(gpr_str):
    gpr_list3 = []
    gpr_list4 = []
    gpr_list6 = []

    # Remove parentheses
    gpr_list1 = gpr_str.split('(')
    for gpr in gpr_list1:
        if ')' in gpr:
            gpr_list2 = gpr.split(')')
            for gpr2 in gpr_list2:
                if len(gpr2) != 0:
                    gpr3 = gpr2.strip()
                    gpr_list3.append(gpr3)
        elif len(gpr) != 0:
            gpr2 = gpr.strip()
            gpr_list3.append(gpr2)
    print 'check1', gpr_list3

    # Remove 'and' and 'or'
    # TODO: May hold this information to retain original GPR structure
    for gpr in gpr_list3:
        if 'and' in gpr:
            gpr_list4 = gpr.split('and')
            # Check possibility of sub-parenthesis
            gpr_list5 = []
            for gpr2 in gpr_list4:
                if len(gpr2) != 0:
                    gpr3 = gpr2.strip()
                    gpr_list5.append(gpr3)
            if len(gpr_list5) != 0:
                gpr_list6.append(gpr_list5)
        elif 'or' in gpr:
            gpr_list4 = gpr.split('or')
            # Check possibility of sub-parenthesis
            gpr_list5 = []
            for gpr2 in gpr_list4:
                if len(gpr2) != 0:
                    gpr3 = gpr2.strip()
                    gpr_list5.append(gpr3)
            if len(gpr_list5) != 0:
                gpr_list6.append(gpr_list5)

    print 'check2', gpr_list6
    return gpr_list6


def swap_locusTag_tempModel(modelPruned, options):

    #Retrieve reactions associated with each homologous gene in template model
    for BBHrxn in modelPruned.reactions:
	booleanList = []
        #Retrieve all the genes associated with a reaction having the homologous gene
        #and transforms String to List
	booleanList = get_gpr_fromString_toList(BBHrxn.gene_reaction_rule)

        modified_booleanList = []
	for tempLocusTag in booleanList:

            #Check if the element itself is List.
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


        #Convert GPR in List to String:
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
	BBHrxn.gene_reaction_rule = booleanstring2

    modelPrunedGPR = copy.deepcopy(modelPruned)
    return modelPrunedGPR

