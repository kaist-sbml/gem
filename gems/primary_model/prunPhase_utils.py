
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


# Retain structures of original GPR associations in the template model:
# e.g., parenthesis structures and Boolean relationships
def swap_locustag_with_homolog(modelPruned, options):

    for locustag in options.temp_target_BBH_dict:
        for rxn in modelPruned.reactions:
            if rxn.gene_reaction_rule and locustag in rxn.gene_reaction_rule:
                if len(options.temp_target_BBH_dict[locustag]) == 1:
                    new_gpr = rxn.gene_reaction_rule.replace(
                            locustag, options.temp_target_BBH_dict[locustag][0])
                    rxn.gene_reaction_rule = new_gpr
                elif len(options.temp_target_BBH_dict[locustag]) > 1:
                    homologs = ' or '.join(options.temp_target_BBH_dict[locustag])
                    new_gpr = rxn.gene_reaction_rule.replace(locustag, '( %s )' %homologs)
                    rxn.gene_reaction_rule = new_gpr

    modelPrunedGPR = copy.deepcopy(modelPruned)
    return modelPrunedGPR
