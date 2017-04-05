
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import ast
import copy
import logging
import pyparsing
from cobra.flux_analysis import single_reaction_deletion


def get_rxn_fate2(locustag_list, temp_target_BBH_dict):
    for locustag in locustag_list:
        if isinstance(locustag, list):
            return get_rxn_fate2(locustag, temp_target_BBH_dict)

    bbh_avail_list = []
    boolop_list = []

    for locustag in locustag_list:
        if locustag == 'AND' or locustag == 'and' or locustag == 'OR' or locustag == 'or':
            boolop_list.append(locustag)
        elif locustag != '1' and locustag != '0':
            # Check BBH availability for a given gene
            # '1': gene with BBHs available
            # '0': gene with no BBHs (nonBBHs) available
            # For ( A and B), if one of them is a gene with nonBBHs,
            #its rxn should be False (subject to removal).
            if locustag in temp_target_BBH_dict:
                bbh_avail_list.append('1')
            else:
                bbh_avail_list.append('0')
        elif locustag == '1' and locustag != '0':
            bbh_avail_list.append('1')
        elif locustag != '1' and locustag == '0':
            bbh_avail_list.append('0')


    boolop_list2 = list(set(boolop_list))

    # e.g., bbh_avail_list = ['1']; no Boolean logic
    if len(boolop_list2) == 0:
        return locustag_list, bbh_avail_list[0]

    elif len(boolop_list2) == 1:
        if 'AND' in boolop_list2 or 'and' in boolop_list2:
            return locustag_list, min(bbh_avail_list)
        elif 'OR' in boolop_list2 or 'or' in boolop_list2:
            return locustag_list, max(bbh_avail_list)

    # gene_reaction_rule contains both AND and OR within a parenthesis
    #and without parentheses around AND. This is a bad practise.
    elif len(boolop_list2) >= 2:

        # NOTE: This issue has not been resolved. OR is returned regardless of the Boolean.
        #See: test_primary_model.py

        #print '1', bbh_avail_list
        #print '2', boolop_list2
        #bbh_avail_list_str = str(bbh_avail_list)
        #bbh_avail_list_str2 = bbh_avail_list_str.replace("'", "")
        #bbh_avail_list_str3 = bbh_avail_list_str2.replace(",", "")
        #bbh_avail_list_str4 = bbh_avail_list_str3.replace("[", "(")
        #bbh_avail_list_str5 = bbh_avail_list_str4.replace("]", ")")
        #print '3', bbh_avail_list_str5

        #gpr_regex = pyparsing.Word(pyparsing.alphanums)
        #and_booleanop = pyparsing.oneOf('AND and')
        #or_booleanop = pyparsing.oneOf('OR or')
        #expr = pyparsing.infixNotation(gpr_regex,
        #        [
        #            (and_booleanop, 2, pyparsing.opAssoc.LEFT),
        #            (or_booleanop, 2, pyparsing.opAssoc.LEFT)
        #        ])
        #bbh_avail_list = expr.parseString(bbh_avail_list_str5)[0].asList()
        #print '4', bbh_avail_list

        logging.warning("Bad 'gene_reaction_rule' description", locustag_list)
        return locustag_list, max(bbh_avail_list)


def get_rxn_fate(locustag_list, temp_target_BBH_dict):
    while True:
        locustag_list_str = str(locustag_list)
        locustag, rxn_fate = get_rxn_fate2(locustag_list, temp_target_BBH_dict)
        rxn_fate = "'" + str(rxn_fate) + "'"
        locustag_list_str = locustag_list_str.replace(str(locustag), str(rxn_fate))
        locustag_list = ast.literal_eval(locustag_list_str)
        if type(locustag_list) != list:
            rxn_fate2 = locustag_list
            break
    return rxn_fate2


def label_rxn_to_remove(model, options):
    rxnToRemove_dict = {}

    for biggRxnid in options.tempModel_biggRxnid_locusTag_dict:
	rxn = model.reactions.get_by_id(biggRxnid)
        #Prevent removal of transport reactions from the template model
	if 'Transport' not in rxn.name and 'transport' not in rxn.name \
            and 'Exchange' not in rxn.name and 'exchange' not in rxn.name:
            rxnToRemove_dict[biggRxnid] = get_rxn_fate(
                    options.tempModel_biggRxnid_locusTag_dict[biggRxnid],
                    options.temp_target_BBH_dict)
            logging.debug('%s; %s', rxn.id, rxnToRemove_dict[biggRxnid])

    options.rxnToRemove_dict = rxnToRemove_dict


def prune_model(model, options):

    for rxnid in options.rxnToRemove_dict:

        logging.debug("Fate of reaction %s: %s", rxnid, options.rxnToRemove_dict[rxnid])
        #Single reaction deletion is performed only for reactions labelled as "False"
        if options.rxnToRemove_dict[rxnid] == '0':

            #Solver argument causes an error in cobrapy 0.5.8
            growth_rate_dict, solution_status_dict = single_reaction_deletion(
                    model, reaction_list=list([rxnid]), method='fba')

            #Check optimality first.
            if str(solution_status_dict.values()[0]) == 'optimal':

                #Check growth rate upon reaction deletion
                if float(growth_rate_dict.values()[0]) >= float(options.cobrapy.non_zero_flux_cutoff):
                    model.remove_reactions(rxnid)
                    logging.debug("Removed reaction: %s; %s; %s; %s"
                            %(rxnid, growth_rate_dict.values()[0],
                            len(model.reactions), len(model.metabolites)))
                else:
                    logging.debug("Retained reaction: %s; %s; %s; %s"
                            %(rxnid, growth_rate_dict.values()[0],
                            len(model.reactions), len(model.metabolites)))

    modelPruned = copy.deepcopy(model)

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

