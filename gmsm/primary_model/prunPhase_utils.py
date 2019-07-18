
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

    elif len(boolop_list2) >= 2:
        bbh_avail_list_str = str(locustag_list)
        bbh_avail_list_str = bbh_avail_list_str.replace("'", "")
        bbh_avail_list_str = bbh_avail_list_str.replace(",", "")
        bbh_avail_list_str = bbh_avail_list_str.replace("[", "(")
        bbh_avail_list_str = bbh_avail_list_str.replace("]", ")")
        gpr_regex = pyparsing.Word(pyparsing.alphanums + '_' + '.')
        and_booleanop = pyparsing.oneOf('AND and')
        or_booleanop = pyparsing.oneOf('OR or')
        expr = pyparsing.infixNotation(gpr_regex,
                [
                    (and_booleanop, 2, pyparsing.opAssoc.LEFT),
                    (or_booleanop, 2, pyparsing.opAssoc.LEFT)
                ])
        bbh_avail_list = expr.parseString(bbh_avail_list_str)[0].asList()
        return locustag_list, bbh_avail_list


def get_rxn_fate(locustag_list, temp_target_BBH_dict):
    while True:
        locustag_list_str = str(locustag_list)
        locustag, rxn_fate = get_rxn_fate2(locustag_list, temp_target_BBH_dict)
        if isinstance(rxn_fate, str):
            rxn_fate = "'" + str(rxn_fate) + "'"
        locustag_list_str = locustag_list_str.replace(str(locustag), str(rxn_fate))
        locustag_list = ast.literal_eval(locustag_list_str)
        if type(locustag_list) != list:
            rxn_fate2 = locustag_list
            break
    return rxn_fate2


def label_rxn_to_remove(model, io_ns, homology_ns, primary_model_ns):
    rxnToRemove_dict = {}

    for biggRxnid in io_ns.tempModel_biggRxnid_locusTag_dict:
	rxn = model.reactions.get_by_id(biggRxnid)
        #Prevent removal of transport reactions from the template model
	if 'Transport' not in rxn.name and 'transport' not in rxn.name \
            and 'Exchange' not in rxn.name and 'exchange' not in rxn.name:
            rxnToRemove_dict[biggRxnid] = get_rxn_fate(
                    io_ns.tempModel_biggRxnid_locusTag_dict[biggRxnid],
                    homology_ns.temp_target_BBH_dict)

    primary_model_ns.rxnToRemove_dict = rxnToRemove_dict


def prune_model(model, config_ns, primary_model_ns):

    for rxnid in primary_model_ns.rxnToRemove_dict:

        logging.debug("Fate of reaction %s: %s", rxnid, primary_model_ns.rxnToRemove_dict[rxnid])
        #Single reaction deletion is performed only for reactions labelled as "False"
        if primary_model_ns.rxnToRemove_dict[rxnid] == '0':

            flux_dist = single_reaction_deletion(
                        model, reaction_list=list([rxnid]), method='fba')

            #Check optimality first.
            if flux_dist.status[rxnid] == 'optimal':

                #Check growth rate upon reaction deletion
                if float(flux_dist.flux[rxnid]) >= \
                        float(config_ns.cobrapy.non_zero_flux_cutoff):
                    model.remove_reactions(rxnid)
                    logging.debug("Removed reaction: %s; %s; %s; %s"
                            %(rxnid, flux_dist.flux[rxnid],
                            len(model.reactions), len(model.metabolites)))
                else:
                    logging.debug("Retained reaction: %s; %s; %s; %s"
                            %(rxnid, flux_dist.flux[rxnid],
                            len(model.reactions), len(model.metabolites)))
            else:
                logging.debug("Reaction not optimal: %s; %s",
                              rxnid, flux_dist.status[rxnid])

    modelPruned = copy.deepcopy(model)

    return modelPruned


# Retain structures of original GPR associations in the template model:
# e.g., parenthesis structures and Boolean relationships
def swap_locustag_with_homolog(modelPruned, homology_ns):

    for locustag in homology_ns.temp_target_BBH_dict:
        for rxn in modelPruned.reactions:
            if rxn.gene_reaction_rule and locustag in rxn.gene_reaction_rule:
                if len(homology_ns.temp_target_BBH_dict[locustag]) == 1:
                    new_gpr = rxn.gene_reaction_rule.replace(
                            locustag, homology_ns.temp_target_BBH_dict[locustag][0])
                    rxn.gene_reaction_rule = new_gpr
                elif len(homology_ns.temp_target_BBH_dict[locustag]) > 1:
                    homologs = ' or '.join(homology_ns.temp_target_BBH_dict[locustag])
                    new_gpr = rxn.gene_reaction_rule.replace(locustag, '( %s )' %homologs)
                    rxn.gene_reaction_rule = new_gpr

    modelPrunedGPR = copy.deepcopy(modelPruned)
    return modelPrunedGPR

