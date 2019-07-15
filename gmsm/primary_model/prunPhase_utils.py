
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
                locustag_candidate_list = \
                    copy.deepcopy(sorted(homology_ns.temp_target_BBH_dict[locustag]))
                if not 'and' in rxn.gene_reaction_rule:
                    target_locustag_list = []
                    for target_locustag in locustag_candidate_list:
                        if not target_locustag in rxn.gene_reaction_rule \
                            or target_locustag == locustag:
                            target_locustag_list.append(target_locustag)
                    if target_locustag_list:
                        homologs = ' or '.join(target_locustag_list)
                        new_gpr = rxn.gene_reaction_rule.replace(locustag, '%s' %homologs)
                    elif not locustag in locustag_candidate_list:
                        if (locustag + ' or ') in rxn.gene_reaction_rule:
                            new_gpr = rxn.gene_reaction_rule.replace((locustag + ' or '), '')
                        elif (' or ' + locustag) in rxn.gene_reaction_rule:
                            new_gpr = rxn.gene_reaction_rule.replace((' or ' + locustag), '')
                    else:
                        new_gpr = rxn.gene_reaction_rule
                    rxn.gene_reaction_rule = new_gpr
                elif not 'or' in rxn.gene_reaction_rule:
                    homologs = ' or '.join(locustag_candidate_list)
                    if len(locustag_candidate_list) == 1:
                        new_gpr = rxn.gene_reaction_rule.replace(locustag, '%s' %homologs)
                    else:
                        new_gpr = rxn.gene_reaction_rule.replace(locustag, '( %s )' %homologs)
                    rxn.gene_reaction_rule = new_gpr
                else:
                    gpr = copy.deepcopy(rxn.gene_reaction_rule)
                    new_gpr = ''
                    letter_loc = 0
                    locustag_loc = 0
                    new_locustag_loc = 0
                    parentheses_loc_dict = {}
                    extracted_gpr_list = []
                    changed_gpr_list = []
                    for letter in gpr:
                        if letter == '(':
                            parentheses_loc_dict['('] = letter_loc
                        if letter == ')':
                            parentheses_loc_dict[')'] = letter_loc
                            if '(' in parentheses_loc_dict:
                                extracted_gpr = gpr[parentheses_loc_dict['(']:letter_loc+1]
                                extracted_gpr_list.append([parentheses_loc_dict['('],letter_loc+1,extracted_gpr])
                                del parentheses_loc_dict['(']
                        letter_loc += 1
                    while gpr.find(locustag,locustag_loc) != -1:
                        locustag_loc = gpr.find(locustag, locustag_loc)
                        for left_loc, right_loc, extracted_gpr in extracted_gpr_list:
                            if left_loc <= locustag_loc and locustag_loc <= right_loc and locustag in extracted_gpr:
                                if not 'and' in extracted_gpr:
                                    target_locustag_list = []
                                    for target_locustag in locustag_candidate_list:
                                        if not target_locustag in extracted_gpr or target_locustag == locustag:
                                            target_locustag_list.append(target_locustag)
                                    if target_locustag_list:
                                        homologs = ' or '.join(target_locustag_list)
                                        new_extracted_gpr = extracted_gpr.replace(locustag, '%s' %homologs)
                                    elif not locustag in locustag_candidate_list:
                                        if (locustag + ' or ') in rxn.gene_reaction_rule:
                                            new_extracted_gpr = extracted_gpr.replace((locustag + ' or '), '')
                                        elif (' or ' + locustag) in rxn.gene_reaction_rule:
                                            new_extracted_gpr = extracted_gpr.replace((' or ' + locustag), '')
                                    else:
                                        new_extracted_gpr = extracted_gpr
                                else:
                                    homologs = ' or '.join(locustag_candidate_list)
                                    if len(locustag_candidate_list) == 1:
                                        new_extracted_gpr = extracted_gpr.replace(locustag, '%s' %homologs)
                                    else:
                                        new_extracted_gpr = extracted_gpr.replace(locustag, '( %s )' %homologs)
                                loc_diff = right_loc - left_loc
                                changed_gpr_list.append([left_loc, right_loc, new_extracted_gpr, loc_diff])
                                new_locustag_loc = right_loc
                                break
                        if locustag_loc < new_locustag_loc:
                            locustag_loc = new_locustag_loc
                            continue
                        homologs = ' or '.join(locustag_candidate_list)
                        locustag_loc_end = locustag_loc + len(locustag)
                        if len(locustag_candidate_list) > 1:
                            if (locustag + ' and ') in gpr[locustag_loc:locustag_loc_end+5] \
                                or (' and ' + locustag) in gpr[locustag_loc-5:locustag_loc_end]:
                                    changed_gpr = gpr[locustag_loc:locustag_loc_end].replace(locustag, '( %s )' %homologs)
                            else:
                                right_loc_for_section = 0
                                left_loc_for_section = 0
                                for left_loc, right_loc, extracted_gpr in extracted_gpr_list:
                                    if locustag_loc > right_loc:
                                        right_loc_for_section = right_loc
                                    if locustag_loc < left_loc:
                                        left_loc_for_section = left_loc
                                        break
                                target_locustag_list = []
                                if right_loc_for_section > left_loc_for_section:
                                    left_loc_for_section = locustag_loc_end
                                for target_locustag in gpr[right_loc_for_section:left_loc_for_section]:
                                    if not target_locustag in rxn.gene_reaction_rule \
                                        or target_locustag == locustag:
                                        target_locustag_list.append(target_locustag)
                                if target_locustag_list:
                                    homologs = ' or '.join(target_locustag_list)
                                    changed_gpr = gpr[right_loc_for_section:left_loc_for_section].replace(locustag, '%s' %homologs)
                                elif not locustag in locustag_candidate_list:
                                    if (locustag + ' or ') in gpr[right_loc_for_section:left_loc_for_section]:
                                        changed_gpr = 4
                                    elif (' or ' + locustag) in gpr[right_loc_for_section:left_loc_for_section]:
                                        changed_gpr = -4
                                else:
                                    changed_gpr = gpr[locustag_loc:locustag_loc_end]
                        else:
                            changed_gpr = gpr[locustag_loc:locustag_loc_end].replace(locustag, '%s' %homologs)
                        if type(changed_gpr) == int:
                            if changed_gpr == 4:
                                changed_gpr = ''
                                locustag_loc -= 4
                            else:
                                changed_gpr = ''
                                locustag_loc_end +=4
                        loc_diff = locustag_loc_end - locustag_loc
                        changed_gpr_list.append([locustag_loc, locustag_loc_end, changed_gpr, loc_diff])
                        new_locustag_loc = locustag_loc_end
                        locustag_loc = new_locustag_loc
                    changed_gpr_list = sorted(changed_gpr_list)
                    for i in range(len(changed_gpr_list)):
                        new_gpr = gpr[0:changed_gpr_list[i][0]] + \
                                  changed_gpr_list[i][2] + gpr[changed_gpr_list[i][1]:]
                        for j in range(len(changed_gpr_list)):
                            changed_gpr_list[j][0] += len(changed_gpr_list[i][2])-changed_gpr_list[i][3]
                            changed_gpr_list[j][1] += len(changed_gpr_list[i][2])-changed_gpr_list[i][3]
                        gpr = new_gpr
                    rxn.gene_reaction_rule = new_gpr


    modelPrunedGPR = copy.deepcopy(modelPruned)
    return modelPrunedGPR