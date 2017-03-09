
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import cobra
import copy
import logging
import os
import pickle
import re
import urllib2
from os.path import isdir, join, abspath, dirname
from gems.util import time_bomb, stabilize_model

#Retrieves a list of reaction IDs using their EC numbers from KEGG
#Input: E.C number in string form (e.g., 4.1.3.6)
#Output: reactionID in list form (e.g., ['R00362'])
def get_rxnid_from_ECNumber(enzymeEC, options):
    url = options.urls.kegg_enzyme + '%s' %enzymeEC
    ecinfo_text = urllib2.urlopen(url).read()

    #Original line also extracted genes in other organisms: R50912; R50345 (NOT rxnid)
    #The HTTP error was solved by putting "\\b" only at the end (not at the front)
    #in order to also include reaction ID followed by "other" in KEGG
    rxnid_set = re.findall(r'\s+R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]'+'\\b', ecinfo_text)
    rxnid_list = []
    for each_set in rxnid_set:
        rxnid = re.findall('R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]+', each_set)
        rxnid_list+=rxnid

        #Removes redundancy
        rxnid_list = list(set(rxnid_list))
    return rxnid_list


#Get reaction information using its ID from KEGG
#Input: KEGG rxnid in string form (e.g., R00362)
#Output: Reaction information for 'Name', 'Definition', and 'Equation' as dictionary form
#{'NAME': 'citrate oxaloacetate-lyase', 'DEFINITION': Citrate <=> Acetate + Oxaloacetate,
#'EQUATION': C00158 <=> C00033 + C00036}
def get_rxnInfo_from_rxnid(rxnid, options):
    url = options.urls.kegg_rn + '%s' %rxnid
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
        #Consider only reactions mapped in pathways
        #Otherwise reactions have unspecified molecules having R groups
        if sptlist[0].strip() == 'PATHWAY':
            PATHWAY = ' '.join(sptlist[1:])
            return {'NAME':NAME, 'DEFINITION':DEFINITION,
                    'EQUATION':EQUATION, 'ENZYME':ENZYME, 'PATHWAY':PATHWAY}


def load_cache(cache_dir, cache_data, options):

    if os.path.isfile(cache_dir):
        # For regular update of the cache
        time_bomb(cache_dir, options)

        try:
            with open(cache_dir, 'rb') as f:
                cache_data = pickle.load(f)
                return cache_data
        except pickle.UnpicklingError as e:
            logging.warning("Could not read '%s': %s", cache_dir, e)
            if '_dict' in cache_data:
                cache_data = {}
            elif '_list' in cache_data:
                cache_data = []
            return cache_data
        except IOError as e:
            logging.warning("Can't open %s in 'rb' mode: %s", cache_dir, e)
            if '_dict' in cache_data:
                cache_data = {}
            elif '_list' in cache_data:
                cache_data = []
            return cache_data
    else:
        logging.debug('No cache exists: %s', cache_dir)
        if '_dict' in cache_data:
            cache_data = {}
        elif '_list' in cache_data:
            cache_data = []
        return cache_data


def create_cache(cache_dir, cache_data):
    try:
        with open(cache_dir, 'wb') as f:
            pickle.dump(cache_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    except pickle.PicklingError as e:
        logging.warning("Error in serializing '%s': %s", cache_data, e)
    except IOError as e:
        logging.warning("Can't open %s in 'wb' mode: %s", cache_data, e)


def get_targetGenome_locusTag_ec_nonBBH_dict(options):
    targetGenome_locusTag_ec_nonBBH_dict = {}

    for locusTag in options.nonBBH_list:
	if locusTag in options.targetGenome_locusTag_ec_dict.keys():
            targetGenome_locusTag_ec_nonBBH_dict[locusTag] = \
            options.targetGenome_locusTag_ec_dict[locusTag]
    options.targetGenome_locusTag_ec_nonBBH_dict = targetGenome_locusTag_ec_nonBBH_dict


def edit_mnxr_kegg_dict(keggid, options):
    for mnxr in options.mnxr_kegg_dict.keys():
        # Remove candidate KEGG rxn IDs from consideration
        if keggid in options.mnxr_kegg_dict[mnxr]:
            cnt = len(options.mnxr_kegg_dict[mnxr])
            options.mnxr_kegg_dict[mnxr].remove(keggid)
            logging.debug('Number of KEGG IDs for %s: %s --> %s',
                          mnxr, cnt, len(options.mnxr_kegg_dict[mnxr]))

            if len(options.mnxr_kegg_dict[mnxr]) == 0:
                del options.mnxr_kegg_dict[mnxr]
                logging.debug('%s removed', mnxr)


def get_rxnid_locusTag_dict(rxnid_locusTag_dict, rxnid, locusTag):

    # Create 'rxnid_locusTag_dict'
    if rxnid not in rxnid_locusTag_dict:
        rxnid_locusTag_dict[rxnid] = [locusTag]
    elif rxnid in rxnid_locusTag_dict.keys():
        rxnid_locusTag_dict[rxnid].append(locusTag)

    return rxnid_locusTag_dict


def get_rxnid_info_dict_from_kegg(options):
    cache_ec_rxn_dict = {}
    cache_rxnid_info_dict = {}
    rxnid_info_dict = {}
    rxnid_locusTag_dict = {}

    cache_dumped_ec_list = []
    cache_dumped_rxnid_list = []

    # Folder for cache files
    primary_model_dir = join(dirname(abspath(__file__)))
    kegg_cache_dir = join(primary_model_dir, 'kegg_cache')

    if not isdir(kegg_cache_dir):
        os.makedirs(kegg_cache_dir)

    # Access KEGG cache files
    cache_ec_rxn_dict_dir = join(kegg_cache_dir, 'cache_ec_rxn_dict.p')
    cache_rxnid_info_dict_dir = join(kegg_cache_dir, 'cache_rxnid_info_dict.p')
    cache_dumped_ec_list_dir = join(kegg_cache_dir, 'cache_dumped_ec_list.p')
    cache_dumped_rxnid_list_dir = join(kegg_cache_dir, 'cache_dumped_rxnid_list.p')

    cache_ec_rxn_dict = load_cache(
            cache_ec_rxn_dict_dir, cache_ec_rxn_dict, options)
    cache_rxnid_info_dict = load_cache(
            cache_rxnid_info_dict_dir, cache_rxnid_info_dict, options)
    cache_dumped_ec_list = load_cache(
            cache_dumped_ec_list_dir, cache_dumped_ec_list, options)
    cache_dumped_rxnid_list = load_cache(
            cache_dumped_rxnid_list_dir, cache_dumped_rxnid_list, options)

    for locusTag in options.targetGenome_locusTag_ec_nonBBH_dict.keys():
	for enzymeEC in options.targetGenome_locusTag_ec_nonBBH_dict[locusTag]:

            #KEGG REST does not accept unspecific EC_number: e.g., 3.2.2.-
            if '-' not in enzymeEC:

                #Check cache file
                if enzymeEC in cache_ec_rxn_dict:
                    if cache_ec_rxn_dict[enzymeEC]:
                        rxnid_list = cache_ec_rxn_dict[enzymeEC]
                        logging.debug("EC_number info from a cache file: %s, %s",
                                      locusTag, enzymeEC)
                else:
                    if enzymeEC not in cache_dumped_ec_list:
                        rxnid_list = get_rxnid_from_ECNumber(enzymeEC, options)
                        logging.debug("EC_number info fetched from KEGG: %s, %s",
                                    locusTag, enzymeEC)
                        if rxnid_list:
                            # Store new info in a cache
                            if enzymeEC not in cache_ec_rxn_dict:
                                cache_ec_rxn_dict[enzymeEC] = rxnid_list
                        elif not rxnid_list and enzymeEC not in cache_dumped_ec_list:
                            cache_dumped_ec_list.append(enzymeEC)

                for rxnid in rxnid_list:

                    #Check cache file
                    if rxnid in cache_rxnid_info_dict:
                        rxnid_info_dict[rxnid] = cache_rxnid_info_dict[rxnid]
                        logging.debug('Reaction info from a cache file: %s, %s, %s',
                                      locusTag, enzymeEC, rxnid)
                    else:
                        if rxnid not in cache_dumped_rxnid_list:
                            rxnid_info = get_rxnInfo_from_rxnid(rxnid, options)
                            logging.debug('Reaction info fetched from KEGG: %s, %s, %s',
                                      locusTag, enzymeEC, rxnid)

                            # Create 'rxnid_info_dict'
                            if rxnid_info != None:
                                rxnid_info_dict[rxnid] = rxnid_info
                                # Store new info in a cache
                                if rxnid not in cache_rxnid_info_dict:
                                    cache_rxnid_info_dict[rxnid] = rxnid_info
                            else:
                                logging.debug('No reaction info available for %s, %s, %s',
                                          locusTag, enzymeEC, rxnid)
                                if rxnid not in cache_dumped_rxnid_list:
                                    cache_dumped_rxnid_list.append(rxnid)

                                edit_mnxr_kegg_dict(rxnid, options)
                        else:
                            edit_mnxr_kegg_dict(rxnid, options)
                    rxnid_locusTag_dict = get_rxnid_locusTag_dict(
                            rxnid_locusTag_dict, rxnid, locusTag)
            else:
                logging.debug("EC_number NOT submitted to KEGG: %s, %s",
                              locusTag, enzymeEC)

    create_cache(cache_ec_rxn_dict_dir, cache_ec_rxn_dict)
    create_cache(cache_rxnid_info_dict_dir, cache_rxnid_info_dict)
    create_cache(cache_dumped_ec_list_dir, cache_dumped_ec_list)
    create_cache(cache_dumped_rxnid_list_dir, cache_dumped_rxnid_list)

    options.rxnid_info_dict = rxnid_info_dict
    options.rxnid_locusTag_dict = rxnid_locusTag_dict


def get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, options):
    modelPrunedGPR_mnxr_list = []

    for j in range(len(modelPrunedGPR.reactions)):
        biggid = modelPrunedGPR.reactions[j].id

        if biggid in options.bigg_mnxr_dict:
            mnxr = options.bigg_mnxr_dict[biggid]
            if mnxr not in modelPrunedGPR_mnxr_list:
                modelPrunedGPR_mnxr_list.append(mnxr)
        else:
            logging.debug('BiGG reaction %s NOT in MNXref', biggid)

    options.modelPrunedGPR_mnxr_list = modelPrunedGPR_mnxr_list


def get_mnxr_to_add_list(options):

    mnxr_to_add_list = []
    for rxnid in options.rxnid_info_dict:
        for mnxr, kegg_list in options.mnxr_kegg_dict.iteritems():
            if rxnid in kegg_list:
                # Check reaction duplicates
                if mnxr not in options.modelPrunedGPR_mnxr_list:
                    if mnxr not in mnxr_to_add_list:
                        if mnxr in options.mnxref.reactions:
                            mnxr_to_add_list.append(mnxr)
                        else:
                            logging.debug('%s (%s) not available in MNXref', rxnid, mnxr)
                else:
                    logging.debug('%s (%s) already in the model', rxnid, mnxr)

    logging.debug('Number of KEGG reactions to be added: %s', len(mnxr_to_add_list))
    options.mnxr_to_add_list = mnxr_to_add_list


def add_nonBBH_rxn(modelPrunedGPR, options):

    for mnxr in options.mnxr_to_add_list:

        # Choose KEGG reaction ID with a greater value for multiple KEGG IDs given to MNXR
        if len(options.mnxr_kegg_dict[mnxr]) > 1:
            keggid_list = []
            keggid_list = options.mnxr_kegg_dict[mnxr]
            keggid_list.sort()
            kegg_id = keggid_list[-1]
        elif len(options.mnxr_kegg_dict[mnxr]) == 1:
            kegg_id = options.mnxr_kegg_dict[mnxr][0]

        logging.debug("--------------------")
        logging.debug("Reaction to be added: %s; %s", mnxr, kegg_id)

        rxn = options.mnxref.reactions.get_by_id(mnxr)
        modelPrunedGPR.add_reaction(rxn)

        #Re-define ID
        rxn.id = kegg_id

        #Re-define Name
        rxn.name = options.rxnid_info_dict[kegg_id]['NAME']

        #'add_reaction' requires writing/reloading of the model
        modelPrunedGPR = stabilize_model(modelPrunedGPR, kegg_id, options)

        rxn = modelPrunedGPR.reactions.get_by_id(kegg_id)

        # NOTE:
        # Writing the reaction notes (e.g., 'GENE ASSOCIATION' and 'SUBSYSTEM')
        #seems to be unstable.
        # The modified reaction notes are missing in resulting SBML file.
        # This bug appears to happen if the new reaction comes from another SBML file.
        # Addition of a newly built reaction using Metabolite object does not seem to have
        #this potential bug.

        #GPR association
        if len(options.rxnid_locusTag_dict[kegg_id]) == 1:
            gpr = '( %s )' %(options.rxnid_locusTag_dict[kegg_id][0])
        else:
            count = 1
            for locusTag in options.rxnid_locusTag_dict[kegg_id]:

                #Check whether the submitted gbk file contains "/product" for CDS
                if locusTag in options.targetGenome_locusTag_prod_dict:

                    #Consider "and" relationship in the GPR association
                    if 'subunit' in options.targetGenome_locusTag_prod_dict[locusTag]:
                        count += 1
            if count == len(options.rxnid_locusTag_dict[kegg_id]):
                gpr = ' and '.join(options.rxnid_locusTag_dict[kegg_id])
            else:
                gpr = ' or '.join(options.rxnid_locusTag_dict[kegg_id])
            gpr = '( %s )' %(gpr)

        rxn.gene_reaction_rule = gpr

        # Subsystem
        rxn.subsystem = options.rxnid_info_dict[kegg_id]['PATHWAY']

        modelPrunedGPR = stabilize_model(modelPrunedGPR, kegg_id, options)

        logging.debug("Number of reactions in the model: %s",
                len(modelPrunedGPR.reactions))

        #Check model prediction consistency
        target_exrxnid_flux_dict = get_exrxnid_flux(modelPrunedGPR,
                options.template_exrxnid_flux_dict)
        exrxn_flux_change_list = check_exrxn_flux_direction(
                options.template_exrxnid_flux_dict,
                target_exrxnid_flux_dict, options)

        if 'F' in exrxn_flux_change_list:
            #'remove_reactions' does not seem to require
            #writing/reloading of the model
            modelPrunedGPR.remove_reactions(rxn)

    target_model = copy.deepcopy(modelPrunedGPR)
    return target_model


#Output: a dictionary file for major Exchange reactions {Exchange reaction ID:flux value}
def get_exrxnid_flux(model, template_exrxnid_flux_dict):

    target_exrxnid_flux_dict = {}
    model.optimize()

    for exrxn_id in template_exrxnid_flux_dict.keys():
        if exrxn_id in model.solution.x_dict:
            target_exrxnid_flux_dict[exrxn_id] = model.solution.x_dict[exrxn_id]
        else:
            continue
    return target_exrxnid_flux_dict


#Output: a list file having either T or F for major Exchange reactions
def check_exrxn_flux_direction(
        template_exrxnid_flux_dict, target_exrxnid_flux_dict, options):

    exrxn_flux_change_list = []

    for exrxn_id in template_exrxnid_flux_dict.keys():
        if exrxn_id in target_exrxnid_flux_dict.keys():
            template_exrxn_flux = template_exrxnid_flux_dict[exrxn_id]
            target_exrxn_flux = target_exrxnid_flux_dict[exrxn_id]
            ratio_exrxn_flux = float(target_exrxn_flux)/float(template_exrxn_flux)

            #Similar species are allowed to uptake nutrients within a decent range
            if float(target_exrxn_flux)*float(template_exrxn_flux) > float(options.cobrapy.non_zero_flux_cutoff) \
                    and float(options.cobrapy.nutrient_uptake_rate) < ratio_exrxn_flux \
                    and ratio_exrxn_flux < float(options.cobrapy.nutrient_uptake_rate):
                exrxn_flux_change_list.append('T')

            #Cause drastic changes in Exchange reaction fluxes
            #(direction and/or magnitude)
            else:
                exrxn_flux_change_list.append('F')

    return exrxn_flux_change_list
