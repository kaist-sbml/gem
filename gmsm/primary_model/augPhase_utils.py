
import cobra
import copy
import logging
import os
import pickle
import re
import urllib2
from cobra import Metabolite, Reaction
from os.path import isdir, join, abspath, dirname
from gmsm import utils

#Retrieves a list of reaction IDs using their EC numbers from KEGG
#Input: E.C number in string form (e.g., 4.1.3.6)
#Output: reactionID in list form (e.g., ['R00362'])
def get_rxnid_from_ECNumber(rxnid_list, enzymeEC, config_ns):
    url = config_ns.urls.kegg_enzyme + '%s' %enzymeEC
    ecinfo_text = urllib2.urlopen(url).read()

    #Original line also extracted genes in other organisms: R50912; R50345 (NOT rxnid)
    #The HTTP error was solved by putting "\\b" only at the end (not at the front)
    #in order to also include reaction ID followed by "other" in KEGG
    rxnid_set = re.findall(r'\s+R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]'+'\\b', ecinfo_text)

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
def get_rxnInfo_from_rxnid(rxnid, config_ns):
    url = config_ns.urls.kegg_rn + '%s' %rxnid
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


def load_cache(cache_dir, cache_data):

    if os.path.isfile(cache_dir):
        # NOTE: Disabled at the moment because Docker images do not have cache files.
        # TODO: As to the motivation behind this function, automatic removal of
        #existing cache files everytime GMSM starts or ends can be an option.
        # For regular update of the cache
        #utils.time_bomb(cache_dir, config_ns)

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


def get_targetGenome_locusTag_ec_nonBBH_dict(io_ns, homology_ns, primary_model_ns):
    targetGenome_locusTag_ec_nonBBH_dict = {}

    for locusTag in homology_ns.nonBBH_list:
	if locusTag in io_ns.targetGenome_locusTag_ec_dict.keys():
            targetGenome_locusTag_ec_nonBBH_dict[locusTag] = \
            io_ns.targetGenome_locusTag_ec_dict[locusTag]
    primary_model_ns.targetGenome_locusTag_ec_nonBBH_dict = targetGenome_locusTag_ec_nonBBH_dict


def edit_mnxr_kegg_dict(keggid, io_ns):
    for mnxr in io_ns.mnxr_kegg_dict.keys():
        # Remove candidate KEGG rxn IDs from consideration
        if keggid in io_ns.mnxr_kegg_dict[mnxr]:
            cnt = len(io_ns.mnxr_kegg_dict[mnxr])
            io_ns.mnxr_kegg_dict[mnxr].remove(keggid)
            logging.debug('Number of KEGG IDs for %s: %s --> %s',
                          mnxr, cnt, len(io_ns.mnxr_kegg_dict[mnxr]))

            if len(io_ns.mnxr_kegg_dict[mnxr]) == 0:
                del io_ns.mnxr_kegg_dict[mnxr]
                logging.debug('%s removed', mnxr)


def get_rxnid_locusTag_dict(rxnid_locusTag_dict, rxnid, locusTag):

    # Create 'rxnid_locusTag_dict'
    if rxnid not in rxnid_locusTag_dict:
        rxnid_locusTag_dict[rxnid] = [locusTag]
    elif rxnid in rxnid_locusTag_dict.keys():
        rxnid_locusTag_dict[rxnid].append(locusTag)

    return rxnid_locusTag_dict


def get_rxnid_info_dict_from_kegg(io_ns, config_ns, primary_model_ns):
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
            cache_ec_rxn_dict_dir, cache_ec_rxn_dict)
    cache_rxnid_info_dict = load_cache(
            cache_rxnid_info_dict_dir, cache_rxnid_info_dict)
    cache_dumped_ec_list = load_cache(
            cache_dumped_ec_list_dir, cache_dumped_ec_list)
    cache_dumped_rxnid_list = load_cache(
            cache_dumped_rxnid_list_dir, cache_dumped_rxnid_list)

    for locusTag in primary_model_ns.targetGenome_locusTag_ec_nonBBH_dict.keys():
	for enzymeEC in primary_model_ns.targetGenome_locusTag_ec_nonBBH_dict[locusTag]:
            # This should be declared in case enzymeEC is not available at KEGG: e.g.,
            #"UnboundLocalError: local variable 'rxnid_list' referenced before assignment"
            rxnid_list = []

            #KEGG REST does not accept unspecific EC_number: e.g., 3.2.2.-
            if '-' not in enzymeEC:

                #Check cache file
                if enzymeEC in cache_ec_rxn_dict:
                    if cache_ec_rxn_dict[enzymeEC]:
                        rxnid_list = cache_ec_rxn_dict[enzymeEC]
                        logging.debug("EC_number info from a cache file: %s, %s",
                                      locusTag, enzymeEC)

                elif enzymeEC not in cache_dumped_ec_list:
                    rxnid_list = get_rxnid_from_ECNumber(rxnid_list, enzymeEC, config_ns)
                    logging.debug("EC_number info fetched from KEGG: %s, %s",
                                  locusTag, enzymeEC)
                    if rxnid_list:
                        # Store new info in a cache
                        if enzymeEC not in cache_ec_rxn_dict:
                            cache_ec_rxn_dict[enzymeEC] = rxnid_list
                    elif not rxnid_list and enzymeEC not in cache_dumped_ec_list:
                        cache_dumped_ec_list.append(enzymeEC)

                else:
                    logging.debug("EC_number NOT available at KEGG: %s, %s",
                                  locusTag, enzymeEC)

                for rxnid in rxnid_list:

                    #Check cache file
                    if rxnid in cache_rxnid_info_dict:
                        rxnid_info_dict[rxnid] = cache_rxnid_info_dict[rxnid]
                        logging.debug('Reaction info from a cache file: %s, %s, %s',
                                          locusTag, enzymeEC, rxnid)
                    else:
                        if rxnid not in cache_dumped_rxnid_list:
                            rxnid_info = get_rxnInfo_from_rxnid(rxnid, config_ns)
                            logging.debug(
                                         'Reaction info fetched from KEGG: %s, %s, %s',
                                         locusTag, enzymeEC, rxnid)

                            # Create 'rxnid_info_dict'
                            if rxnid_info != None:
                                rxnid_info_dict[rxnid] = rxnid_info
                                # Store new info in a cache
                                if rxnid not in cache_rxnid_info_dict:
                                    cache_rxnid_info_dict[rxnid] = rxnid_info
                            else:
                                logging.debug(
                                            'No reaction info available for %s, %s, %s',
                                            locusTag, enzymeEC, rxnid)
                                if rxnid not in cache_dumped_rxnid_list:
                                    cache_dumped_rxnid_list.append(rxnid)

                                edit_mnxr_kegg_dict(rxnid, io_ns)
                        else:
                            edit_mnxr_kegg_dict(rxnid, io_ns)

                    rxnid_locusTag_dict = get_rxnid_locusTag_dict(
                                rxnid_locusTag_dict, rxnid, locusTag)

    create_cache(cache_ec_rxn_dict_dir, cache_ec_rxn_dict)
    create_cache(cache_rxnid_info_dict_dir, cache_rxnid_info_dict)
    create_cache(cache_dumped_ec_list_dir, cache_dumped_ec_list)
    create_cache(cache_dumped_rxnid_list_dir, cache_dumped_rxnid_list)

    primary_model_ns.rxnid_info_dict = rxnid_info_dict
    primary_model_ns.rxnid_locusTag_dict = rxnid_locusTag_dict


def get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, io_ns, primary_model_ns):
    modelPrunedGPR_mnxr_list = []

    for j in range(len(modelPrunedGPR.reactions)):
        biggid = modelPrunedGPR.reactions[j].id

        if biggid in io_ns.bigg_mnxr_dict:
            mnxr = io_ns.bigg_mnxr_dict[biggid]
            if mnxr not in modelPrunedGPR_mnxr_list:
                modelPrunedGPR_mnxr_list.append(mnxr)
        else:
            logging.debug('BiGG reaction %s NOT in MNXref', biggid)

    primary_model_ns.modelPrunedGPR_mnxr_list = modelPrunedGPR_mnxr_list


def get_mnxr_to_add_list(io_ns, primary_model_ns):

    mnxr_to_add_list = []
    for rxnid in primary_model_ns.rxnid_info_dict:
        for mnxr, kegg_list in io_ns.mnxr_kegg_dict.iteritems():
            if rxnid in kegg_list:
                # Check reaction duplicates
                if mnxr not in primary_model_ns.modelPrunedGPR_mnxr_list:
                    if mnxr not in mnxr_to_add_list:
                        if mnxr in io_ns.mnxref.reactions:
                            mnxr_to_add_list.append(mnxr)
                        else:
                            logging.debug('%s (%s) not available in MNXref', rxnid, mnxr)
                else:
                    logging.debug('%s (%s) already in the model', rxnid, mnxr)

    logging.debug('Number of KEGG reactions to be added (based on EC numbers): %s',
                  len(mnxr_to_add_list))
    primary_model_ns.mnxr_to_add_list = mnxr_to_add_list


def add_nonBBH_rxn(modelPrunedGPR, io_ns, config_ns, primary_model_ns):

    for mnxr in primary_model_ns.mnxr_to_add_list:

        kegg_id = utils.get_keggid_from_mnxr(mnxr, io_ns, primary_model_ns)

        logging.debug("--------------------")
        logging.debug("Reaction to be added: %s; %s", mnxr, kegg_id)

        rxn = io_ns.mnxref.reactions.get_by_id(mnxr)
        modelPrunedGPR.add_reaction(rxn)

        #Re-define ID
        # Assignment of the same reaction ID causes ValueError.
        #This happens if the model has KEGG reaction IDs (e.g., 'R00227' in 'nsal').
        try:
            rxn.id = kegg_id
        except ValueError as e:
            logging.error(e)

        #Re-define Name
        rxn.name = primary_model_ns.rxnid_info_dict[kegg_id]['NAME']

        #'add_reaction' requires writing/reloading of the model
        modelPrunedGPR = utils.stabilize_model(
                modelPrunedGPR, io_ns.outputfolder5, kegg_id)

        rxn = modelPrunedGPR.reactions.get_by_id(kegg_id)

        # NOTE:
        # Writing the reaction notes (e.g., 'GENE ASSOCIATION' and 'SUBSYSTEM')
        #seems to be unstable.
        # The modified reaction notes are missing in resulting SBML file.
        # This bug appears to happen if the new reaction comes from another SBML file.
        # Addition of a newly built reaction using Metabolite object does not seem to have
        #this potential bug.

        #GPR association
        if len(primary_model_ns.rxnid_locusTag_dict[kegg_id]) == 1:
            gpr = '( %s )' %(primary_model_ns.rxnid_locusTag_dict[kegg_id][0])
        else:
            count = 1
            for locusTag in primary_model_ns.rxnid_locusTag_dict[kegg_id]:

                #Check whether the submitted gbk file contains "/product" for CDS
                if locusTag in io_ns.targetGenome_locusTag_prod_dict:

                    #Consider "and" relationship in the GPR association
                    if 'subunit' in io_ns.targetGenome_locusTag_prod_dict[locusTag]:
                        count += 1
            if count == len(primary_model_ns.rxnid_locusTag_dict[kegg_id]):
                gpr = ' and '.join(primary_model_ns.rxnid_locusTag_dict[kegg_id])
            else:
                gpr = ' or '.join(primary_model_ns.rxnid_locusTag_dict[kegg_id])
            gpr = '( %s )' %(gpr)

        rxn.gene_reaction_rule = gpr

        # Subsystem
        rxn.subsystem = primary_model_ns.rxnid_info_dict[kegg_id]['PATHWAY']

        modelPrunedGPR = utils.stabilize_model(
                modelPrunedGPR, io_ns.outputfolder5, kegg_id)

        logging.debug("Number of reactions in the model: %s",
                len(modelPrunedGPR.reactions))

        #Check model prediction consistency
        target_exrxnid_flux_dict = utils.get_exrxnid_flux(modelPrunedGPR,
                io_ns.template_exrxnid_flux_dict)
        exrxn_flux_change_list = utils.check_exrxn_flux_direction(
                io_ns.template_exrxnid_flux_dict,
                target_exrxnid_flux_dict, config_ns)

        if 'F' in exrxn_flux_change_list:
            #'remove_reactions' does not seem to require
            #writing/reloading of the model
            modelPrunedGPR.remove_reactions(rxn)

    target_model = copy.deepcopy(modelPrunedGPR)
    return target_model


def get_rxn_newComp_list_from_model(model, io_ns):

    rxn_newComp_list = []

    for locustag in io_ns.locustag_comp_dict:
        for j in range(len(model.reactions)):
            rxn = model.reactions[j]

            if locustag in str(rxn.genes) and rxn.id not in rxn_newComp_list:
                rxn_newComp_list.append(rxn.id)

    return rxn_newComp_list


def create_rxn_newComp(rxn_newComp_list, model, io_ns, primary_model_ns):

    primary_model_ns.rxn_newComp_fate_dict = {}
    added_rxn_newComp_list = []

    for rxnid in rxn_newComp_list:
        rxn = model.reactions.get_by_id(rxnid)

        for locustag in io_ns.locustag_comp_dict:
            if locustag in str(rxn.genes):
                for i in range(len(io_ns.locustag_comp_dict[locustag])):
                    rxn_newCompt_dict = {}
                    rxn_comp_list = []

                    for metab in rxn.metabolites:
                        if metab.compartment not in rxn_comp_list:
                            rxn_comp_list.append(metab.compartment)

                    if len(rxn_comp_list) > 1:

                        logging.debug(
                            "Reaction %s for %s has metabolites with multiple compartments ('%s'): manual addition needed", rxnid, locustag, rxn_comp_list)

                        primary_model_ns.rxn_newComp_fate_dict[rxnid] = \
                                "Locus tag: %s; Exisiting compartment: %s; New comparment: %s; Not added - metabolites with multiple compartments" %(locustag, rxn_comp_list, io_ns.locustag_comp_dict[locustag])

                    elif len(rxn_comp_list) == 1 and \
                            rxn_comp_list[0] == io_ns.locustag_comp_dict[locustag][i]:

                        logging.debug(
                        "Reaction %s for %s with the compartment ('%s'): already exists",
                        rxnid, locustag, io_ns.locustag_comp_dict[locustag][i])

                        primary_model_ns.rxn_newComp_fate_dict[rxnid] = \
                                "Locus tag: %s; Exisiting compartment: %s; New comparment: %s; Not added - already exists" %(locustag, rxn_comp_list, io_ns.locustag_comp_dict[locustag])

                    elif len(rxn_comp_list) == 1 and \
                            rxn_comp_list[0] != io_ns.locustag_comp_dict[locustag][i]:

                        new_rxn_id = \
                                ''.join([rxn.id, io_ns.locustag_comp_dict[locustag][i]])

                        # rxn.metabolites extracts metabolites and their coeff's
                        #from the corresponding reaction
                        for metab in rxn.metabolites:
                            if metab.id in model.metabolites:
                                new_metab_id = '_'.join(
                                            [metab.id[:-2],
                                            io_ns.locustag_comp_dict[locustag][i]])
                                new_metab = Metabolite(
                                        new_metab_id,
                                        name = metab.name,
                                        formula = metab.formula,
                                        compartment = io_ns.locustag_comp_dict[locustag][i]
                                        )
                                rxn_newCompt_dict[new_metab] = \
                                        float(rxn.metabolites[metab])

                        new_rxn_id = \
                                ''.join([rxn.id, io_ns.locustag_comp_dict[locustag][i]])
                        rxn_newComp = Reaction(new_rxn_id)
                        rxn_newComp.name = rxn.name
                        rxn_newComp.subsystem = rxn.subsystem
                        rxn_newComp.lower_bound = rxn.lower_bound
                        rxn_newComp.upper_bound = rxn.upper_bound
                        rxn_newComp.reversibility = rxn.reversibility
                        rxn_newComp.gene_reaction_rule = rxn.gene_reaction_rule
                        rxn_newComp.add_metabolites(rxn_newCompt_dict)

                        res = utils.check_duplicate_rxn(model, rxn_newComp)

                        if res == 'unique' and rxn_newComp.id not in model.reactions:
                            #'add_reaction' requires writing/reloading of the model
                            model.add_reactions(rxn_newComp)
                            model = utils.stabilize_model(
                                    model, io_ns.outputfolder5, rxn_newComp.id)
                            added_rxn_newComp_list.append(rxn_newComp.id)

                            logging.debug(
                            "Reaction %s for %s with the compartment ('%s'): added ",
                            rxnid, locustag, io_ns.locustag_comp_dict[locustag][i])

                            primary_model_ns.rxn_newComp_fate_dict[rxnid] = \
                                "Locus tag: %s; Exisiting compartment: %s; New comparment: %s; Added" %(locustag, rxn_comp_list, io_ns.locustag_comp_dict[locustag])

                        elif res == 'duplicate' or rxn_newComp.id in model.reactions:
                            logging.debug(
                                    "Reaction %s for %s with the compartment ('%s'): already exists", rxnid, locustag, io_ns.locustag_comp_dict[locustag][i])

                            primary_model_ns.rxn_newComp_fate_dict[rxnid] = \
                                    "Locus tag: %s; Exisiting compartment: %s; New comparment: %s; Not added - already exists" %(locustag, rxn_comp_list, io_ns.locustag_comp_dict[locustag])

    return model, added_rxn_newComp_list


def remove_inactive_rxn_newComp(added_rxn_newComp_list, model, io_ns, primary_model_ns):
    rxn_newComp_list2 = []
    primary_model_ns.inactive_rxn_newComp_list = []

    for rxn in added_rxn_newComp_list:
        rxn_newComp_list2.append(model.reactions.get_by_id(rxn))

    primary_model_ns.inactive_rxn_newComp_list = \
            cobra.flux_analysis.variability.find_blocked_reactions(
                                                        model,
                                                        reaction_list=rxn_newComp_list2)

    model.remove_reactions(primary_model_ns.inactive_rxn_newComp_list)
    model = utils.stabilize_model(model, io_ns.outputfolder5, '')

    logging.debug("Following reactions with new compartments have been removed from the model as they carry no fluxes: %s", primary_model_ns.inactive_rxn_newComp_list)

    return model

