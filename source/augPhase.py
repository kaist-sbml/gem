'''
2014-2015
Hyun Uk Kim, Tilmann Weber, Jae Yong Ryu and Kyu-Sang Hwang
'''

import copy
import logging
import os
import pickle
import re
import urllib2
from Bio import SeqIO
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
from cobra.manipulation.delete import prune_unused_metabolites

#Converts ncbi_geneid(e.g. 944762) to species_geneid (e.g. b0031) in KEGG database by using KEGGAPI
#Input: ncbi_gene_id (string) (e.g., 1217641)
#Return: Locus_tag ID in KEGG database (e.g., sma:SAV_7535)
def get_species_locusTag(ncbi_geneid):
    url = "http://rest.kegg.jp/conv/genes/ncbi-gi:%s"%(ncbi_geneid)
    #Open and read data for the results of query in url address
    data = urllib2.urlopen(url).read()
    sptlist = data.strip().split()
    logging.debug("sptlist: %s" %sptlist)
    species_locusTag = sptlist[1]
    #print "species_locusTag:", species_locusTag
    return species_locusTag


#Retrieves EC numbers relevant to their  species_locusTag from KEGG
#Input: Species_locusTag in KEGG database (e.g., sma:SAV_7535)
#Output: A set of EC numbers in list form (e.g., ['4.1.3.6 ', '4.1.3.34]'])
def get_ECNumberList_from_locusTag(species_locusTag):
    url = "http://rest.kegg.jp/get/%s"%(species_locusTag)
    gene_info_text = urllib2.urlopen(url).read()
    split_text = gene_info_text.strip().split('\n')
    for line in split_text:
        sptlist = line.split()
        if sptlist[0] == 'ORTHOLOGY':
            ECNumbers = re.findall(r'[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+[\s\]]', line)
            return ECNumbers
    return []


#Retrieves a list of reaction IDs using their EC numbers from KEGG
#Input: E.C number in string form (e.g., 4.1.3.6)
#Output: reactionID in list form (e.g., ['R00362'])
def get_rxnid_from_ECNumber(enzymeEC):
    url = "http://rest.kegg.jp/get/enzyme:%s"%(enzymeEC)
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
        #Otherwise reactions have unspecified molecules having R groups
        if sptlist[0].strip() == 'PATHWAY': 
            PATHWAY = ' '.join(sptlist[1:])
	    return {'NAME':NAME, 'DEFINITION':DEFINITION, 'EQUATION':EQUATION, 'ENZYME':ENZYME, 'PATHWAY':PATHWAY}


def get_targetGenome_locusTag_ec_nonBBH_dict(options):
    targetGenome_locusTag_ec_nonBBH_dict = {}

    for locusTag in options.nonBBH_list:
	if locusTag in options.targetGenome_locusTag_ec_dict.keys():
	    targetGenome_locusTag_ec_nonBBH_dict[locusTag] = options.targetGenome_locusTag_ec_dict[locusTag] 
    options.targetGenome_locusTag_ec_nonBBH_dict = targetGenome_locusTag_ec_nonBBH_dict


#Two nested function calling four functions above
def make_all_rxnInfo_fromRefSeq(options):
    rxnid_info_dict = {}
    rxnid_locusTag_dict = {}
 
    for locusTag in options.targetGenome_locusTag_ec_nonBBH_dict.keys():
	for enzymeEC in options.targetGenome_locusTag_ec_nonBBH_dict[locusTag]:
	    logging.debug("EC_number for locusTag: %s, %s" %(locusTag, enzymeEC))

            #KEGG REST does not accept unspecific EC_number: e.g., 3.2.2.-
	    if '-' not in enzymeEC:
                rxnid_list = get_rxnid_from_ECNumber(enzymeEC)
                for rxnid in rxnid_list:
                    rxnid_info_dict[rxnid] = get_rxnInfo_from_rxnid(rxnid)

		    if rxnid not in rxnid_locusTag_dict.keys():
		        rxnid_locusTag_dict[rxnid] = [(locusTag)]

                    #Appends additional different genes to the same reaction ID
		    elif rxnid in rxnid_locusTag_dict.keys():
		        rxnid_locusTag_dict[rxnid].append((locusTag))
                    #print locusTag, rxnid, rxnid_info_dict[rxnid], "\n"

    options.rxnid_info_dict = rxnid_info_dict
    options.rxnid_locusTag_dict = rxnid_locusTag_dict


#Output: a list of MNXRs available in the modelPrunedGPR
def get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, bigg_mnxr_dict):
    modelPrunedGPR_mnxr_list = []

    index_last = len(modelPrunedGPR.reactions)
    index_last = index_last - 1
    index = 0

    while index <= index_last:
        rxn = modelPrunedGPR.reactions[index].id
        if rxn in bigg_mnxr_dict.keys():
            modelPrunedGPR_mnxr_list.append(bigg_mnxr_dict[rxn])
        index+=1

    options.modelPrunedGPR_mnxr_list = modelPrunedGPR_mnxr_list


def check_existing_rxns(kegg_mnxr_dict, modelPrunedGPR_mnxr_list, options):
    rxnid_to_add_list =[]

    for rxnid in options.rxnid_info_dict.keys():
        #Consider only reactions mapped in pathways
	if rxnid in kegg_mnxr_dict.keys():
            kegg_mnxr = kegg_mnxr_dict[rxnid]

            #Check with reactions in the template model through MNXref
            if kegg_mnxr not in modelPrunedGPR_mnxr_list and rxnid not in rxnid_to_add_list:
                rxnid_to_add_list.append(rxnid)

    rxnid_to_add_list = list(sorted(set(rxnid_to_add_list)))
    options.rxnid_to_add_list = rxnid_to_add_list


#Output: MNXR for the reactions to add, converted from KEGG rxnid
def get_mnxr_using_kegg(options, kegg_mnxr_dict):
    mnxr_to_add_list = []
    for rxnid in options.rxnid_to_add_list:
	if rxnid in kegg_mnxr_dict.keys():
	    mnxr_to_add_list.append(kegg_mnxr_dict[rxnid])

    options.mnxr_to_add_list = mnxr_to_add_list


def get_correct_metab_coeff(converted_metab_id, metab_coeff, metab_type, mnxm_coeff_dict, mnxm_metab_list):

    #If the same metabolite appears multiple times as either substrates or products,
    #their stoichiometric coeff's are all added up
    if converted_metab_id in mnxm_metab_list:
        overlap_metab_coeff = float(mnxm_coeff_dict[converted_metab_id])
        mnxm_coeff_dict[converted_metab_id] = overlap_metab_coeff+float(metab_coeff)*-1
    else:
        if metab_type == 'substrate':
            mnxm_coeff_dict[converted_metab_id] = float(metab_coeff)*-1
        elif metab_type == 'product':
            mnxm_coeff_dict[converted_metab_id] = float(metab_coeff)

    return mnxm_coeff_dict


#Check if the same metabolite appears as a substrate and a product
def check_overlap_subs_prod(mnxm_subs_list, mnxm_prod_list):

    for each_substrate in mnxm_subs_list:
        if each_substrate in mnxm_prod_list:
            overlap_check = True
            break
        else:
            overlap_check = False

    return overlap_check


#Creating: e.g., {'R03232': {'f1p': -1.0, 'C04261': 1.0, 'fru': 1.0, 'C00615': -1.0}}
#Metabolites are presented primarily with bigg, otherwise with KEGG
def extract_rxn_mnxm_coeff(options, mnxr_rxn_dict, mnxm_bigg_compound_dict, mnxm_kegg_compound_dict, mnxr_kegg_dict):
    rxnid_mnxm_coeff_dict = {}
    mnxm_coeff_dict = {}

    for mnxr in options.mnxr_to_add_list:
	unparsed_equation = mnxr_rxn_dict[mnxr]
	logging.debug(unparsed_equation)

        #"substrates" and "products" contain stoichiometric coeff of each compound
	sptReaction = unparsed_equation.split('=')
	substrates = sptReaction[0].strip()
	products = sptReaction[1].strip()

        #Discards polymerization reactions with undefinite coeff's
        #e.g., 1 MNXM9 + (n+2) MNXM90033 = 1 MNXM5617 + (n) MNXM90033
	if '(' not in substrates and '(' not in products:
            #Creating:
            #e.g., {bigg compoundID:(-1)coeff}, {kegg compoundID:(-1)coeff}
            #or {mnxm:(-1)coeff}
	    substrates = substrates.split(' + ')
	    mnxm_coeff_dict = {}
            mnxm_subs_list = []
            mnxm_prod_list = []

	    for substrate in substrates:
                metab_type = 'substrate'
	        substrate = substrate.split()

	        if substrate[1] in mnxm_bigg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_bigg_compound_dict[substrate[1]], substrate[0], metab_type, mnxm_coeff_dict, mnxm_subs_list)
                    mnxm_subs_list.append(mnxm_bigg_compound_dict[substrate[1]])

	        elif substrate[1] in mnxm_kegg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_kegg_compound_dict[substrate[1]], substrate[0], metab_type, mnxm_coeff_dict, mnxm_subs_list)
                    mnxm_subs_list.append(mnxm_kegg_compound_dict[substrate[1]])

	        else:
                    mnxm_coeff_dict = get_correct_metab_coeff(substrate[1], substrate[0], metab_type, mnxm_coeff_dict, mnxm_subs_list)
                    mnxm_subs_list.append(substrate[1])

            #Creating:
            #e.g., {bigg compoundID:coeff}, {kegg compoundID:coeff} or {mnxm:coeff}
	    products = products.split(' + ')
	    for product in products:
                metab_type = 'product'
	        product = product.split()

	        if product[1] in mnxm_bigg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_bigg_compound_dict[product[1]], product[0], metab_type, mnxm_coeff_dict, mnxm_prod_list)
                    mnxm_prod_list.append(mnxm_bigg_compound_dict[product[1]])

	        elif product[1] in mnxm_kegg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_kegg_compound_dict[product[1]], product[0], metab_type, mnxm_coeff_dict, mnxm_prod_list)
                    mnxm_prod_list.append(mnxm_kegg_compound_dict[product[1]])

	        else:
                    mnxm_coeff_dict = get_correct_metab_coeff(product[1], product[0], metab_type, mnxm_coeff_dict, mnxm_prod_list)
                    mnxm_prod_list.append(product[1])

            #Check overlapping metabolites as a substrate and a product
            #e.g., ATP + ADP <=> ADP + ATP 
            overlap_check = check_overlap_subs_prod(mnxm_subs_list, mnxm_prod_list)
            if overlap_check == True:
                continue
            else:
                #Creating: 
                #e.g., {'R03232': {'f1p': -1.0, 'C04261': 1.0, 'fru': 1.0, 'C00615': -1.0}}
	        rxnid_mnxm_coeff_dict[mnxr_kegg_dict[mnxr]] = mnxm_coeff_dict

    options.rxnid_mnxm_coeff_dict = rxnid_mnxm_coeff_dict


def add_nonBBH_rxn(modelPrunedGPR, bigg_mnxm_compound_dict, kegg_mnxm_compound_dict, mnxm_compoundInfo_dict, template_exrxnid_flux_dict, options):

    for rxnid in options.rxnid_mnxm_coeff_dict.keys():
	logging.debug(rxnid)

        if options.rxnid_info_dict[rxnid] != None:

            #ID
	    rxn = Reaction(rxnid)

            #Name
            #Some reaction IDs do not have NAME despite the presence of PATHWAY
	    rxn.name = options.rxnid_info_dict[rxnid]['NAME']
  
            #Reversibility / Lower and upper bounds
	    rxn.lower_bound = -1000
	    rxn.uppwer_bound = 1000

            #Metabolites and their stoichiometric coeff's
	    for metab in options.rxnid_mnxm_coeff_dict[rxnid]:
	        metab_compt = '_'.join([metab,'c'])

                #Adding metabolites already in the model
	        if metab_compt in modelPrunedGPR.metabolites:
		    rxn.add_metabolites({modelPrunedGPR.metabolites.get_by_id(metab_compt):options.rxnid_mnxm_coeff_dict[rxnid][metab]})

                #Adding metabolites with bigg compoundID, but not in the model
	        elif metab in bigg_mnxm_compound_dict.keys():
		    mnxm = bigg_mnxm_compound_dict[metab]
		    metab_compt = Metabolite(metab, formula = mnxm_compoundInfo_dict[mnxm][1], name = mnxm_compoundInfo_dict[mnxm][0], compartment='c')
		    rxn.add_metabolites({metab_compt:options.rxnid_mnxm_coeff_dict[rxnid][metab]})

                #Adding metabolites with MNXM and not in the model
	        else:
		    mnxm = kegg_mnxm_compound_dict[metab]
		    metab_compt = Metabolite(mnxm, formula = mnxm_compoundInfo_dict[mnxm][1], name = mnxm_compoundInfo_dict[mnxm][0], compartment='c')
		    rxn.add_metabolites({metab_compt:options.rxnid_mnxm_coeff_dict[rxnid][metab]})

            #GPR association
	    if len(options.rxnid_locusTag_dict[rxnid]) == 1:
	        gpr = '( %s )' %(options.rxnid_locusTag_dict[rxnid][0])
	    else:
	        count = 1
	        for locusTag in options.rxnid_locusTag_dict[rxnid]:

                    #Check whether the submitted gbk file contains "/product" for CDS
                    if locusTag in options.targetGenome_locusTag_prod_dict:

                        #Considers "and" relationship in the GPR association
		        if 'subunit' in options.targetGenome_locusTag_prod_dict[locusTag]:
		            count += 1
	        if count == len(options.rxnid_locusTag_dict[rxnid]):
	            gpr = ' and '.join(options.rxnid_locusTag_dict[rxnid])
 	        else:
	            gpr = ' or '.join(options.rxnid_locusTag_dict[rxnid])
	        gpr = '( %s )' %(gpr)
	    rxn.add_gene_reaction_rule(gpr)

            #Subsystem
	    rxn.subsystem = options.rxnid_info_dict[rxnid]['PATHWAY']

            #E.C. number: not available feature in COBRApy
            #Objective coeff: default
	    rxn.objective_coefficient = 0

            #Add a reaction to the model if it does not affect Exchange reaction flux direction
	    modelPrunedGPR.add_reaction(rxn)

            write_cobra_model_to_sbml_file(modelPrunedGPR, "./%s/3_temp_models/modelPrunedGPR.xml" %options.output)
            modelPrunedGPR = create_cobra_model_from_sbml_file("./%s/3_temp_models/modelPrunedGPR.xml" %options.output)

            target_exrxnid_flux_dict = get_exrxnid_flux(modelPrunedGPR, template_exrxnid_flux_dict)
            exrxn_flux_change_list = check_exrxn_flux_direction(template_exrxnid_flux_dict, target_exrxnid_flux_dict)

            if 'F' in exrxn_flux_change_list:
	        modelPrunedGPR.remove_reactions(rxn)
                
                write_cobra_model_to_sbml_file(modelPrunedGPR, "./%s/3_temp_models/modelPrunedGPR.xml" %options.output)
                modelPrunedGPR = create_cobra_model_from_sbml_file("./%s/3_temp_models/modelPrunedGPR.xml" %options.output)

    prune_unused_metabolites(modelPrunedGPR)
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
def check_exrxn_flux_direction(template_exrxnid_flux_dict, target_exrxnid_flux_dict):

    exrxn_flux_change_list = []

    for exrxn_id in template_exrxnid_flux_dict.keys():
        if exrxn_id in target_exrxnid_flux_dict.keys():
            template_exrxn_flux = template_exrxnid_flux_dict[exrxn_id]
            target_exrxn_flux = target_exrxnid_flux_dict[exrxn_id]
            ratio_exrxn_flux = float(target_exrxn_flux)/float(template_exrxn_flux)

            #Similar species are allowed to uptake nutrients within a decent range
            if float(target_exrxn_flux)*float(template_exrxn_flux) > 0.0 and 0.2 < ratio_exrxn_flux and ratio_exrxn_flux < 2.0:
                exrxn_flux_change_list.append('T')

            #Causing drastic changes in Exchange reaction fluxes (direction and/or magnitude)
            else:
                exrxn_flux_change_list.append('F')

    return exrxn_flux_change_list
