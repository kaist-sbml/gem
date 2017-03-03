
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

#This file generates metabolic reactions for the genes newly annotated to be present in the secondary metabolite-biosynthetic gene cluster from antiSMASH.

import logging
import os
import pickle
from Bio import SeqIO
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from antismash_monomer_info import get_std_id_from_antismash_id


def get_cluster_location(cluster_nr, options):

    for feature in options.seq_record.features:

        if feature.type == 'cluster':
            cluster_number = 'Cluster number: %s' %cluster_nr
            options.cluster_number = cluster_number

            if options.cluster_number in feature.qualifiers.get('note'):
                options.cluster_loc1 = feature.location.start
                options.cluster_loc2 = feature.location.end

#Exract all the information associated with a particular locus_tag for the selected cluster
def get_cluster_info_from_seq_record(options):

    cluster_info_dict = {}

    for feature in options.seq_record.features:

        if feature.type == 'CDS':
            if feature.location.start >= options.cluster_loc1 \
                    and feature.location.end <= options.cluster_loc2:
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                if feature.qualifiers.get('sec_met'):
                    cluster_info_dict[qualifier_locus_tag] = \
                            feature.qualifiers.get('sec_met')

    #logging.debug('cluster_info_dict: %s' %cluster_info_dict)
    options.cluster_info_dict = cluster_info_dict


def get_cluster_product(cluster_nr, options):

    for feature in options.seq_record.features:

        #Retrieving "Cluster number"
        if feature.type == 'cluster':
            qualifier_cluster = feature.qualifiers.get('note')
            if options.cluster_number in qualifier_cluster:
                product = feature.qualifiers.get('product')[0]

    #Handle legacy problem
    product2 = product.replace("-","_")

    if float(cluster_nr) < 10:
        product3 = "Cluster0"+str(cluster_nr)+"_"+product2
    else:
        product3 = "Cluster"+str(cluster_nr)+"_"+product2

    #print product, "\n"
    options.product = product3


#TODO: Check how this information can be used in future
# Addition of currency metabolites for each domain
#rather than previous way of domain combinations

# domain information(nrps) :
# Condensation    Condensation domain
# Condensation_DCL    Condensation domain that links L-amino acid to peptide ending with D    -amino acid
# Condensation_LCL    Condensation domain that links L-amino acid to peptide ending with L    -amino acid
# Condensation_Dual     Dual condensation / epimerization domain
# Condensation_Starter    Starter condensation domain
# CXglyc    Putatively inactive glycopeptide condensation-like domain
# Cglyc    Glycopeptide condensation domain
# Heterocyclization    Heterocyclization domain
# Epimerization    Epimerization domain
# AMP-binding    Adenylation domain
# A-OX     Adenylation domain with integrated oxidase
# PCP     Peptidyl-carrier protein domain
# ACP    4'-phosphopantetheinyl transferase
# NRPS-COM_Nterm    NRPS COM domain Nterminal
# NRPS-COM_Cterm    NRPS COM domain Cterminal

# domain information(pks) :
# AT    acyltransferase
# KS    ketosynthase
# ACP    acyl carrier protein
# KR    ketoreductase
# DH    dehydratase
# ER    enolase
# cMT    methyltransferase
# TD    thiolesterase domain
def get_cluster_domain(options):

    locustag_domain_dict = {}
    locustag_kr_dict = {}

    for each_gene in options.cluster_info_dict.keys():

        domain_count = 0
        sec_met_domain_list = []
        kr_domain_info_dict = {}

        for sec_met_info in options.cluster_info_dict[each_gene]:

            if "NRPS/PKS Domain" in sec_met_info:
                #e.g., 'NRPS/PKS Domain: Condensation_LCL (36-335)'
                domain_info1 = sec_met_info.split('. ')[0]

                #e.g., 'Condensation_LCL (36-335)'
                domain_info2 = domain_info1.split(':')[1]

                #e.g., 'Condensation_LCL'
                each_sec_met_domain = domain_info2.split()[0]

                # Count the number of domains
                if float(domain_count) < 10:
                    domain_number = "_D0" + str(domain_count)
                else:
                    domain_number = "_D" + str(domain_count)

                domain_number = each_sec_met_domain + domain_number
                sec_met_domain_list.append(domain_number)
                domain_count += 1

                #Collects KR activity information
                if each_sec_met_domain == "PKS_KR":
                    sec_met_info_list = sec_met_info.split('; ')
                    for kr_chars in sec_met_info_list:
                        if 'Predicted KR activity' in kr_chars:
                            kr_info_list = kr_chars.split(': ')
                    dm_kr_activity = kr_info_list[1]
                    kr_domain_info_dict[domain_number] = dm_kr_activity

        locustag_domain_dict[each_gene] = sec_met_domain_list
        locustag_kr_dict[each_gene] = kr_domain_info_dict

    #logging.debug('locustag_domain_dict: %s' %locustag_domain_dict)
    options.locustag_domain_dict = locustag_domain_dict
    options.locustag_kr_dict = locustag_kr_dict


#Output: e.g., {'SAV_943_M1':['mmal', 'Ethyl_mal', 'pk']}
def get_cluster_monomers(options):

    locustag_monomer_dict = {}
    for each_gene in options.cluster_info_dict.keys():
        module_count = 0

        for sec_met_info in options.cluster_info_dict[each_gene]:

            if 'Substrate specificity predictions' in sec_met_info:
                #type(sec_met_info) is string
                #Convert 'sec_met_info' into a list
                sec_met_info_list = sec_met_info.split(';')
                for each_sec_met_info in sec_met_info_list:
                    if 'Substrate specificity predictions' in each_sec_met_info:
                        pred_monomer_list = []
                        #Following statement produces a list with 2 elements:
                        #e.g., [' Substrate specificity predictions',
                        #' gly (NRPSPredictor2 SVM), gly (Stachelhaus code),
                        #gly (Minowa), gly (consensus)']
                        monomer_list = each_sec_met_info.split(':')
                        for monomer in monomer_list:
                            #Make sure to include a space ''
                            #for monomers predicted from different engines:
                            #e.g., 'orn,lys,arg (NRPSPredictor2 SVM),
                            #lys (Stachelhaus code), leu (Minowa), nrp (consensus)'
                            if 'Substrate specificity predictions' not in monomer \
                                    and ', ' in monomer:
                                monomer_list2 = monomer.split(', ')
                                for monomer2 in monomer_list2:
                                    monomer2_list = monomer2.split('(')
                                    monomer3 = monomer2_list[0].strip()
                                    pred_monomer_list.append(monomer3)

                module_number = each_gene + '_M' + str(module_count)
                locustag_monomer_dict[module_number] = pred_monomer_list
                module_count += 1

    options.locustag_monomer_dict = locustag_monomer_dict


#Add stoichiometric coeff's of monomers
#Output: e.g., {'mmalcoa': -4, 'malcoa': -7}
def get_all_metab_coeff(options):

    metab_coeff_dict = {}

    for each_module in options.locustag_monomer_dict.keys():
        #locustag_monomer_dict[each_module] for nrps
        #Position [0]: NRPSPredictor2 SVM
        #Position [1]: Stachelhaus code
        #Position [2]: Minowa
        #Position [3]: consensus
        if len(options.locustag_monomer_dict[each_module]) == 4:

            #In case "consensus" is not reached:
            if options.locustag_monomer_dict[each_module][3] == 'nrp':
                #From NRPSPredictor2 SVM
                #Not considered: e.g., NRPSPredictor2 SVM: val,leu,ile,abu,iva
                #Checked by ',' in aSid_met2
                aSid_met2 = options.locustag_monomer_dict[each_module][0]
                if aSid_met2 != 'hydrophobic-aliphatic' \
                        and aSid_met2 != 'hydrophilic' \
                        and aSid_met2 != 'hydrophobic-aromatic' \
                        and aSid_met2 != 'N/A' \
                        and ',' not in aSid_met2:
                    biggid_met2 = get_std_id_from_antismash_id(aSid_met2)

                    #In case of non-consensus, NRPSPredictor2 SVM is considered
                    if biggid_met2 not in metab_coeff_dict:
                        metab_coeff_dict[biggid_met2] = -1
                    else:
                        metab_coeff_dict[biggid_met2] -= 1

                elif aSid_met2 == 'hydrophobic-aliphatic' \
                        or aSid_met2 == 'hydrophilic' \
                        or aSid_met2 == 'hydrophobic-aromatic' \
                        or aSid_met2 == 'N/A' \
                        or ',' in aSid_met2:
                    #If NRPSPredictor2 SVM has invalid monomer, then Minowa is considered
                    aSid_met4 = options.locustag_monomer_dict[each_module][2]
                    if aSid_met4 != 'hydrophobic-aliphatic' \
                            and aSid_met4 != 'hydrophilic' \
                            and aSid_met4 != 'hydrophobic-aromatic' \
                            and aSid_met4 != 'N/A':
                        biggid_met4 = get_std_id_from_antismash_id(aSid_met4)

                        if biggid_met4 not in metab_coeff_dict:
                            metab_coeff_dict[biggid_met4] = -1
                        else:
                            metab_coeff_dict[biggid_met4] -= 1

                    #If Minowa has invalid monomer, then Stachelhaus code is considered
                    elif aSid_met4 == 'hydrophobic-aliphatic' \
                            or aSid_met4 == 'hydrophilic' \
                            or aSid_met4 == 'hydrophobic-aromatic' \
                            or aSid_met4 == 'N/A':
                        aSid_met3 = options.locustag_monomer_dict[each_module][1]
                        biggid_met3 = get_std_id_from_antismash_id(aSid_met3)

                        if biggid_met3 not in metab_coeff_dict:
                            metab_coeff_dict[biggid_met3] = -1
                        else:
                            metab_coeff_dict[biggid_met3] -= 1

            #In case "consensus" is reached:
            elif options.locustag_monomer_dict[each_module][3] != 'nrp':
                aSid_met5 = options.locustag_monomer_dict[each_module][3]
                biggid_met5 = get_std_id_from_antismash_id(aSid_met5)
                #print "aSid_met5", aSid_met5, biggid_met5

                if biggid_met5 not in metab_coeff_dict:
                    metab_coeff_dict[biggid_met5] = -1
                else:
                    metab_coeff_dict[biggid_met5] -= 1

        #locustag_monomer_dict[each_module] for pks
        #Position [0]: PKS signature
        #Position [1]: Minowa
        #Position [2]: consensus
        elif len(options.locustag_monomer_dict[each_module]) == 3:

            if len(options.locustag_monomer_dict[each_module]) < 3:
                continue

            #In case "consensus" is not reached:
            if options.locustag_monomer_dict[each_module][2] == 'pk':

                #From PKS signature
                #In case of non-consensus, PKS signature is considered
                aSid_met6 = options.locustag_monomer_dict[each_module][0]
                if aSid_met6 != 'N/A' and aSid_met6 != 'mal_or_prop':
                    biggid_met6 = get_std_id_from_antismash_id(aSid_met6)
                    #print "aSid_met6", aSid_met6, biggid_met6

                    if biggid_met6 not in metab_coeff_dict:
                        metab_coeff_dict[biggid_met6] = -1
                    else:
                        metab_coeff_dict[biggid_met6] -= 1

                #If PKS signature has invalid monomer, then Minowa is considered
                else:
                    aSid_met7 = options.locustag_monomer_dict[each_module][1]
                    if aSid_met7 != 'inactive':
                        biggid_met7 = get_std_id_from_antismash_id(aSid_met7)

                        if biggid_met7 not in metab_coeff_dict:
                            metab_coeff_dict[biggid_met7] = -1
                        else:
                            metab_coeff_dict[biggid_met7] -= 1

            #In case "consensus" is reached:
            elif options.locustag_monomer_dict[each_module][2] != 'pk':
                aSid_met8 = options.locustag_monomer_dict[each_module][2]

                #Original monomers are considered for those reduced by KR, DH and/or ER
                if aSid_met8 == 'ohmal' or aSid_met8 == 'ccmal' or aSid_met8 == 'redmal':
                    aSid_met8 = 'mal'
                elif aSid_met8 == 'ohmmal' or aSid_met8 == 'ccmmal' \
                        or aSid_met8 == 'redmmal':
                    aSid_met8 = 'mmal'
                elif aSid_met8 == 'ohmxmal' or aSid_met8 == 'ccmxmal' \
                        or aSid_met8 == 'redmxmal':
                    aSid_met8 = 'mxmal'
                elif aSid_met8 == 'ohemal' or aSid_met8 == 'ccemal' \
                        or aSid_met8 == 'redemal':
                    aSid_met8 = 'emal'

                biggid_met8 = get_std_id_from_antismash_id(aSid_met8)
                #print "aSid_met8", aSid_met8, biggid_met8

                if biggid_met8 not in metab_coeff_dict:
                    metab_coeff_dict[biggid_met8] = -1
                else:
                    metab_coeff_dict[biggid_met8] -= 1

    #Add secondary metabolite product to the reaction
    metab_coeff_dict[options.product] = 1

    logging.debug('metab_coeff_dict: %s' %metab_coeff_dict)
    options.metab_coeff_dict = metab_coeff_dict


def get_pickles(options):

    if not hasattr(options, 'mnxref'):
        mnxref = pickle.load(open('./gems/io/data/input2/MNXref.p','rb'))
        options.mnxref = mnxref

    if not hasattr(options, 'mnxm_compoundInfo_dict'):
        mnxm_compoundInfo_dict = pickle.load(
                open('./gems/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
        options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict


def add_sec_met_rxn(target_model, options):

    #ID
    rxn = Reaction(options.product)

    #Reversibility / Lower and upper bounds
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Metabolites and their stoichiometric coeff's
    for metab in options.metab_coeff_dict.keys():

        #Consider only metabolites consumed or produced
        if options.metab_coeff_dict[metab] != 0:
            metab_compt = '_'.join([metab,'c'])

            #Adding metabolites already in the model
            if metab_compt in target_model.metabolites:
                logging.debug('Metabolite %s already present in the model', metab_compt)
                rxn.add_metabolites({target_model.metabolites.get_by_id(
                    metab_compt):options.metab_coeff_dict[metab]})

            #Adding metabolites with bigg compoundID, but not in the model
            elif metab_compt in options.mnxref.metabolites:
                logging.debug('Metabolite %s available in the MNXref', metab_compt)
                metab_compt = options.mnxref.metabolites.get_by_id(metab_compt)
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

            elif 'Cluster' in metab:
                logging.debug("Secondary metabolite ('Cluster') %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

            #Add metabolite MNXM having no bigg ID to the model
            else:
                logging.debug("Metabolite (MNXM ID) %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt,
                        formula = options.mnxm_compoundInfo_dict[metab][1],
                        name = options.mnxm_compoundInfo_dict[metab][0],
                        compartment='c')
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

    #GPR association
    gpr_count = 0
    for each_gene in options.cluster_info_dict.keys():
        if gpr_count == 0:
            gpr_list = each_gene
            gpr_count += 1
        else:
            gpr_list = ' and '.join([gpr_list, each_gene])

    gpr_list = '( %s )' %(gpr_list)

    rxn.gene_reaction_rule = gpr_list

    #Adding the new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("Secondary metabolite reaction: %s" %rxn)
    logging.debug(rxn.reaction)

    ##############################
    #Creating a transport reaction
    #Creating reaction ID
    rxn = Reaction("Transport_" + options.product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Adding a substrate metabolite
    #print cobra_model.metabolites.get_by_id(str(product_c))
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(options.product+'_c')):-1})

    #Adding product metabolite(s)
    product_e = Metabolite(options.product+"_e", name='', compartment='e')
    rxn.add_metabolites({product_e:1})

    #Adding the new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("Transport reaction: %s" %rxn)
    logging.debug(rxn.reaction)

    ##############################
    #Creating an exchange reaction
    #Creating reaction ID
    rxn = Reaction("Ex_"+options.product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible 0: irreversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Adding a substrate metabolite
    #print cobra_model.metabolites.get_by_id(str(product_c))
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(product_e)):-1})

    #Adding the new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("Exchange reaction: %s" %rxn)
    logging.debug(rxn.reaction)

    return target_model


def check_producibility_sec_met(target_model, options):
    for rxn in target_model.reactions:
        rxn.objective_coefficient = 0

    #target_model.reactions.get_by_id('Biomass_SCO').objective_coefficient = 0
    target_model.reactions.get_by_id("Ex_"+options.product).objective_coefficient = 1

    #Model reloading and overwrtting are necessary for model stability
    #Without these, model does not produce an accurate prediction
    write_cobra_model_to_sbml_file(target_model,
            options.outputfolder5 + os.sep + 'target_model_%s.xml'
            %options.product, use_fbc_package=False)
    target_model = create_cobra_model_from_sbml_file(
            options.outputfolder5 + os.sep + 'target_model_%s.xml'
            %options.product)

    target_model.optimize()
    logging.debug("Flux: %s" %target_model.solution.f)

    target_model.reactions.get_by_id('Biomass_SCO').objective_coefficient = 1
    target_model.reactions.get_by_id("Ex_"+options.product).objective_coefficient = 0

    #return target_model.solution.f, product
    return target_model


def get_monomers_nonprod_sec_met(options):

    nonprod_sec_met_metab_list = []

    for metab in options.metab_coeff_dict.keys():
        #Exclude currency metabolites
        if options.metab_coeff_dict[metab] <0 \
                and metab != 'atp' \
                and metab != 'amp' \
                and metab != 'ppi' \
                and metab != 'amet' \
                and metab != 'ahcys' \
                and metab != 'fmn' \
                and metab != 'fmnh2' \
                and metab != 'nadp' \
                and metab != 'nadph' \
                and metab != 'h' \
                and metab != 'h2o' \
                and metab != 'hco3' \
                and metab != 'coa':
            nonprod_sec_met_metab_list.append(metab)

    return nonprod_sec_met_metab_list


def get_monomers_prod_sec_met(options):

    prod_sec_met_metab_list = []

    for metab in options.metab_coeff_dict.keys():
        #Exclude currency metabolites
        if options.metab_coeff_dict[metab] <0 \
                and metab != 'atp' \
                and metab != 'amp' \
                and metab != 'ppi' \
                and metab != 'amet' \
                and metab != 'ahcys' \
                and metab != 'fmn' \
                and metab != 'fmnh2' \
                and metab != 'nadp' \
                and metab != 'nadph' \
                and metab != 'h' \
                and metab != 'h2o' \
                and metab != 'hco3' \
                and metab != 'coa':
            prod_sec_met_metab_list.append(metab)

    return prod_sec_met_metab_list

