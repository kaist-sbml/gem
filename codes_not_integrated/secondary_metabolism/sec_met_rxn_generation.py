'''
2015 Kyu-Sang Hwang
2014-2015 Hyun Uk Kim

This file generates metabolic reactions for the genes newly annotated to be present in the secondary metabolite-biosynthetic gene cluster from antiSMASH.
'''

from Bio import SeqIO
from sets import Set
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file,write_cobra_model_to_sbml_file
from MNX_checker2 import COBRA_TemplateModel_checking_MNXref_metabolites, fix_legacy_id
from general_sec_met_info import determine_module, get_biggid_from_aSid, get_metab_coeff_dict
import pickle
import copy


#Exracts all the information associated wiht a particular locus_tag
#def get_locustag_info_from_cluster_gbk(gbkFile, FileType, locustag_product_monomer_dict):
def get_cluster_info_from_cluster_gbk(gbkFile, FileType):

    cluster_info_dict = {}
 
    #Read GenBank file
    record = SeqIO.read(gbkFile, FileType)
    
    for feature in record.features:
        sec_met_info_list = []
        if feature.type == 'CDS':
            #MIGHT NEED TO ADD 'function' for the qualifier
            #Prevent extracting genes without specific information
            #Necessary not to have unidentified monomers in subsequent functions, cuasing error
            if 'hypothetical' not in feature.qualifiers.get('product')[0]: 
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                if feature.qualifiers.get('sec_met'):
                
                    qualifier_sec_met = feature.qualifiers.get('sec_met')
                    #print qualifier_sec_et
                    sec_met_info_list.append(qualifier_sec_met)

                    cluster_info_dict[qualifier_locus_tag] = sec_met_info_list

    #print 'cluster_info_dict'
    #print cluster_info_dict, '\n'
    for i in cluster_info_dict.keys():
        print i
        print cluster_info_dict[i], '\n'
    return cluster_info_dict


#Output: e.g.
#Cluster number: 2
#Product: nrps
#NC021055_Cluster_02_nrps
def get_product_from_cluster_gbk(gbkFile, FileType):
 
    #Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)
    
    for feature in record.features:

        #Retrieving "Cluster number"
        if feature.type == 'cluster':

            qualifier_cluster = feature.qualifiers.get('note')
            qualifier_cluster = qualifier_cluster[0].split(':')
            clusterNo = qualifier_cluster[1].strip()

            #Retrieving "product"
            product = feature.qualifiers.get('product')
            product = product[0]
            
    if float(clusterNo) < 10:
        product = "Cluster_0"+clusterNo+"_"+product
    else:
        product = "Cluster_"+clusterNo+"_"+product

    print product, "\n"
    return product


#Output: e.g.,
#['SAV_938'] = ['PKS_AT', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_KR', 'ACP', 'PKS_KS']
def get_cluster_domain(cluster_info_dict):
    
    locustag_domain_dict = {}
        
    for each_gene in cluster_info_dict.keys():
        
        sec_met_info_list = cluster_info_dict[each_gene][0]

        domain_count = 0
        sec_met_domain_list = []

        for each_sub_set in sec_met_info_list:

            if "NRPS/PKS Domain" in each_sub_set:
                sptline1 = each_sub_set.split('. ')
                crude_domain_info = sptline1[0] 

                #Extract the information of KR activity and AT substrate specificity
                sptline2 = each_sub_set.split('; ')                

                spt_domain_info = crude_domain_info.split(':')
                whole_domain_info = spt_domain_info[1]

                spt_list_domain_info = whole_domain_info.split()
                spt_list_domain_info.append(each_gene)

                each_sec_met_domain = spt_list_domain_info[0]
                sec_met_domain_list.append(each_sec_met_domain)
                
                #Domain information
                domain_number = each_gene + '_DM' + str(domain_count)
                domain_count = domain_count + 1
                            
        locustag_domain_dict[each_gene] = sec_met_domain_list

    print 'locustag_domain_dict'
    print locustag_domain_dict, '\n'

    return locustag_domain_dict


#Two nested functions: extract_nrp_monomers, extract_pk_monomers
#Output: e.g., {'SAV_943_M1':['mmal', 'Ethyl_mal', 'pk']}
def get_cluster_monomers(cluster_info_dict):

    locustag_monomer_dict = {}
    for each_gene in cluster_info_dict.keys():
        module_count = 0
        sec_met_info_list =  cluster_info_dict[each_gene][0]

        for each_sub_set in sec_met_info_list:
            discriminator = "true"

            if "Substrate specificity predictions" in each_sub_set and "AMP-binding" in each_sub_set:
                pred_monomer_list, discriminator = extract_nrp_monomers(each_sub_set, discriminator)
                module_number = each_gene + '_M' + str(module_count)
                locustag_monomer_dict[module_number] = pred_monomer_list
                module_count += 1
                #print "check", module_number, pred_monomer_list

            if "Substrate specificity predictions" in each_sub_set and "A-OX" in each_sub_set:
                pred_monomer_list, discriminator = extract_nrp_monomers(each_sub_set, discriminator)
                module_number = each_gene + '_M' + str(module_count)
                locustag_monomer_dict[module_number] = pred_monomer_list
                module_count += 1
                #print "check", module_number, pred_monomer_list

            if "Substrate specificity predictions" in each_sub_set and "PKS_AT" in each_sub_set:
                pred_monomer_list = extract_pk_monomers(each_sub_set)
                module_number = each_gene + '_M' + str(module_count)
                locustag_monomer_dict[module_number] = pred_monomer_list
                module_count += 1
                print "check", module_number, pred_monomer_list

            if discriminator == "false":
                continue
 
    print 'locustag_monomer_dict'
    #print locustag_monomer_dict, '\n'
    return locustag_monomer_dict


def extract_nrp_monomers(each_sub_set, discriminator):
    sptline2 = each_sub_set.split(';')
    whole_substrate_info = sptline2[1]
    pred_monomer_list = []

    predicted_monomers = whole_substrate_info.split(':')
    sptSubstrates = predicted_monomers[1]

    if ', ' not in sptSubstrates:
        print "Insufficient substrate_info"
        discriminator = "false"
        return pred_monomer_list

    substrates = sptSubstrates.split(', ')

    for each_substrate in substrates:
        sptSubstrate = each_substrate.split('(')
        predicted_monomer = sptSubstrate[0].strip()

        pred_monomer_list.append(predicted_monomer)

    return pred_monomer_list, discriminator


def extract_pk_monomers(each_sub_set):
    sptline2 = each_sub_set.split(';')
    whole_substrate_info = sptline2[1]

    predicted_monomers = whole_substrate_info.split(':')
    sptSubstrates = predicted_monomers[1]
    substrates = sptSubstrates.split(', ')

    pred_monomer_list = []

    for each_substrate in substrates:
        sptSubstrate = each_substrate.split('(')
        predicted_monomer = sptSubstrate[0].strip()

        pred_monomer_list.append(predicted_monomer)

    return pred_monomer_list


#Output: e.g., {'SAV_943_M1': ['PKS_KS', 'PKS_AT', 'ACP']}
def get_cluster_module(locustag_domain_dict):
    
    locustag_module_domain_dict = {}

    for locustag in locustag_domain_dict.keys():
        
        count = 0
        list_module_info = []
        number_of_list = len(locustag_domain_dict[locustag])

        for each_domain in locustag_domain_dict[locustag]:
            list_module_info.append(each_domain)
            number_of_list -= 1
        
            if each_domain == 'PCP' or each_domain == 'ACP':
                module_number = locustag + '_M' + str(count)
                locustag_module_domain_dict[module_number] = list_module_info

                list_module_info = []
                count += 1

            elif each_domain == 'Epimerization':
                count -= 1
                module_number = locustag + '_M' + str(count)
                list_module_info = locustag_module_domain_dict[module_number]
                list_module_info.append('Epimerization')
                locustag_module_domain_dict[module_number] = list_module_info

                list_module_info = []
                count += 1

            #'ACP' was inserted, otherwise it causes an error
            #by having a locus_tag with unspecified monomer
            elif each_domain == 'Thioesterase' and 'ACP' in locustag_domain_dict[locustag]:
                count -= 1
                module_number = locustag + '_M' + str(count)
                list_module_info = locustag_module_domain_dict[module_number]
                list_module_info.append('Thioesterase')
                locustag_module_domain_dict[module_number] = list_module_info

            elif list_module_info.count('PKS_KS') == 2:
                module_number = t1pks_gene + '_M' + str(count)
                poped_domain = list_module_info.pop()
                dic_pksnrps_module[module_number] = list_module_info

                list_module_info = []
                list_module_info.append(poped_domain)
                count += 1

            elif list_module_info.count('Condensation_DCL') == 2 or list_module_info.count('Condensation_LCL') == 2 or list_module_info.count('Condensation_LCL') + list_module_info.count('Condensation_DCL') == 2:
                module_number = locustag + '_M' + str(count)
                list_module_info.pop()
                locustag_module_domain_dict[module_number] = list_module_info

                list_module_info = []
                count += 1

            elif float(number_of_list) == 0:
                module_number = locustag + '_M' + str(count)
                locustag_module_domain_dict[module_number] = list_module_info

                list_module_info = []
                count += 1

    print 'locustag_module_domain_dict'
    print locustag_module_domain_dict, '\n'
    return locustag_module_domain_dict


#Ouput: e.g., {'SAV_943_M0':{'coa': 1, 'nadph': -1, 'nadp': 1, 'hco3': 1, 'h': -1}
def get_currency_metabolites(locustag_module_domain_dict):

    module_currency_metab_dict = {}

    for each_module in locustag_module_domain_dict:
        domain_comb = locustag_module_domain_dict[each_module]

        each_module_substrates = {}
        discriminant = determine_module(domain_comb)

        if discriminant == 'None':
            #print "Discriminant not defined : %s" % (domain_comb)
            continue

        if discriminant == 'A' or discriminant == 'Aox':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['h2o'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates 
            #print 'A or Aox:', each_module_substrates 

        elif discriminant == 'A_PCP' or discriminant == 'Cs_A_PCP' or discriminant == 'C_A_PCP' or discriminant == 'Cdcl_A_PCP' or discriminant == 'Clcl_A_PCP' or discriminant == 'Clcl_A_PCP' or discriminant == 'Cglyc_A_PCP' or discriminant == 'CXglyc_A_PCP':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['h2o'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'A_PCP or C_A_PCP:', each_module_substrates

        elif discriminant == 'MT_A_PCP' or discriminant == 'Cs_A_MT_PCP' or discriminant == 'C_A_MT_PCP' or discriminant == 'Cdcl_A_MT_PCP' or discriminant == 'Clcl_A_MT_PCP' or discriminant == 'Cglyc_A_MT_PCP' or discriminant == 'CXglyc_A_MT_PCP':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            each_module_substrates['h2o'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'MT-A-PCP or C-A-MT-PCP:', each_module_substrates

        elif discriminant == 'Cs_A_PCP_E' or discriminant == 'C_A_PCP_E' or discriminant == 'Cdcl_A_PCP_E' or discriminant == 'Clcl_A_PCP_E' or discriminant == 'Cd_A_PCP' or discriminant == 'Cglyc_A_PCP_E' or discriminant == 'CXglyc_A_PCP_E':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['h2o'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'C_A_PCP_E:', each_module_substrates

        elif discriminant == 'Cs_A_MT_E_PCP' or discriminant == 'C_A_MT_E_PCP' or discriminant == 'Cdcl_A_MT_E_PCP' or discriminant == 'Clcl_A_MT_E_PCP' or discriminant == 'Cglyc_A_MT_E_PCP' or discriminant == 'CXglyc_A_MT_E_PCP' or discriminant == 'Cd_A_MT_PCP':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            each_module_substrates['h2o'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'C_A_MT_E_PCP:', each_module_substrates

        elif discriminant == 'HC_A_PCP' :
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['h2o'] = 2
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'HC_A_PCP:', each_module_substrates

        elif discriminant == 'HC_Aox_PCP':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            each_module_substrates['fmn'] = 1
            each_module_substrates['fmnh2'] = -1
            each_module_substrates['h2o'] = 2
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'HC_Aox_PCP:', each_module_substrates

        elif discriminant == 'HC_A_MT_PCP' or discriminant == 'HC_A_MT_E_PCP':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['ppi'] = 1
            each_module_substrates['h2o'] = 2
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'HC_A_MT_PCP:', each_module_substrates

        elif discriminant == 'HC_Aox_MT_PCP' or discriminant == 'HC_Aox_MT_E_PCP':
            each_module_substrates['atp'] = -1
            each_module_substrates['amp'] = 1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            each_module_substrates['fmn'] = 1
            each_module_substrates['fmnh2'] = -1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            each_module_substrates['h2o'] = 2
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'HC_Aox_MT_PCP:', each_module_substrates

        elif discriminant == 'AT_ACP' or discriminant == 'AT_KR(inactive)_ACP' or discriminant == 'AT_DH_KR(inactive)_ACP' or discriminant == 'AT_DH_ER_KR(inactive)_ACP':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'AT-ACP:', each_module_substrates

        elif discriminant == 'AT_KS_ACP' or discriminant == 'AT_KS' or discriminant == 'AT_KS_KR(inactive)_ACP' or discriminant == 'AT_KS_KR(inactive)' or discriminant == 'AT_KS_DH_KR(inactive)_ACP' or discriminant == 'AT_KS_DH_KR(inactive)' or discriminant == 'AT_KS_DH_KR(inactive)_ACP' or discriminant == 'AT_KS_DH_ER_KR(inactive)':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'AT-KS-ACP:', each_module_substrates

        elif discriminant == 'AT_KS_KR_ACP' or discriminant == 'AT_KR_ACP' or discriminant == 'AT_KS_KR' or discriminant == 'AT_ER_KR_ACP' or discriminant == 'AT_KS_ER_KR' or discriminant == 'AT_KS_ER_KR_ACP':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            each_module_substrates['nadp'] = 1
            each_module_substrates['nadph'] = -1
            each_module_substrates['h'] = -1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'KS-AT-KR-ACP:', each_module_substrates

        elif discriminant == 'AT_KS_DH_KR_ACP' or discriminant == 'AT_DH_KR_ACP' or discriminant == 'AT_KS_DH_KR':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            each_module_substrates['h2o'] = 1
            each_module_substrates['nadp'] = 1
            each_module_substrates['nadph'] = -1
            each_module_substrates['h'] = -1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'KS-AT-DH-KR-ACP:', each_module_substrates

        elif discriminant == 'AT_KS_DH_ER_KR_ACP' or discriminant == 'AT_DH_ER_KR_ACP' or discriminant == 'AT_KS_DH_ER_KR':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            each_module_substrates['h2o'] = 1
            each_module_substrates['nadp'] = 2
            each_module_substrates['nadph'] = -2
            each_module_substrates['h'] = -2
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'KS-AT-DH-ER-KR-ACP:', each_module_substrates

        elif discriminant == 'AT_KS_cMT_ACP' or discriminant == 'AT_KS_cMT':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'AT_KS_cMT_ACP:', each_module_substrates

        elif discriminant == 'AT_KS_KR_cMT_ACP' or discriminant == 'AT_KS_KR_cMT':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            each_module_substrates['nadp'] = 1
            each_module_substrates['nadph'] = -1
            each_module_substrates['h'] = -1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'AT_KS_KR_cMT_ACP:', each_module_substrates

        elif discriminant == 'AT_KS_DH_KR_cMT_ACP' or discriminant == 'AT_KS_DH_KR_cMT':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            each_module_substrates['h2o'] = 1
            each_module_substrates['nadp'] = 1
            each_module_substrates['nadph'] = -1
            each_module_substrates['h'] = -1
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'AT_KS_DH_KR_cMT_ACP:', each_module_substrates

        elif discriminant == 'AT_KS_DH_KR_cMT_ER_ACP' or discriminant == 'AT_KS_DH_KR_cMT_ER':
            each_module_substrates['coa'] = 1
            each_module_substrates['hco3'] = 1
            each_module_substrates['h2o'] = 1
            each_module_substrates['nadp'] = 2
            each_module_substrates['nadph'] = -2
            each_module_substrates['h'] = -2
            each_module_substrates['amet'] = -1
            each_module_substrates['ahcys'] = 1
            module_currency_metab_dict[each_module] = each_module_substrates
            #print 'AT_KS_DH_KR_cMT_ER_ACP:', each_module_substrates

        elif discriminant == 'TD':
            each_module_substrates['nadp'] = 1
            each_module_substrates['nadph'] = -1
            each_module_substrates['h'] = -1
            module_currency_metab_dict[each_module] = each_module_substrates                
            #print 'reaction 8: HC_Aox_MT_PCP', each_module_substrates 
            
        elif discriminant == 'ACP':
            continue
        elif discriminant == 'PCP':
            continue

    print "module_currency_metab_dict"
    print module_currency_metab_dict, '\n'
    return module_currency_metab_dict


#Output: {'nadph': 0, 'fmnh2': 0, 'h': 0, 'ppi': 1, 'ahcys': 0}
def get_total_currency_metab_coeff(module_currency_metab_dict):

    currency_metab_coeff_dict = get_metab_coeff_dict()

    for each_module in module_currency_metab_dict.keys():
        for each_metabolite in module_currency_metab_dict[each_module]:
            metab_coeff = module_currency_metab_dict[each_module][each_metabolite]

            if module_currency_metab_dict[each_module][each_metabolite] > 0:
                currency_metab_coeff_dict[each_metabolite] += metab_coeff

            else:
                currency_metab_coeff_dict[each_metabolite] += metab_coeff

    print 'currency_metab_coeff_dict' 
    print currency_metab_coeff_dict, '\n'
    return currency_metab_coeff_dict


#Coeff data of major monomers are added in the same dict file used for currency metabolites
#Output: e.g.,
#{'coa': 13, 'mmalcoa': -4, 'h': -10, 'malcoa': -7,     'hco3': 13, 'nadph': -10, 'h2o': 5, 'nadp': 10}
def get_all_metab_coeff(locustag_monomer_dict, metab_coeff_dict, product):
    
    for each_module in locustag_monomer_dict.keys():
        #locustag_monomer_dict[each_module] for nrps
        #Position [0]: NRPSPredictor2 SVM
        #Position [1]: Stachelhaus code
        #Position [2]: Minowa
        #Position [3]: consensus
        if len(locustag_monomer_dict[each_module]) == 4:

            sptlist1 = locustag_monomer_dict[each_module][0].split(',')
            #print "CHECK", sptlist1, len(sptlist1)

            #In case "consensus" is not reached:
            if locustag_monomer_dict[each_module][3] == 'nrp':
                #From NRPSPredictor2 SVM 
                aSid_met2 = locustag_monomer_dict[each_module][0]
                biggid_met2 = get_biggid_from_aSid(aSid_met2)
                #print "aSid_met2", aSid_met2, biggid_met2

                #In case of non-consensus, NRPSPredictor2 SVM is considered 
                metab_coeff_dict[biggid_met2] -= 1

            #In case "consensus" is reached:
            elif locustag_monomer_dict[each_module][3] != 'nrp':
                aSid_met5 = locustag_monomer_dict[each_module][3]
                biggid_met5 = get_biggid_from_aSid(aSid_met5)
                #print "aSid_met5", aSid_met5, biggid_met5
                metab_coeff_dict[biggid_met5] -= 1

        #locustag_monomer_dict[each_module] for pks
        #Position [0]: PKS signature
        #Position [1]: Minowa
        #Position [2]: consensus
        elif len(locustag_monomer_dict[each_module]) == 3:

            if len(locustag_monomer_dict[each_module]) < 3:
                continue

            #In case "consensus" is not reached:
            if locustag_monomer_dict[each_module][2] == 'pk':

                #From PKS signature 
                aSid_met6 = locustag_monomer_dict[each_module][0]
                biggid_met6 = get_biggid_from_aSid(aSid_met6)
                #print "aSid_met6", aSid_met6, biggid_met6

                #In case of non-consensus, PKS signature is considered
                metab_coeff_dict[biggid_met6] -= 1

            #In case "consensus" is reached:
            elif locustag_monomer_dict[each_module][2] != 'pk':
                aSid_met8 = locustag_monomer_dict[each_module][2]
                biggid_met8 = get_biggid_from_aSid(aSid_met8)
                #print "aSid_met8", aSid_met8, biggid_met8
                metab_coeff_dict[biggid_met8] -= 1

    #Add secondary metabolite product to the reaction
    metab_coeff_dict[product] = 1

    print 'metab_coeff_dict'
    print metab_coeff_dict, '\n'
    return metab_coeff_dict


def converting_MNXMID_to_biggid(MetID):
    converted_MNXMID = []
     
    if MetID == 'ala_DASH_L':
        converted_MNXMID = ['MNXM32', 'ala_DASH_L', 'L-alanine', 'C00041']
    elif MetID == 'arg_DASH_L':
        converted_MNXMID = ['MNXM70', 'arg_DASH_L', 'L-Arginine', 'C00062']
    elif MetID == 'asn_DASH_L':
        converted_MNXMID = ['MNXM147', 'asn_DASH_L', 'L-asparagine', 'C00152'] 
    elif MetID == 'asp_DASH_L':
        converted_MNXMID = ['MNXM42', 'asp_DASH_L', 'L-Aspartate' ,'C00049']
    elif MetID == 'cys_DASH_L':
        converted_MNXMID = ['MNXM55', 'cys_DASH_L', 'L-Cysteine', 'C00097']
    elif MetID == 'gln_DASH_L': 
        converted_MNXMID = ['MNXM37', 'gln_DASH_L', 'L-Glutamine', 'C00064']
    elif MetID == 'glu_DASH_L':
        converted_MNXMID = ['MNXM89557', 'glu_DASH_L', 'L-Glutamate', 'C00025']
    elif MetID == 'gly':
        converted_MNXMID = ['MNXM29', 'gly', 'glycine', 'C00037']
    elif MetID == 'his_DASH_L':
        converted_MNXMID = ['MNXM569', 'his_DASH_L', 'L-Histidine', 'C01033']
    elif MetID == 'leu_DASH_L':
        converted_MNXMID = ['MNXM134', 'leu_DASH_L', 'L-Leucine', 'C00135']
    elif MetID == 'lys_DASH_L':
        converted_MNXMID = ['MNXM78', 'lys_DASH_L', 'L-Lysine', 'C00047']
    elif MetID == 'met_DASH_L':
        converted_MNXMID = ['MNXM61', 'met_DASH_L', 'L-Methionine', 'C00073']
    elif MetID == 'phe_DASH_L':
        converted_MNXMID = ['MNXM97', 'phe_DASH_L', 'L-Phenylalanine', 'C00079']
    elif MetID == 'pro_DASH_L':
        converted_MNXMID = ['MNXM114', 'pro_DASH_L', 'L-Proline', 'C00148']
    elif MetID == 'ser_DASH_L':
        converted_MNXMID = ['MNXM53', 'ser_DASH_L', 'L-Serine', 'C00065']
    elif MetID == 'thr_DASH_L':
        converted_MNXMID = ['MNXM142', 'thr_DASH_L', 'L-Threonine', 'C00188']
    elif MetID == 'trp_DASH_L':
        converted_MNXMID = ['MNXM94', 'trp_DASH_L', 'L-Tryptophan', 'C00078']
    elif MetID == 'tyr_DASH_L':
        converted_MNXMID = ['MNXM76', 'tyr_DASH_L', 'L-Tyrosine', 'C00082']
    elif MetID == 'ile_DASH_L':
        converted_MNXMID = ['MNXM231', 'ile_DASH_L', 'L-Isoleucine', 'C00407']
    elif MetID == 'phg_DASH_L':    
        converted_MNXMID = ['MNXM59292', 'phg_DASH_L', 'L-phenylglycine', 'C18623']
    elif MetID == 'bht_DASH_L':    
        converted_MNXMID = ['N/A', 'bht_DASH_L', 'beta-hydroxy-tyrocine', 'N/A']
    elif MetID == 'orn':    
        converted_MNXMID = ['MNXM89689', 'orn', 'ornitine', 'C01602']
    elif MetID == 'abu':    
        converted_MNXMID = ['MNXM17054', 'abu', 'D-2-Aminobutyric acid', 'C02261']
    elif MetID == 'iva':
        converted_MNXMID = ['MNXM34821', 'iva', '2-Amino-2-methylbutanoate', 'C03571'] 
    elif MetID == 'L2aadp':
        converted_MNXMID = ['MNXM268', 'L2aadp', 'L-2-Aminoadipic acid', 'C00956'] 
    elif MetID == 'hpg':
        converted_MNXMID = ['N/A', 'hpg', 'D-4-Hydroxyphenylglycine', 'C03493'] 
    elif MetID == '23dhb':
        converted_MNXMID = ['MNXM455', '23dhb', '2,3-Dihydroxybenzoic acid', 'C00196']
    elif MetID == 'dhpg':
        converted_MNXMID = ['MNXM9962', 'dhpg', '3,5-Dihydroxy-phenylglycine', 'C12026']
    elif MetID == 'hty':
        converted_MNXMID = ['MNXM59438', 'hty', 'L-Homotyrosine', 'C18622']
    elif MetID == 'citr_DASH_L':
        converted_MNXMID = ['MNXM211', 'citr_DASH_L', 'L-Citruline', 'C00327']       
    elif MetID == 'Lpipecol':
        converted_MNXMID = ['MNXM684', 'Lpipecol', 'L-pipecolate', 'C00408']       
    elif MetID == 'ala_DASH_B':
        converted_MNXMID = ['MNXM144', 'ala_DASH_B', 'beta-Alanine', 'C00099']  
    elif MetID == '24dab':
        converted_MNXMID = ['MNXM840', '24dab', 'L-2,4-diazaniumylbutyrate', 'C03283']
    elif MetID == 'pac':
        converted_MNXMID = ['NXM497', 'pac', 'phenylacetate', 'C00548']  
    elif MetID == 'tcl':
        converted_MNXMID = ['MNXM37380', 'tcl', '4-Chlorothreonine', 'N/A'] 
    elif MetID == 'qa':
        converted_MNXMID = ['MNXM4797', 'qa', 'Quinaldinic acid', 'C06325']

    elif MetID == 'malcoa':
        converted_MNXMID = ['MNXM40', 'malcoa', 'malonyl-CoA', 'C00083']
    elif MetID == 'mmcoa_DASH_S':
        converted_MNXMID = ['MNXM89955', 'mmcoa_DASH_S', 'methylmaloncyl-CoA', 'C02557']
    elif MetID == '2mbcoa':
        converted_MNXMID = ['MNXM569', '2mbcoa', '2-methylbutanoyl-CoA', 'C01033']
    elif MetID == 'emcoa_DASH_S':
        converted_MNXMID = ['MNXM2043', 'emcoa_DASH_S', 'ethylmalonyl-CoA', 'C18026']
    elif MetID == 'ibcoa':
        converted_MNXMID = ['MNXM470', 'ibcoa', '2-Methylpropanoyl-CoA', 'C00630']
    elif MetID == 'accoa':
        converted_MNXMID = ['MNXM21', 'accoa', 'Acetyl-CoA', 'C00024']
    elif MetID == 'ppcoa':
        converted_MNXMID = ['MNXM86', 'ppcoa', 'Propionyl-CoA', 'C00100']
    elif MetID == 'ivcoa':
        converted_MNXMID = ['MNXM471', 'ivcoa', '3-Methylbutanoyl-CoA', 'C02939']
    elif MetID == 'mxmalacp':
        converted_MNXMID = ['MNXM61686', 'mxmalacp', 'Methoxymalonyl-[acp]', 'C18616']
    elif MetID == 'chccoa':
        converted_MNXMID = ['MNXM5111', 'chccoa', 'cyclohexane-1-carboxyl-CoA', 'C09823']
 
    elif MetID == 'nadp':
        converted_MNXMID = ['MNXM5', 'nadp', 'NADP+', 'C00006'] 
    elif MetID == 'nadph':
        converted_MNXMID = ['MNXM6', 'nadph', 'NADPH' ,'C00005']
    elif MetID == 'h':
        converted_MNXMID = ['MNXM1', 'h', 'H+', 'C00080']     
    elif MetID == 'atp':
        converted_MNXMID = ['MNXM3', 'atp', 'Adenosine-5-prime-triphosphate', 'C00002']
    elif MetID == 'amp':
        converted_MNXMID = ['MNXM14', 'amp', 'Adenosine-5-prime-phosphate', 'C00020']    
    elif MetID == 'ppi':
        converted_MNXMID = ['MNXM11', 'ppi', 'diphosphate', 'C00013'] 
    elif MetID == 'amet':
        converted_MNXMID = ['MNXM16', 'amet', 'S-adenosyl-L-methionine', 'C00019']      
    elif MetID == 'ahcys':
        converted_MNXMID = ['MNXM19', 'ahcys', 'S-adenosyl-L-homocysteine', 'C00021']
    elif MetID == 'fmn':
        converted_MNXMID = ['MNXM119', 'fmn', 'Riboflavin-5-phosphate', 'C00061']
    elif MetID == 'fmnh2':
        converted_MNXMID = ['MNXM208', 'fmnh2', 'Reduced Riboflavin-5-phosphate', 'C01847'] 
    elif MetID == 'h2o':
        converted_MNXMID = ['MNXM2', 'h2o', 'water', 'C00001']                
    else:
        converted_MNXMID = ['', MetID, '', '']
     
    return converted_MNXMID

def add_sec_met_rxn(cobra_model, product, locustag_product_monomer_dict, metab_coeff_dict):
    product_count = 1
    
    list_reaction_name_SM = []
    list_novel_secondary_metabolite_reactions = []
    
    new_product_name = product + '_' +  str(product_count)
    list_reaction_name_SM.append(new_product_name)
    list_novel_secondary_metabolite_reactions.append(each_integrated_reaction)

    product_count += 1

#Creating reaction ID
    reaction = Reaction(new_product_name)

#Setting bounds
    reaction.lower_bound = 0
    reaction.upper_bound = 999999

#Adding substrate metabolites
    for each_metabolite in each_integrated_reaction:

        converted_MNXMID = converting_MNXMID_to_biggid(each_metabolite)

        Met_MNXMID = converted_MNXMID[0]
        abbr_MetID = converted_MNXMID[1]
        Met_Name = converted_MNXMID[2]
        Met_KEGG_ID = converted_MNXMID[3]

        if Met_MNXMID in metab_MNXM_dict and metab_MNXM_dict[Met_MNXMID][-1:] == 'c':
            reaction.add_metabolites({cobra_model.metabolites.get_by_id(metab_MNXM_dict[Met_MNXMID]):each_integrated_reaction[abbr_MetID]})
        elif each_integrated_reaction[abbr_MetID] == 0:
            continue
        else:
            obj_each_metabolite = Metabolite(abbr_MetID+'_c', name=Met_Name, compartment='c')
            reaction.add_metabolites({obj_each_metabolite:each_integrated_reaction[abbr_MetID]})

#Setting GPR association
    gpr_count = 0
    for each_gene in locustag_product_monomer_dict:
        if gpr_count == 0:
            gpr_list = each_gene 
            gpr_count += 1
        else:
            gpr_list = gpr_list + ' AND ' + each_gene
     
    print gpr_list
    reaction.add_gene_reaction_rule(gpr_list)

#Adding the new reaction to the model
    cobra_model.add_reaction(reaction)

    reaction_name = reaction.id
    strain_name = reaction_name.split("_")
    strain_name = strain_name[0].strip()

    print "\n", "Cluster reaction:", reaction
    print "Cluster genes:", reaction.gene_reaction_rule
    print reaction.reaction

#Creating a transport reaction
#Creating reaction ID
    reaction = Reaction("Transport_" + new_product_name )

#Setting bounds
    reaction.reversibility = 0 # 1: reversible
    reaction.lower_bound = 0
    reaction.upper_bound = 999999

#Adding a substrate metabolite
#    print cobra_model.metabolites.get_by_id(str(product_c))
    reaction.add_metabolites({cobra_model.metabolites.get_by_id(str(new_product_name+'_c')):-1})

#Adding product metabolite(s)
    new_product_name_e = Metabolite(new_product_name+"_e", name='', compartment='e')
    reaction.add_metabolites({new_product_name_e:1})

#Adding the new reaction to the model
    cobra_model.add_reaction(reaction)

    print "\n", "Transport reaction:", reaction
    print reaction.reaction

#Creating an exchange reaction
#Creating reaction ID
    reaction = Reaction("Ex_"+new_product_name)

#Setting bounds
    reaction.reversibility = 0 # 1: reversible 0: irreversible
    reaction.lower_bound = 0
    reaction.upper_bound = 999999

#Adding a substrate metabolite
#    print cobra_model.metabolites.get_by_id(str(product_c))
    reaction.add_metabolites({cobra_model.metabolites.get_by_id(str(new_product_name_e)):-1})


#Adding the new reaction to the model
    cobra_model.add_reaction(reaction)

    print "\n", "Exchange reaction:", reaction
    print reaction.reaction

