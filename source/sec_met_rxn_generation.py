'''
2015 Kyu-Sang Hwang
2014-2015 Hyun Uk Kim

This file generates metabolic reactions for the genes newly annotated to be present in the secondary metabolite-biosynthetic gene cluster from antiSMASH.
'''

from Bio import SeqIO
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from general_sec_met_info import (
    determine_module,
    determine_kr_activity,
    get_module_currency_metab_dict,
    get_biggid_from_aSid,
    get_metab_coeff_dict,
    add_sec_met_mnxm_having_no_biggid_to_model
)


#Exracts all the information associated wiht a particular locus_tag
#def get_locustag_info_from_cluster_gbk(gbkFile, FileType, locustag_product_monomer_dict):
def get_cluster_info_from_cluster_gbk(gbkFile, FileType):

    cluster_info_dict = {}
 
    #Read GenBank file
    record = SeqIO.read(gbkFile, FileType)
    
    for feature in record.features:
        sec_met_info_list = []
        if feature.type == 'CDS':
            qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
            if feature.qualifiers.get('sec_met'):
                qualifier_sec_met = feature.qualifiers.get('sec_met')
                sec_met_info_list.append(qualifier_sec_met)

                cluster_info_dict[qualifier_locus_tag] = sec_met_info_list

    #print 'cluster_info_dict'
    #print cluster_info_dict, '\n'
    return cluster_info_dict, record


#Output: e.g.
#Cluster number: 2
#Product: nrps
#NC021055_Cluster_02_nrps
def get_product_from_cluster_gbk(record):
 
    for feature in record.features:

        #Retrieving "Cluster number"
        if feature.type == 'cluster':

            qualifier_cluster = feature.qualifiers.get('note')
            qualifier_cluster = qualifier_cluster[0].split(':')
            clusterNo = qualifier_cluster[1].strip()

            #Retrieving "product"
            product = feature.qualifiers.get('product')
            product = product[0]

    #Handle legacy problem
    product = product.replace("-","_")

    if float(clusterNo) < 10:
        product = "Cluster0"+clusterNo+"_"+product
    else:
        product = "Cluster"+clusterNo+"_"+product

    #print product, "\n"
    return product


#Output: e.g.,
#['SAV_938'] = ['PKS_AT', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_KR', 'ACP', 'PKS_KS']
def get_cluster_domain(cluster_info_dict):
    
    locustag_domain_dict = {}
    locustag_kr_dict = {}
    
    for each_gene in cluster_info_dict.keys():
        
        sec_met_info_list = cluster_info_dict[each_gene][0]

        domain_count = 0
        sec_met_domain_list = []
        kr_domain_info_dict = {}

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

                #Domain information
                #Correction by HKS; counts the number of domains
                if float(domain_count) < 10:
                    domain_number = "_DM0"+str(domain_count)
                else:
                    domain_number = "_DM"+str(domain_count)

                domain_number = each_sec_met_domain + domain_number
                sec_met_domain_list.append(domain_number)
                domain_count += 1

                #Collects KR activity information
                if domain_number[:-5] == "PKS_KR":   
                    kr_info_list = sptline2[1].split(': ')
                    dm_kr_activity = kr_info_list[1]
                    kr_domain_info_dict[domain_number] = dm_kr_activity
 
        locustag_domain_dict[each_gene] = sec_met_domain_list
        locustag_kr_dict[each_gene] = kr_domain_info_dict

    #print 'locustag_domain_dict'
    #print locustag_domain_dict, '\n'

    return locustag_domain_dict, locustag_kr_dict


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
                #print "check", module_number, pred_monomer_list

            if discriminator == "false":
                continue
 
    #print 'locustag_monomer_dict'
    #print locustag_monomer_dict, '\n'
    return locustag_monomer_dict


def extract_nrp_monomers(each_sub_set, discriminator):
    sptline2 = each_sub_set.split(';')

    for element in sptline2:
        if 'Substrate specificity predictions' in element:
            pred_monomer_list = []

            predicted_monomers = element.split(':')
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
        module_info_list = []
        module_info_trunc_list = []
        number_of_list = len(locustag_domain_dict[locustag])

        for each_domain in locustag_domain_dict[locustag]:
            domain_name = each_domain[:-5]
            module_info_list.append(each_domain)
            module_info_trunc_list.append(domain_name)
            number_of_list -= 1

            #Check final domain of a module: core domain
            if domain_name == 'PCP' or domain_name == 'ACP':
                module_number = get_locustag_module_number(locustag, count)
                locustag_module_domain_dict[module_number] = module_info_list

                module_info_list = []
                module_info_trunc_list = []
                count += 1

            #Check final domain of a module: optional domain
            #Presence of 'PCP' was considered in case Epimerization appears first before 'PCP'
            #e.g., 'M271_46685' in cluster 46 of 'src' (KEGG organism code)
            elif domain_name == 'Epimerization' and 'PCP' in module_info_list:
                count -= 1

                module_number = get_locustag_module_number(locustag, count)
                module_info_list = locustag_module_domain_dict[module_number]
                module_info_list.append(str(each_domain))
                locustag_module_domain_dict[module_number] = module_info_list

                module_info_list = []
                module_info_trunc_list = []
                count += 1

            #Check final domain of a module: linker domain
            elif domain_name == 'PKS_Docking_Cterm' and 'ACP' in locustag_domain_dict[locustag]:
                count -= 1
                module_number = locustag + '_M' + str(count)
                module_info_list = locustag_module_domain_dict[module_number]
                module_info_list.append('PKS_Docking_Cterm')
                locustag_module_domain_dict[module_number] = module_info_list

                module_info_list = []
                module_info_trunc_list = []
                count += 1

            #Check final domain of a module: final point of the carbon scaffold
            #'ACP' was inserted, otherwise it causes an error
            #by having a locus_tag with unspecified monomer
            elif domain_name == 'Thioesterase' and 'ACP' in locustag_domain_dict[locustag]:
                count -= 1
                module_number = get_locustag_module_number(locustag, count)
                module_info_list = locustag_module_domain_dict[module_number]
                module_info_list.append(str(each_domain))
                locustag_module_domain_dict[module_number] = module_info_list

            elif module_info_list.count('PKS_KS') == 2 or module_info_list.count('Condensation') == 2 or module_info_list.count('Condensation_DCL') == 2 or module_info_list.count('Condensation_LCL') == 2 or module_info_list.count('Condensation_LCL') + module_info_list.count('Condensation_DCL') == 2 or module_info_list.count('Condensation_Dual') == 2 or module_info_list.count('Cglyc') == 2 or module_info_list.count('CXglyc') == 2 or module_info_list.count('Heterocyclization') == 2:
                module_number = get_locustag_module_number(locustag, count)
                poped_domain = module_info_list.pop()
                abbr_poped_domain = module_info_trunc_list.pop()
                locustag_module_domain_dict[module_number] = module_info_list

                module_info_list = []
                module_info_trunc_list = []
                module_info_list.append(poped_domain)
                module_info_trunc_list.append(abbr_poped_domain)
                count += 1

            elif float(number_of_list) == 0:
                module_number = get_locustag_module_number(locustag, count)
                locustag_module_domain_dict[module_number] = module_info_list

                module_info_list = []
                module_info_trunc_list = []
                count += 1

    #print 'locustag_module_domain_dict'
    #print locustag_module_domain_dict, '\n'
    return locustag_module_domain_dict


def get_locustag_module_number(locustag, count):
    if float(count) < 10:
        module_number = locustag + '_M0' + str(count)
    else:
        module_number = locustag + '_M' + str(count)

    return module_number


#Ouput: e.g., {'SAV_943_M0':{'coa': 1, 'nadph': -1, 'nadp': 1, 'hco3': 1, 'h': -1}
def get_currency_metabolites(locustag_module_domain_dict, locustag_kr_dict):

    module_currency_metab_dict = {}

    for each_module in locustag_module_domain_dict:
        each_locustag = each_module[:-4]
        domain_comb = locustag_module_domain_dict[each_module]

        each_module_substrates = {}
        domain_trunc_list = []

        for each_domain in domain_comb:
            abbr_domain = each_domain[:-5]
            domain_trunc_list.append(abbr_domain)

        discriminant = determine_module(domain_trunc_list)
        f_discriminant = determine_kr_activity(each_locustag, domain_comb, locustag_kr_dict, discriminant)

        if f_discriminant == 'None':
            continue

        elif f_discriminant == 'ACP':
            continue

        elif f_discriminant == 'PCP':
            continue

        else:
            module_currency_metab_dict = get_module_currency_metab_dict(f_discriminant, each_module, each_module_substrates, module_currency_metab_dict)

    #print "module_currency_metab_dict"
    #print module_currency_metab_dict, '\n'
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

    #print 'currency_metab_coeff_dict' 
    #print currency_metab_coeff_dict, '\n'
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

            #In case "consensus" is not reached:
            if locustag_monomer_dict[each_module][3] == 'nrp':
                #From NRPSPredictor2 SVM
                #Not considered: e.g., NRPSPredictor2 SVM: val,leu,ile,abu,iva
                #Checked by ',' in aSid_met2
                aSid_met2 = locustag_monomer_dict[each_module][0]
                if aSid_met2 != 'hydrophobic-aliphatic' and aSid_met2 != 'hydrophilic' and aSid_met2 != 'hydrophobic-aromatic' and aSid_met2 != 'N/A' and ',' not in aSid_met2:
                    biggid_met2 = get_biggid_from_aSid(aSid_met2)

                    #In case of non-consensus, NRPSPredictor2 SVM is considered 
                    metab_coeff_dict[biggid_met2] -= 1

                elif aSid_met2 == 'hydrophobic-aliphatic' or aSid_met2 == 'hydrophilic' or aSid_met2 == 'hydrophobic-aromatic' or aSid_met2 == 'N/A' or ',' in aSid_met2:
                    #If NRPSPredictor2 SVM has invalid monomer, then Minowa is considered
                    aSid_met4 = locustag_monomer_dict[each_module][2]
                    print "check: aSid_met4:", aSid_met4
                    if aSid_met4 != 'hydrophobic-aliphatic' and aSid_met4 != 'hydrophilic' and aSid_met4 != 'hydrophobic-aromatic' and aSid_met4 != 'N/A':
                        biggid_met4 = get_biggid_from_aSid(aSid_met4)
                        metab_coeff_dict[biggid_met4] -= 1

                    #If Minowa has invalid monomer, then Stachelhaus code is considered
                    elif aSid_met4 == 'hydrophobic-aliphatic' or aSid_met4 == 'hydrophilic' or aSid_met4 == 'hydrophobic-aromatic' or aSid_met4 == 'N/A':
                        aSid_met3 = locustag_monomer_dict[each_module][1]
                        biggid_met3 = get_biggid_from_aSid(aSid_met3)
                        metab_coeff_dict[biggid_met3] -= 1

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
                #In case of non-consensus, PKS signature is considered
                aSid_met6 = locustag_monomer_dict[each_module][0]
                print "check: aSid_met6", aSid_met6
                if aSid_met6 != 'N/A' and aSid_met6 != 'mal_or_prop':
                    biggid_met6 = get_biggid_from_aSid(aSid_met6)
                    #print "aSid_met6", aSid_met6, biggid_met6

                    metab_coeff_dict[biggid_met6] -= 1

                #If PKS signature has invalid monomer, then Minowa is considered
                else:
                    aSid_met7 = locustag_monomer_dict[each_module][1]
                    if aSid_met7 != 'inactive':
                        biggid_met7 = get_biggid_from_aSid(aSid_met7)
                        metab_coeff_dict[biggid_met7] -= 1

            #In case "consensus" is reached:
            elif locustag_monomer_dict[each_module][2] != 'pk':
                aSid_met8 = locustag_monomer_dict[each_module][2]

                #Original monomers are considered for those reduced by KR, DH and/or ER
                if aSid_met8 == 'ohmal' or aSid_met8 == 'ccmal' or aSid_met8 == 'redmal':
                    aSid_met8 = 'mal'
                elif aSid_met8 == 'ohmmal' or aSid_met8 == 'ccmmal' or aSid_met8 == 'redmmal':
                    aSid_met8 = 'mmal'
                elif aSid_met8 == 'ohmxmal' or aSid_met8 == 'ccmxmal' or aSid_met8 == 'redmxmal':
                    aSid_met8 = 'mxmal'
                elif aSid_met8 == 'ohemal' or aSid_met8 == 'ccemal' or aSid_met8 == 'redemal':
                    aSid_met8 = 'emal'

                biggid_met8 = get_biggid_from_aSid(aSid_met8)
                #print "aSid_met8", aSid_met8, biggid_met8
                metab_coeff_dict[biggid_met8] -= 1

    #Add secondary metabolite product to the reaction
    metab_coeff_dict[product] = 1

    #print 'metab_coeff_dict'
    #print metab_coeff_dict, '\n'
    return metab_coeff_dict


def add_sec_met_rxn(target_model, metab_coeff_dict, product, bigg_mnxm_compound_dict, mnxm_compoundInfo_dict, cluster_info_dict):
    
    #ID
    rxn = Reaction(product)

    #Reversibility / Lower and upper bounds
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Metabolites and their stoichiometric coeff's
    for metab in metab_coeff_dict.keys():

        #Consider only metabolites consumed or produced
        if metab_coeff_dict[metab] != 0:
            metab_compt = '_'.join([metab,'c'])

            #Adding metabolites already in the model
            if metab_compt in target_model.metabolites:
                rxn.add_metabolites({target_model.metabolites.get_by_id(metab_compt):metab_coeff_dict[metab]})

            #Adding metabolites with bigg compoundID, but not in the model
            elif metab in bigg_mnxm_compound_dict.keys():
                mnxm = bigg_mnxm_compound_dict[metab]
                metab_compt = Metabolite(metab_compt, formula = mnxm_compoundInfo_dict[mnxm][1], name = mnxm_compoundInfo_dict[mnxm][0], compartment='c')
                rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

            elif 'Cluster' in metab:
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

            #Add metabolite MNXM having no bigg ID to the model
            else:
                metab_compt = add_sec_met_mnxm_having_no_biggid_to_model(metab, metab_compt, mnxm_compoundInfo_dict)
                rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    #GPR association
    gpr_count = 0
    for each_gene in cluster_info_dict.keys():
        if gpr_count == 0:
            gpr_list = each_gene 
            gpr_count += 1
        else:
            gpr_list = gpr_list + ' AND ' + each_gene
     
    rxn.add_gene_reaction_rule(gpr_list)

    #Adding the new reaction to the model
    target_model.add_reaction(rxn)

    print "\n", "Secondary metabolite reaction:", rxn
    print rxn.reaction

    ##############################
    #Creating a transport reaction
    #Creating reaction ID
    rxn = Reaction("Transport_" + product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Adding a substrate metabolite
    #print cobra_model.metabolites.get_by_id(str(product_c))
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(product+'_c')):-1})

    #Adding product metabolite(s)
    product_e = Metabolite(product+"_e", name='', compartment='e')
    rxn.add_metabolites({product_e:1})

    #Adding the new reaction to the model
    target_model.add_reaction(rxn)

    print "\n", "Transport reaction:", rxn
    print rxn.reaction

    ##############################
    #Creating an exchange reaction
    #Creating reaction ID
    rxn = Reaction("Ex_"+product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible 0: irreversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Adding a substrate metabolite
    #print cobra_model.metabolites.get_by_id(str(product_c))
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(product_e)):-1})

    #Adding the new reaction to the model
    target_model.add_reaction(rxn)

    print "\n", "Exchange reaction:", rxn
    print rxn.reaction

    return target_model


def check_producibility_sec_met(target_model, product, dirname):
    for rxn in target_model.reactions:
        rxn.objective_coefficient = 0

    #target_model.reactions.get_by_id('Biomass_SCO').objective_coefficient = 0
    target_model.reactions.get_by_id("Ex_"+product).objective_coefficient = 1

    #Model reloading and overwrtting are necessary for model stability
    #Without these, model does not produce an accurate prediction
    write_cobra_model_to_sbml_file(target_model, dirname+'3_temp_models/'+"target_model_%s.xml" %product)
    target_model = create_cobra_model_from_sbml_file(dirname+'3_temp_models/'+"target_model_%s.xml" %product)

    target_model.optimize()
    print "Flux:", target_model.solution.f

    target_model.reactions.get_by_id('Biomass_SCO').objective_coefficient = 1 
    target_model.reactions.get_by_id("Ex_"+product).objective_coefficient = 0

    #return target_model.solution.f, product
    return target_model, product


def get_monomers_nonprod_sec_met(metab_coeff_dict):

    nonprod_sec_met_metab_list = []

    for metab in metab_coeff_dict.keys():
        #Exclude currency metabolites
        if metab_coeff_dict[metab] <0 and metab != 'atp' and metab != 'amp' and metab != 'ppi' and metab != 'amet' and metab != 'ahcys' and metab != 'fmn' and metab != 'fmnh2' and metab != 'nadp' and metab != 'nadph' and metab != 'h' and metab != 'h2o' and metab != 'hco3' and metab != 'coa':
            nonprod_sec_met_metab_list.append(metab)

    return nonprod_sec_met_metab_list


def get_monomers_prod_sec_met(metab_coeff_dict):

    prod_sec_met_metab_list = []

    for metab in metab_coeff_dict.keys():
        #Exclude currency metabolites
        if metab_coeff_dict[metab] <0 and metab != 'atp' and metab != 'amp' and metab != 'ppi' and metab != 'amet' and metab != 'ahcys' and metab != 'fmn' and metab != 'fmnh2' and metab != 'nadp' and metab != 'nadph' and metab != 'h' and metab != 'h2o' and metab != 'hco3' and metab != 'coa':
            prod_sec_met_metab_list.append(metab)

    return prod_sec_met_metab_list

