'''
Created by Kyu-Sang Hwang, Sept 2014

This file generates metabolic reactions for the genes newly annotated to be present in the secondary metabolite-biosynthetic gene cluster from antiSMASH.
'''

from Bio import SeqIO
from sets import Set
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file
from MNX_checker2 import COBRA_TemplateModel_checking_MNXref_metabolites
from MNX_checker2 import fix_legacy_id
import pickle
import copy

def monomers_list_for_nrps(inputFile):
    fp1 = open(inputFile,"r")
    monomer_mnx_dic = {}

    monomer = fp1.readline()

    while monomer:
        print monomer
        monomer = monomer.split("\t")
        monomer[0] = monomer[0].strip()
        monomer[1] = monomer[1].strip()
        monomer[2] = monomer[2].strip('\n')
        monomer_mnx_dic[monomer[0]] = [monomer[1],monomer[2]]
        monomer = fp1.readline()

    print "\n", "List of secondary metabolic monomers:"
    print monomer_mnx_dic
    fp1.close()
    return monomer_mnx_dic

def second_metab_monomers(gbkFile, FileType, monomer_mnx_dic):
#     fp1 = open('Output_second_metab_total_monomers.txt','w')

#Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)

    for feature in record.features:

#Identifies the feature "cluster"
        if feature.type == 'cluster':

#Identifies "Monomers prediction"
            qualifier_monomers = feature.qualifiers.get('note')
            qualifier_monomers = qualifier_monomers[2].split(':')
            qualifier_monomers = qualifier_monomers[1].strip()
            total_monomer_order = qualifier_monomers
            qualifier_monomers = qualifier_monomers.split(' + ')
            
            second_total_monomers = qualifier_monomers

#Modifies elements in list  
#second_total_monomers = ['pk-mmal-mal', 'mal-mal-mal-mmal', 'mmal-mal-mmal', 'mal-pk-mal']          
            count = 0
            
            for each_module_substrate in second_total_monomers:
                each_module_substrate = each_module_substrate.replace('(','')
                each_module_substrate = each_module_substrate.replace(')','')
                each_module_substrate = each_module_substrate.strip()
                second_total_monomers[count] = each_module_substrate
                
                count = count + 1

#Saves results as text file            
            count = 1
            
            for each_module_monomer in second_total_monomers:
#                 print >>fp1, "%s\t%s\t%s" % (total_monomer_order, count, each_module_monomer)
                
                count = count + 1
    
    print second_total_monomers
    
#     fp1.close()
    return second_total_monomers



def second_metab_reaction_product_names(gbkFile, FileType):
 
#Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)
    
    for feature in record.features:

#Retrieving "Cluster number"
        if feature.type == 'cluster':

            qualifier_cluster = feature.qualifiers.get('note')
            qualifier_cluster = qualifier_cluster[0].split(':')
            clusterNo = qualifier_cluster[1].strip()
            print "\n", "Cluster number:", clusterNo

#Retrieving "product"
            product = feature.qualifiers.get('product')
            product = product[0]
            print "Product:", product
            
            gene_strain = record.id
            gene_strain = gene_strain.split('.')
            gene_strain = gene_strain[0].strip()
            gene_strain = gene_strain.replace('_','') 
#             if gene_strain != None:
#                 gene_strain = gene_strain[0].split(':')
#                 gene_strain = gene_strain[1].split('.')
#                 gene_strain = gene_strain[0].strip()
#                 print "gene strain:", gene_strain
#             else:
#                 gene_strain = 'unknown'
#                 print "gene strain:", gene_strain
                  

    if float(clusterNo) < 10:
        product = gene_strain+"_"+"Cluster_0"+clusterNo+"_"+product
    else:
        product = gene_strain+"_"+"Cluster_"+clusterNo+"_"+product
    
    return product  

def second_metab_genes(gbkFile, FileType):
    fp1 = open('Output_second_metab_gene.txt','w')

    dic_pksnrps_gene = {}
    list_pksnrps_gene = []

#Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)

    for feature in record.features:
        
#Reads cluster_name from genebank file        
        if feature.type == 'cluster':
            gene_strain = record.id
            gene_strain = gene_strain.split('.')
            gene_strain = gene_strain[0].strip()
            gene_strain = gene_strain.replace('_','')
            print gene_strain
            cluster_info = feature.qualifiers.get('note')
            cluster_info = cluster_info[0].split(':')
            cluster_number = cluster_info[1].strip()
            
            substrate_order_info = feature.qualifiers.get('note')
            substrate_order_info = substrate_order_info[2].split(':')
            order_of_substrates = substrate_order_info[1].strip()
            print cluster_number, order_of_substrates
            
            if float(cluster_number) < 10:
                whole_cluster_name = gene_strain+".c00"+cluster_number
            else:
                whole_cluster_name = gene_strain+".c0"+cluster_number
            
            print whole_cluster_name
#             if cluster_info != None:
#                 spt_cluster_name = cluster_info[0].split(':')
#                 whole_cluster_name = spt_cluster_name[1]
#             else:
#                 gene_strain = record.id
#                 gene_strain = gene_strain.split('.')
#                 gene_strain = gene_strain[0].strip()
#                 gene_strain = gene_strain.replace('_','')

#                 sub_cluster_info = feature.qualifiers.get('note')
#                 sub_spt_cluster_name = sub_cluster_info[0].split(':')
#                 cluster_number = sub_spt_cluster_name[1].strip()
#                 whole_cluster_name = 'cluster'+'_'+cluster_number
            
#Feature CDS sometimes does not have "gene" in the qualifier.
        if feature.type == 'CDS':
            
            qualifier_product = feature.qualifiers.get('product')
            product = qualifier_product[0]
            
            qualifier_locus_tag = feature.qualifiers.get('locus_tag')
            locus_tag = qualifier_locus_tag[0]

#Feature CDS sometimes does not have "gene" in the qualifier. Therefore, this function is removed from this algorithm           
#             if 'aSProdPred' in feature.qualifiers:
#                 qualifier_aSProdPred = feature.qualifiers.get('aSProdPred')
#                 aSProdPred = qualifier_aSProdPred[0]

# Extract information of cluster_number, locus_tag, domain_inform, locus_tag, substrate_inform(aSProdPred) from query genebankfile            
            qualifier_sec_met = feature.qualifiers.get('sec_met')
            if qualifier_sec_met == None or len(qualifier_sec_met) < 4:
                continue
            
            qualifier_np_type = qualifier_sec_met[0]

            if qualifier_np_type[:5] == 'Type:':
                dic_pksnrps_gene[locus_tag] = [product, order_of_substrates]
                print >>fp1, "%s\t%s\t%s\t%s" % (whole_cluster_name, locus_tag, product, order_of_substrates)

    fp1.close()
    return dic_pksnrps_gene

def extracting_sub_set_met_info_from_genebank(gbkFile, FileType, dic_pksnrps_gene):

    dic_info_of_bundle_set_met = {}
     
#Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)
    
    for dic_gene_key in dic_pksnrps_gene:
        
        list_pksnrps_domain = []
        
        for feature in record.features: 
        
            if feature.type == 'CDS':
        
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')
                        
                if qualifier_locus_tag[0] == dic_gene_key:    

                    qualifier_set_met = feature.qualifiers.get('sec_met')
                    list_pksnrps_domain.append(qualifier_set_met)
        
        dic_info_of_bundle_set_met[dic_gene_key] = list_pksnrps_domain
                        
    return dic_info_of_bundle_set_met

def second_metab_domain(dic_info_of_bundle_set_met, dic_pksnrps_gene):
    fp1 = open('Output_second_metab_gene_domain.txt','w')
    fp2 = open('Output_second_metab_gene_substrate.txt','w')
    fp3 = open('Output_second_metab_gene_KR_activity.txt','w')
    
    dic_pksnrps_domain = {}
    dic_pksnrps_gene_domain = {}
    dic_pksnrps_gene_substrate = {}
    dic_t1pks_PKS_KR_activity = {}
        
    for each_gene in dic_pksnrps_gene:
        list_set_met = dic_info_of_bundle_set_met[each_gene][0]
        
        print each_gene
        print list_set_met

        domain_count = 0
        list_pksnrps_domain = []
        list_pksnrps_gene_substrate = []
        list_t1pks_gene_substrate = []
        list_t1pks_PKS_KR_activity = []
    
        for each_sub_set in list_set_met:

            if "NRPS/PKS Domain" in each_sub_set:
                            
                sptline1 = each_sub_set.split('. ')
                crude_domain_info = sptline1[0] 

                # to extract the information of KR activity and AT substrate specificity
                sptline2 = each_sub_set.split('; ')                
                
                if "AMP-binding" in crude_domain_info:
                    spt_NRPS_AMP = sptline2[1].split(':')
                    spt_NRPS_AMP_met = spt_NRPS_AMP[1].split(', ')
                    spt_NRPS_AMP_met = spt_NRPS_AMP_met[3].split(' (')
                    spt_NRPS_AMP_met = spt_NRPS_AMP_met[0].strip()
                    print spt_NRPS_AMP_met
                    
                    list_pksnrps_gene_substrate.append(spt_NRPS_AMP_met)
                    
                if "PKS_AT" in crude_domain_info:
                    spt_PKS_AT = sptline2[1].split(':')
                    spt_PKS_AT_met = spt_PKS_AT[1].split(', ')
                    spt_PKS_AT_met = spt_PKS_AT_met[2].split(' (')
                    spt_PKS_AT_met = spt_PKS_AT_met[0].strip()
                    print spt_PKS_AT_met
                    
                    list_pksnrps_gene_substrate.append(spt_PKS_AT_met)
                
                if "PKS_KR" in crude_domain_info:
                    spt_PKS_KR = sptline2[1].split(': ')
                    spt_PKS_KR = spt_PKS_KR[1].strip()
                     
                    list_t1pks_PKS_KR_activity.append(spt_PKS_KR)
                    
                spt_domain_info = crude_domain_info.split(':')
                whole_domain_info = spt_domain_info[1]

                spt_list_domain_info = whole_domain_info.split()
                spt_list_domain_info.append(each_gene)
                            
                each_nrps_domain = spt_list_domain_info[0]
                list_pksnrps_domain.append(each_nrps_domain)
                
                # domain information
                            
                domain_number = each_gene + '_DM' + str(domain_count)
                domain_count = domain_count + 1
                            
                dic_pksnrps_domain[domain_number] = spt_list_domain_info
                    
        dic_pksnrps_gene_domain[each_gene] = list_pksnrps_domain
        dic_pksnrps_gene_substrate[each_gene] = list_pksnrps_gene_substrate
        dic_t1pks_PKS_KR_activity[each_gene] = list_t1pks_PKS_KR_activity
        
        print each_gene, list_pksnrps_domain
        print each_gene, list_pksnrps_gene_substrate
        print each_gene, list_t1pks_PKS_KR_activity

        print >>fp1, "%s\t%s" % (each_gene, list_pksnrps_domain)
        print >>fp2, "%s\t%s" % (each_gene, list_pksnrps_gene_substrate)
        print >>fp3, "%s\t%s" % (each_gene, list_t1pks_PKS_KR_activity)

    fp1.close()
    fp2.close()
    fp3.close()      
    return dic_pksnrps_domain, dic_pksnrps_gene_domain, dic_pksnrps_gene_substrate, dic_t1pks_PKS_KR_activity

def second_metab_substrate(dic_info_of_bundle_set_met, dic_pksnrps_gene):
    fp1 = open('Output_second_metab_substrate_with_domain.txt','w')

    dic_pksnrps_domain_substrate = {}
    
    for each_gene in dic_pksnrps_gene:

        module_count = 0
        list_set_met =  dic_info_of_bundle_set_met[each_gene][0]
        
        for each_sub_set in list_set_met:
            
            discriminator = "true"
            print each_sub_set
     
            if "Substrate specificity predictions" in each_sub_set and "AMP-binding" in each_sub_set:
                
                list_participated_sustrate, discriminator = extract_substrate_information_nrps(each_sub_set, discriminator)
                
                module_number = each_gene + '_M' + str(module_count)
                                
                dic_pksnrps_domain_substrate[module_number] = list_participated_sustrate
                
                print module_number, list_participated_sustrate
                print >>fp1, "%s\t%s" % (module_number, list_participated_sustrate)
                        
                module_count = module_count + 1
                
            if "Substrate specificity predictions" in each_sub_set and "A-OX" in each_sub_set:
                
                list_participated_sustrate, discriminator = extract_substrate_information_nrps(each_sub_set, discriminator)
                
                module_number = each_gene + '_M' + str(module_count)
                                
                dic_pksnrps_domain_substrate[module_number] = list_participated_sustrate
                
                print module_number, list_participated_sustrate
                print >>fp1, "%s\t%s" % (module_number, list_participated_sustrate)

                module_count = module_count + 1
                
            if "Substrate specificity predictions" in each_sub_set and "PKS_AT" in each_sub_set:
                
                list_participated_sustrate = extract_substrate_information_pks(each_sub_set)
                
                module_number = each_gene + '_M' + str(module_count)
                                
                dic_pksnrps_domain_substrate[module_number] = list_participated_sustrate
                
                print module_number, list_participated_sustrate
                print >>fp1, "%s\t%s" % (module_number, list_participated_sustrate)
                
                module_count = module_count + 1
            
            if discriminator == "false":
                continue
    
    fp1.close()
    return dic_pksnrps_domain_substrate

def extract_substrate_information_nrps(each_sub_set, discriminator):
    
    sptline2 = each_sub_set.split(';')
    whole_substrate_info = sptline2[1]
    print each_sub_set
    print whole_substrate_info
    list_participated_sustrate = []
                                
    participated_substrates = whole_substrate_info.split(':')
    sptSubstrates = participated_substrates[1]
    print sptSubstrates
                                
    if ', ' not in sptSubstrates:
        print "insufficient substrate_info"
        discriminator = "false"
        return list_participated_sustrate
                
    substrates = sptSubstrates.split(', ')
    print substrates
                                
    for each_substrate in substrates:
                                               
        sptSubstrate = each_substrate.split('(')
        participated_substrate = sptSubstrate[0].strip()
        print participated_substrate
                                    
        list_participated_sustrate.append(participated_substrate)
    
    return list_participated_sustrate, discriminator

def extract_substrate_information_pks(each_sub_set):
    
    sptline2 = each_sub_set.split(';')
    whole_substrate_info = sptline2[1]
                                
    participated_substrates = whole_substrate_info.split(':')
    sptSubstrates = participated_substrates[1]
                                
    substrates = sptSubstrates.split(', ')
                                
    list_participated_sustrate = []
                                
    for each_substrate in substrates:
                                               
        sptSubstrate = each_substrate.split('(')
        participated_substrate = sptSubstrate[0].strip()
                                    
        list_participated_sustrate.append(participated_substrate)
    
    return list_participated_sustrate


def second_metab_module(dic_pksnrps_gene ,dic_pksnrps_gene_domain, dic_t1pks_PKS_KR_activity):
    fp1 = open('Output_second_metab_module.txt','w')
    
    dic_pksnrps_module = {}
    dic_pksnrps_module_KR_activity = {}
    
    for pksnrps_gene in dic_pksnrps_gene:
        
        count = 0
  
        if pksnrps_gene in dic_pksnrps_gene_domain:
            
            list_each_nrps_domain = dic_pksnrps_gene_domain[pksnrps_gene]
            print list_each_nrps_domain

            list_KR_activity = dic_t1pks_PKS_KR_activity[pksnrps_gene]
        
            list_module_info = []
            
            number_of_list = len(list_each_nrps_domain)

            KR_number = 0
            
            for each_domain in list_each_nrps_domain:    
                
#                 if each_domain == 'PKS_Docking_Nterm' or each_domain == 'PKS_Docking_Cterm':
#                     number_of_list = number_of_list - 1
#                     continue
                
                list_module_info.append(each_domain)
                number_of_list = number_of_list - 1
                
                if each_domain == 'PCP' or each_domain == 'ACP':
                    
                    module_number = pksnrps_gene + '_M' + str(count)
                    
                    # this algorithm should be pruned or make the new function.
                    if 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'active':
                        each_module_KR_activity = 1
                        KR_number = KR_number + 1
            
                    elif 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'inactive':
                        each_module_KR_activity = 2
                        KR_number = KR_number + 1
        
                    else:
                        each_module_KR_activity = 0
                                        
                    dic_pksnrps_module[module_number] = list_module_info
                    dic_pksnrps_module_KR_activity[module_number] = each_module_KR_activity
                    
                    print >>fp1, "%s\t%s\t%s\t%s" % (pksnrps_gene, module_number, list_module_info, each_module_KR_activity)
                    print module_number, list_module_info, each_module_KR_activity

                    list_module_info = []
                    
                    count = count + 1
                    
                elif each_domain == 'Epimerization':

                    count = count - 1
                    
                    module_number = pksnrps_gene + '_M' + str(count)
                    
                    # this algorithm should be pruned or make the new function.
                    if 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'active':
                        each_module_KR_activity = 1
                        KR_number = KR_number + 1
            
                    elif 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'inactive':
                        each_module_KR_activity = 2
                        KR_number = KR_number + 1
        
                    else:
                        each_module_KR_activity = 0
                                        
                    list_module_info = dic_pksnrps_module[module_number]

                    list_module_info.append('Epimerization')

                    A = dic_pksnrps_module.pop(module_number)
                    
                    dic_pksnrps_module[module_number] = list_module_info
                    dic_pksnrps_module_KR_activity[module_number] = each_module_KR_activity
                    
                    print >>fp1, "%s\t%s\t%s\t%s" % (pksnrps_gene, module_number, list_module_info, each_module_KR_activity)
                    print module_number, list_module_info, each_module_KR_activity
                    
                    list_module_info = []
                    
                    count = count + 1
                    
                elif each_domain == 'Thioesterase':
                    
                    count = count - 1
                    
                    module_number = pksnrps_gene + '_M' + str(count)
                    
                    # this algorithm should be pruned or make the new function.
                    if 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'active':
                        each_module_KR_activity = 1
                        KR_number = KR_number + 1
            
                    elif 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'inactive':
                        each_module_KR_activity = 2
                        KR_number = KR_number + 1
        
                    else:
                        each_module_KR_activity = 0
                                        
                    list_module_info = dic_pksnrps_module[module_number]
                    print list_module_info
                    
                    list_module_info.append('Thioesterase')
                    print list_module_info
                    
                    A = dic_pksnrps_module.pop(module_number)
                    
                    print A

                    dic_pksnrps_module[module_number] = list_module_info
                    dic_pksnrps_module_KR_activity[module_number] = each_module_KR_activity
                                        
                    print >>fp1, "%s\t%s\t%s\t%s" % (pksnrps_gene, module_number, list_module_info, each_module_KR_activity)
                    print module_number, list_module_info, each_module_KR_activity       
                
                elif list_module_info.count('PKS_KS') == 2:
                    
                    module_number = t1pks_gene + '_M' + str(count)
                    
                    poped_domain = list_module_info.pop()
                    
                    # this algorithm should be pruned or make the new function.
                    if 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'active':
                        each_module_KR_activity = 1
                        KR_number = KR_number + 1
            
                    elif 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'inactive':
                        each_module_KR_activity = 2
                        KR_number = KR_number + 1
        
                    else:
                        each_module_KR_activity = 0
                    #########
                    
                    dic_pksnrps_module[module_number] = list_module_info
                    dic_pksnrps_module_KR_activity[module_number] = each_module_KR_activity
                    
                    print >>fp1, "%s\t%s\t%s\t%s" % (pksnrps_gene, module_number, list_module_info, each_module_KR_activity)
                    print module_number, list_module_info, each_module_KR_activity

                    list_module_info = []

                    list_module_info.append(poped_domain)   
                    
                    count = count + 1              
                
                ## check required##  
                elif list_module_info.count('Condensation_DCL') == 2 or list_module_info.count('Condensation_LCL') == 2 or list_module_info.count('Condensation_LCL') + list_module_info.count('Condensation_DCL') == 2:
                    
                    module_number = pksnrps_gene + '_M' + str(count)
                    
                    list_module_info.pop()
                    
                    dic_pksnrps_module[module_number] = list_module_info

                    list_module_info = []
                    
                    count = count + 1
                    
                elif float(number_of_list) == 0:
                    
                    module_number = pksnrps_gene + '_M' + str(count)
                    
                    # this algorithm should be pruned or make the new function.
                    if 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'active':
                        each_module_KR_activity = 1
                        KR_number = KR_number + 1
            
                    elif 'PKS_KR' in list_module_info and list_KR_activity[KR_number] == 'inactive':
                        each_module_KR_activity = 2
                        KR_number = KR_number + 1
        
                    else:
                        each_module_KR_activity = 0
                    
                    dic_pksnrps_module[module_number] = list_module_info
                    dic_pksnrps_module_KR_activity[module_number] = each_module_KR_activity
                    
                    print >>fp1, "%s\t%s\t%s\t%s" % (pksnrps_gene, module_number, list_module_info, each_module_KR_activity)
                    print module_number, list_module_info, each_module_KR_activity

                    list_module_info = []
          
                    count = count + 1
                
    fp1.close()

    return dic_pksnrps_module, dic_pksnrps_module_KR_activity

##################################

def generating_each_module_of_backbone_biosynthesis_for_pksnrps(dic_pksnrps_module, dic_pksnrps_module_KR_activity):
    fp1 = open('Output_backbone_biosynthesis_of_each_module.txt','w')
    
    dic_converted_metabolic_reaction_without_substrate = {}

    for each_module in dic_pksnrps_module:
        
        domain_comb = dic_pksnrps_module[each_module]
        
        print domain_comb

        each_module_substrates = {}
        
        discriminant = module_discriminator(domain_comb)
        discriminant_with_KRact = Identifier_KR_activity(discriminant, each_module, dic_pksnrps_module_KR_activity)
        
        if discriminant == 'None':
            print "this discriminant is not defined : %s" % (domain_comb)
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, 'None')
            continue
        
        if discriminant == 'A' or discriminant == 'Aox': 
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates     
            print 'reaction 0: A or Aox', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
        
        elif discriminant == 'A_PCP' or discriminant == 'Cs_A_PCP' or discriminant == 'C_A_PCP' or discriminant == 'Cdcl_A_PCP' or discriminant == 'Clcl_A_PCP' or discriminant == 'Clcl_A_PCP' or discriminant == 'Cglyc_A_PCP' or discriminant == 'CXglyc_A_PCP': 
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates     
            print 'reaction 1: A_PCP or C_A_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'MT_A_PCP' or discriminant == 'Cs_A_MT_PCP' or discriminant == 'C_A_MT_PCP' or discriminant == 'Cdcl_A_MT_PCP' or discriminant == 'Clcl_A_MT_PCP' or discriminant == 'Cglyc_A_MT_PCP' or discriminant == 'CXglyc_A_MT_PCP': 
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates     
            print 'reaction 2: MT-A-PCP or C-A-MT-PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
        
        elif discriminant == 'Cs_A_PCP_E' or discriminant == 'C_A_PCP_E' or discriminant == 'Cdcl_A_PCP_E' or discriminant == 'Clcl_A_PCP_E' or discriminant == 'Cd_A_PCP' or discriminant == 'Cglyc_A_PCP_E' or discriminant == 'CXglyc_A_PCP_E':
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 3: C_A_PCP_E', each_module_substrates
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates) 
             
        elif discriminant == 'Cs_A_MT_E_PCP' or discriminant == 'C_A_MT_E_PCP' or discriminant == 'Cdcl_A_MT_E_PCP' or discriminant == 'Clcl_A_MT_E_PCP' or discriminant == 'Cglyc_A_MT_E_PCP' or discriminant == 'CXglyc_A_MT_E_PCP' or discriminant == 'Cd_A_MT_PCP':
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 4: C_A_MT_E_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
                  
        elif discriminant == 'HC_A_PCP' :
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['h2o'] = 2 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 5: HC_A_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'HC_Aox_PCP':
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            each_module_substrates['fmn'] = 1 #'Riboflavin-5-phosphate(FMN)', (C00061), 'MNXM119' 
            each_module_substrates['fmnh2'] = -1  # Reduced FMN (FMNH2)', (C01847), 'MNXM208'
            each_module_substrates['h2o'] = 2 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 6: HC_Aox_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'HC_A_MT_PCP' or discriminant == 'HC_A_MT_E_PCP':
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['ppi'] = 1 #'diphosphate', (C00013), 'MNXM11'
            each_module_substrates['h2o'] = 2 #H2O (C00001)
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 7: HC_A_MT_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'HC_Aox_MT_PCP' or discriminant == 'HC_Aox_MT_E_PCP':
            each_module_substrates['atp'] = -1 #'ATP' (C00002), 'MNXM3'
            each_module_substrates['amp'] = 1 #'AMP', (C00020), 'MNXM14'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            each_module_substrates['fmn'] = 1 #'Riboflavin-5-phosphate(FMN)', (C00061), 'MNXM119' 
            each_module_substrates['fmnh2'] = -1  # Reduced FMN (FMNH2)', (C01847), 'MNXM208'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            each_module_substrates['h2o'] = 2 #H2O (C00001)
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 8: HC_Aox_MT_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'AT_ACP' or discriminant == 'AT_KR(inactive)_ACP' or discriminant == 'AT_DH_KR(inactive)_ACP' or discriminant == 'AT_DH_ER_KR(inactive)_ACP': 
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates     
            print 'reaction 1: AT-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'AT_KS_ACP' or discriminant == 'AT_KS' or discriminant == 'AT_KS_KR(inactive)_ACP' or discriminant == 'AT_KS_KR(inactive)' or discriminant == 'AT_KS_DH_KR(inactive)_ACP' or discriminant == 'AT_KS_DH_KR(inactive)' or discriminant == 'AT_KS_DH_KR(inactive)_ACP' or discriminant == 'AT_KS_DH_ER_KR(inactive)': 
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates     
            print 'reaction 2: AT-KS-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
             
        elif discriminant == 'AT_KS_KR_ACP' or discriminant == 'AT_KR_ACP' or discriminant == 'AT_KS_KR' or discriminant == 'AT_ER_KR_ACP' or discriminant == 'AT_KS_ER_KR' or discriminant == 'AT_KS_ER_KR_ACP':
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 3: KS-AT-KR-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
             
        elif discriminant == 'AT_KS_DH_KR_ACP' or discriminant == 'AT_DH_KR_ACP' or discriminant == 'AT_KS_DH_KR':
            each_module_substrates['coa'] = 1 #'Coenzyme A', (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 4: KS-AT-DH-KR-ACP', each_module_substrates
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates) 
             
        elif discriminant == 'AT_KS_DH_ER_KR_ACP' or discriminant == 'AT_DH_ER_KR_ACP' or discriminant == 'AT_KS_DH_ER_KR':
            each_module_substrates['coa'] = 1 #'Coenzyme A', (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['h2o'] = 1 #'H2O', (C00001), 'MNXM2'
            each_module_substrates['nadp'] = 2 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -2 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -2 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 5: KS-AT-DH-ER-KR-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'AT_KS_cMT_ACP' or discriminant == 'AT_KS_cMT':
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 6: AT_KS_cMT_ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'AT_KS_KR_cMT_ACP' or discriminant == 'AT_KS_KR_cMT':
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 7: AT_KS_KR_cMT_ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'AT_KS_DH_KR_cMT_ACP' or discriminant == 'AT_KS_DH_KR_cMT':
            each_module_substrates['coa'] = 1 #'Coenzyme A', (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 8: AT_KS_DH_KR_cMT_ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'AT_KS_DH_KR_cMT_ER_ACP' or discriminant == 'AT_KS_DH_KR_cMT_ER':
            each_module_substrates['coa'] = 1 #'Coenzyme A', (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['h2o'] = 1 #'H2O', (C00001), 'MNXM2'
            each_module_substrates['nadp'] = 2 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -2 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -2 #'H+', (C00080), 'MNXM1'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 8: AT_KS_DH_KR_cMT_ER_ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'TD':
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 8: HC_Aox_MT_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'ACP':
            continue
            
        elif discriminant == 'PCP':
            continue
    
    fp1.close()
    
    print dic_converted_metabolic_reaction_without_substrate
                
    return dic_converted_metabolic_reaction_without_substrate

def module_discriminator(domain_comb):
## domain information(nrps) :
## Condensation    Condensation domain
## Condensation_DCL    Condensation domain that links L-amino acid to peptide ending with D-amino acid
## Condensation_LCL    Condensation domain that links L-amino acid to peptide ending with L-amino acid 
## Condensation_Dual     Dual condensation / epimerization domain 
## Condensation_Starter    Starter condensation domain
## CXglyc    Putatively inactive glycopeptide condensation-like domain
## Cglyc    Glycopeptide condensation domain 
## Heterocyclization    Heterocyclization domain
## Epimerization    Epimerization domain 
## AMP-binding    Adenylation domain
## A-OX     Adenylation domain with integrated oxidase
## PCP     Peptidyl-carrier protein domain 
## ACPS    4'-phosphopantetheinyl transferase
## NRPS-COM_Nterm    NRPS COM domain Nterminal
## NRPS-COM_Cterm    NRPS COM domain Cterminal

## domain information(pks) :
## AT    acyltransferase
## KS    ketosynthase
## ACP    acyl carrier protein
## KR    ketoreductase
## DH    dehydratase
## ER    enolase
## cMT    methyltransferase
## TD    thiolesterase domain

    print domain_comb

# Exceptionsal cases
    if 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'A'
    
# Starter unit_AMP-binding+PCP
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'A_PCP' #
        
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'A_PCP'
        
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'MT_A_PCP' ##
        
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'MT_A_PCP'
        
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'MT_A_PCP'
        
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'MT_A_PCP'
              
#####
# Condensation starter domain : C-domain
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_PCP' #
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_PCP' #    
       
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##    
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_PCP_E' ####
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_PCP_E' ####
         
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP' ###
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP' ###

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP'
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP'
        
# Condensation domain that links D-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_PCP' #
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_PCP' #
         
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##
    
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##
          
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_PCP_E' ####
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_PCP_E' ####
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###     
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###     

# Condensation domain that links D-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_PCP' #
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_PCP' #
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_PCP_E' ####
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_PCP_E' ####
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###

# Condensation domain that links L-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_PCP' #
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_PCP' #
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_PCP_E' ####
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_PCP_E' ####
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###

# Condensation_dual : C-domain
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cd_A_PCP' #### + epimerization
        
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cd_A_PCP' #### + epimerization
        
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_com:
        discriminant = 'Cd_A_MT_PCP' ###
        
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cd_A_MT_PCP' ###

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cd_A_MT_PCP' ###
        
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cd_A_MT_PCP' ###
        
# Glycopeptide condensation domain (O) : C-domain
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_PCP' #
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_PCP' #
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_PCP_E' ####
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_PCP_E' ####
          
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###


# Glycopeptide condensation domain (X) : C-domain
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_PCP' #
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_PCP' #
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_PCP_E' ####
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_PCP_E' ####
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###

        
# Heterocyclization : C-domain
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_PCP' ###
          
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'  
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'  
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'  

# t1pks modules
# Starter units 
    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KR_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_DH_KR_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_DH_ER_KR_ACP'
        
# Exceptional case
    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_ER_KR_ACP'
        
    elif 'PKS_AT' not in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'ACP'

# Extension unit    
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS'
    
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_KR_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_KR'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_KR_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_KR'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_ER_KR_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_ER_KR'
        
# exeptional case
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_ER_KR_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_ER_KR'
    
# Methytransferase
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_cMT_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_cMT'
    
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_KR_cMT_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_KR_cMT'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT_ACP'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT'
        
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT_ER_ACP'
    
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT_ER'


# Terminal reductase domain : terminal domain       
    elif 'TD' in domain_comb and 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'TD'
    
    else:
        discriminant = 'None'
        
    return discriminant

def Identifier_KR_activity(discriminant, each_module, dic_pksnrps_module_KR_activity):
     
    print discriminant 
    print each_module
    print dic_pksnrps_module_KR_activity
    
    each_module_KR_activity = dic_pksnrps_module_KR_activity[each_module]
    print each_module_KR_activity
 
    if each_module_KR_activity == 1 or each_module_KR_activity == 0:
        discriminant = discriminant
        return discriminant
     
    elif each_module_KR_activity == 2:
        discriminant = discriminant.replace('KR','KR(inactive)')      
        return discriminant

def integrated_metabolic_reaction1(participated_cofactor_info):
    fp1 = open('Output_participated_cofactors.txt','w')
    
    dic_integrated_metabolic_reaction = {}
    dic_integrated_metabolic_reaction['atp'] = 0 #'ATP' (C00002), 'MNXM3'
    dic_integrated_metabolic_reaction['amp'] = 0 #'AMP', (C00020), 'MNXM14'
    dic_integrated_metabolic_reaction['ppi'] = 0 #'diphosphate', (C00013), 'MNXM11'
    dic_integrated_metabolic_reaction['amet'] = 0 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
    dic_integrated_metabolic_reaction['ahcys'] = 0 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
    dic_integrated_metabolic_reaction['fmn'] = 0 #'Riboflavin-5-phosphate(FMN)', (C00061), 'MNXM119' 
    dic_integrated_metabolic_reaction['fmnh2'] = 0  # Reduced FMN (FMNH2)', (C01847), 'MNXM208'
    dic_integrated_metabolic_reaction['amet'] = 0 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
    dic_integrated_metabolic_reaction['ahcys'] = 0 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
    dic_integrated_metabolic_reaction['nadp'] = 0 #'NADP+', (C00006), 'MNXM5'
    dic_integrated_metabolic_reaction['nadph'] = 0 #'NADPH' ,(C00005), 'MNXM6'
    dic_integrated_metabolic_reaction['h'] = 0 #'H+', (C00080), 'MNXM1'
    dic_integrated_metabolic_reaction['h2o'] = 0 #'H2O', (C00001), 'MNXM2'
    dic_integrated_metabolic_reaction['hco3'] = 0
    dic_integrated_metabolic_reaction['coa'] = 0
    
    for each_module in participated_cofactor_info:
        
        for each_metabolite in participated_cofactor_info[each_module]:
            met_coff = participated_cofactor_info[each_module][each_metabolite]
                   
            if participated_cofactor_info[each_module][each_metabolite] > 0:
                dic_integrated_metabolic_reaction[each_metabolite] += met_coff
                print >>fp1, "participated_cofactors:\t%s\t%s\t%s" % (each_module, each_metabolite, met_coff)
                
            else:
                dic_integrated_metabolic_reaction[each_metabolite] += met_coff
                print >>fp1, "participated_cofactors:\t%s\t%s\t%s" % (each_module, each_metabolite, met_coff)
    

    print >>fp1, "#####"
    for each_participated_cofactor in dic_integrated_metabolic_reaction:
        print each_participated_cofactor, dic_integrated_metabolic_reaction[each_participated_cofactor]
        print >>fp1, "total_participated_cofactor:\t%s\t%s" % (each_participated_cofactor, dic_integrated_metabolic_reaction[each_participated_cofactor])

    fp1.close()
    
    return dic_integrated_metabolic_reaction

def integrated_metabolic_reaction2(dic_pksnrps_domain_substrate, dic_integrated_metabolic_reaction_without_cofactors):
    fp1 = open('Output_participated_substrates_from_uniformed_prediction.txt','w')
        
    list_of_dismatched_substrate = []
    
    dic_integrated_metabolic_reaction_without_cofactors['ala_DASH_L'] = 0 #'L-alanine', 'C00041', 'MNXM32'
    dic_integrated_metabolic_reaction_without_cofactors['arg_DASH_L'] = 0 #'L-Arginine', 'C00062', 'MNXM70'
    dic_integrated_metabolic_reaction_without_cofactors['asn_DASH_L'] = 0 #'L-asparagine', 'C00152', 'MNXM147'
    dic_integrated_metabolic_reaction_without_cofactors['asp_DASH_L'] = 0 #'L-Aspartate', 'C00049', 'MNXM42'
    dic_integrated_metabolic_reaction_without_cofactors['cys_DASH_L'] = 0  #'L-Cysteine', 'C00097', 'MNXM55'
    dic_integrated_metabolic_reaction_without_cofactors['gln_DASH_L'] = 0 #'L-Glutamine', 'C00064', 'MNXM37'
    dic_integrated_metabolic_reaction_without_cofactors['glu_DASH_L'] = 0 #'L-Glutamate', 'C00025', 'MNXM89557'
    dic_integrated_metabolic_reaction_without_cofactors['gly'] = 0 #'glycine', 'C00037', 'MNXM29'
    dic_integrated_metabolic_reaction_without_cofactors['his_DASH_L'] = 0 #'L-Histidine', 'C00135', 'MNXM134'
    dic_integrated_metabolic_reaction_without_cofactors['leu_DASH_L'] = 0 #'L-Leucine', 'C00123', 'MNXM140'
    dic_integrated_metabolic_reaction_without_cofactors['lys_DASH_L'] = 0 #'L-Lysine', 'C00047 ', 'MNXM78'
    dic_integrated_metabolic_reaction_without_cofactors['met_DASH_L'] = 0 #'L-Methionine', 'C00073', 'MNXM61'
    dic_integrated_metabolic_reaction_without_cofactors['phe_DASH_L'] = 0 #'L-Phenylalanine', 'C00079', 'MNXM97'
    dic_integrated_metabolic_reaction_without_cofactors['pro_DASH_L'] = 0 #'L-Proline', 'C00148', 'MNXM114'
    dic_integrated_metabolic_reaction_without_cofactors['ser_DASH_L'] = 0 #'L-Serine', 'C00065', 'MNXM53'
    dic_integrated_metabolic_reaction_without_cofactors['thr_DASH_L'] = 0 #'L-Threonine', 'C00188', 'MNXM142'
    dic_integrated_metabolic_reaction_without_cofactors['trp_DASH_L'] = 0 #'L-Tryptophan', 'C00078', 'MNXM94'
    dic_integrated_metabolic_reaction_without_cofactors['tyr_DASH_L'] = 0 #'L-Tyrosine', 'C00082', 'MNXM76'
    dic_integrated_metabolic_reaction_without_cofactors['val_DASH_L'] = 0 #'L-Valine', 'C00183', 'MNXM199'
    dic_integrated_metabolic_reaction_without_cofactors['ile_DASH_L'] = 0 #'L-Isoleucine', 'C00407', 'MNXM231'
    dic_integrated_metabolic_reaction_without_cofactors['phg_DASH_L'] = 0 #'phenylglycine', 'C18623', 'MNXM59292'
    dic_integrated_metabolic_reaction_without_cofactors['bht_DASH_L'] = 0 #'beta-hydroxyn-tyrosine', 'N/A', 'N/A'
    dic_integrated_metabolic_reaction_without_cofactors['orn'] = 0 #'Ornithine', 'C01602', 'MNXM89689'
    dic_integrated_metabolic_reaction_without_cofactors['abu'] = 0 #'D-2-Aminobutyric acid', 'C02261', 'MNXM17054'
    dic_integrated_metabolic_reaction_without_cofactors['iva'] = 0 #'2-Amino-2-methylbutanoate', 'C03571', 'MNXM34821'
    dic_integrated_metabolic_reaction_without_cofactors['L2aadp'] = 0 #'L-2-Aminoadipic acid', 'C00956', 'MNXM268'
    dic_integrated_metabolic_reaction_without_cofactors['hpg'] = 0 #'D-4-Hydroxyphenylglycine', 'C03493', 'N/A'
    dic_integrated_metabolic_reaction_without_cofactors['23dhb'] = 0 #'2,3-Dihydroxybenzoic acid', 'C00196', 'MNXM455'
    dic_integrated_metabolic_reaction_without_cofactors['dhpg'] = 0 #'3,5-Dihydroxy-phenylglycine', 'C12026', 'MNXM9962'
    dic_integrated_metabolic_reaction_without_cofactors['hty'] = 0 #'L-Homotyrosine', 'C18622', 'MNXM59438'
    dic_integrated_metabolic_reaction_without_cofactors['citr_DASH_L'] = 0 #'L-citruline', 'C00327', 'MNXM211'
    dic_integrated_metabolic_reaction_without_cofactors['Lpipecol'] = 0 #'L-pipecolate', 'C00408', 'MNXM684'
    dic_integrated_metabolic_reaction_without_cofactors['ala_DASH_B'] = 0 #'beta-alanine zwitterion', 'C00099', 'MNXM144'
    dic_integrated_metabolic_reaction_without_cofactors['24dab'] = 0 #'L-2,4-diazaniumylbutyrate', 'C03283', 'MNXM840'
    dic_integrated_metabolic_reaction_without_cofactors['pac'] = 0 #'phenylacetate', 'C00548', 'MNXM497'
    dic_integrated_metabolic_reaction_without_cofactors['tcl'] = 0 #'triclosan', 'C12059', 'MNXM87879'
    dic_integrated_metabolic_reaction_without_cofactors['qa'] = 0 #'quinoxaline', 'C18575','MNXM80501'
    
    dic_integrated_metabolic_reaction_without_cofactors['malcoa'] = 0 #'malonyl-CoA', 'C00083', 'MNXM40'
    dic_integrated_metabolic_reaction_without_cofactors['mmcoa_DASH_S'] = 0 #'(S)-methylmalonyl-CoA(5-)','C00683', 'MNXM190', 'not detected in bigg database'
    dic_integrated_metabolic_reaction_without_cofactors['2mbcoa'] = 0 #'2-methylbutanoyl-CoA', C01033,'MNXM569'
    dic_integrated_metabolic_reaction_without_cofactors['emcoa_DASH_S'] = 0 #'ethylmalonyl-CoA','C18026', 'MNXM2043', 'not detected in bigg database'
    dic_integrated_metabolic_reaction_without_cofactors['ibcoa'] = 0 #'2-Methylpropanoyl-CoA', 'C00630', 'MNXM470'
    dic_integrated_metabolic_reaction_without_cofactors['accoa'] = 0 #'Acetyl-CoA', 'C00024', 'MNXM21'
    dic_integrated_metabolic_reaction_without_cofactors['ppcoa'] = 0 #'Propionyl-CoA', 'C00100', 'MNXM86'
    dic_integrated_metabolic_reaction_without_cofactors['ivcoa'] = 0 #'3-Methylbutanoyl-CoA', 'C02939', 'MNXM471'
    dic_integrated_metabolic_reaction_without_cofactors['mxmalacp'] = 0 #'Methoxymalonyl-[acp]A', 'C18616', 'MNXM61686'
    dic_integrated_metabolic_reaction_without_cofactors['chccoa'] = 0 #'cyclohexane-1-carboxyl-CoA', 'C09823', 'MNXM5111'
    
    for each_module in dic_pksnrps_domain_substrate:

        print dic_pksnrps_domain_substrate
        print len(dic_pksnrps_domain_substrate[each_module])

        if len(dic_pksnrps_domain_substrate[each_module]) == 4:
            
            sptlist1 = dic_pksnrps_domain_substrate[each_module][0].split(',')
            temp_list = []
            
            if dic_pksnrps_domain_substrate[each_module][3] == 'nrp':
            
                if len(sptlist1) >= 2:
                    for each_substrate in sptlist1:
                        converted_met1 = converting_substrate_to_met_coeff(each_substrate)
                        temp_list.append(converted_met1)
                elif dic_pksnrps_domain_substrate[each_module][0] == 'hydrophobic-aliphatic' or dic_pksnrps_domain_substrate[each_module][0] == 'hydrophilic':
                    each_substrate = 'N/A'
                    temp_list.append(each_substrate)    
                else:
                    query_met2 = dic_pksnrps_domain_substrate[each_module][0]
                    converted_met2 = converting_substrate_to_met_coeff(query_met2)
                    temp_list.append(converted_met2)
                
                query_met3 = dic_pksnrps_domain_substrate[each_module][1]
                converted_met3 = converting_substrate_to_met_coeff(query_met3)
                temp_list.append(converted_met3)
            
                query_met4 = dic_pksnrps_domain_substrate[each_module][2]
                converted_met4 = converting_substrate_to_met_coeff(query_met4)
                temp_list.append(converted_met4)
            
                list_of_dismatched_substrate.append(temp_list)           
        
            elif dic_pksnrps_domain_substrate[each_module][3] != 'nrp':
                query_met5 = dic_pksnrps_domain_substrate[each_module][3]
                converted_met5 = converting_substrate_to_met_coeff(query_met5)
                dic_integrated_metabolic_reaction_without_cofactors[converted_met5] -= 1
            
        elif len(dic_pksnrps_domain_substrate[each_module]) == 3:
            
            if len(dic_pksnrps_domain_substrate[each_module]) < 3:
                continue
            # seperating between unified prediction and not
            
            if dic_pksnrps_domain_substrate[each_module][2] == 'pk':
                temp_list2 = []
                
                query_met6 = dic_pksnrps_domain_substrate[each_module][0]
                converted_met6 = converting_substrate_to_met_coeff(query_met6)
                temp_list2.append(converted_met6)
                
                query_met7 = dic_pksnrps_domain_substrate[each_module][1]
                converted_met7 = converting_substrate_to_met_coeff(query_met7)
                temp_list2.append(converted_met7)
                
                list_of_dismatched_substrate.append(temp_list2)
        
            elif dic_pksnrps_domain_substrate[each_module][2] != 'pk':
                query_met8 = dic_pksnrps_domain_substrate[each_module][2]
                converted_met8 = converting_substrate_to_met_coeff(query_met8)
                dic_integrated_metabolic_reaction_without_cofactors[converted_met8] -= 1
             
    for each_participated_substrate in dic_integrated_metabolic_reaction_without_cofactors:
        print >>fp1, "total_participated_cofactor_with_substrates1:\t%s\t%s" % (each_participated_substrate, dic_integrated_metabolic_reaction_without_cofactors[each_participated_substrate])

    print dic_integrated_metabolic_reaction_without_cofactors
    print list_of_dismatched_substrate
        
    fp1.close()
    
    return dic_integrated_metabolic_reaction_without_cofactors, list_of_dismatched_substrate


def converting_substrate_to_met_coeff(each_substrate):

# nrps substrates 
    if each_substrate == 'ala':
        met_name = 'ala_DASH_L'
        
    elif each_substrate == 'arg':
        met_name = 'arg_DASH_L' 
        
    elif each_substrate == 'asn':
        met_name ='asn_DASH_L'
            
    elif each_substrate == 'asp':
        met_name = 'asp_DASH_L'
            
    elif each_substrate == 'cys':
        met_name = 'cys_DASH_L'
            
    elif each_substrate == 'gln':
        met_name = 'gln_DASH_L'
            
    elif each_substrate == 'glu':
        met_name = 'glu_DASH_L'
            
    elif each_substrate == 'gly':
        met_name = 'gly'
            
    elif each_substrate == 'his':
        met_name = 'his_DASH_L'
            
    elif each_substrate == 'leu':
        met_name = 'leu_DASH_L'
            
    elif each_substrate == 'lys':
        met_name = 'lys_DASH_L'
            
    elif each_substrate == 'met':
        met_name = 'met_DASH_L'
            
    elif each_substrate == 'phe':
        met_name = 'phe_DASH_L'
            
    elif each_substrate == 'pro':
        met_name = 'pro_DASH_L'
            
    elif each_substrate == 'ser':
        met_name = 'ser_DASH_L'
            
    elif each_substrate == 'thr':
        met_name = 'thr_DASH_L'
             
    elif each_substrate == 'trp':
        met_name = 'trp_DASH_L'
             
    elif each_substrate == 'tyr':
        met_name ='tyr_DASH_L'
             
    elif each_substrate == 'val':
        met_name = 'val_DASH_L'
        
    elif each_substrate == 'ile':
        met_name = 'ile_DASH_L'
        
    elif each_substrate == 'phg':
        met_name = 'phg_DASH_L'
        
    elif each_substrate == 'bht':
        met_name = 'bht_DASH_L'
        
    elif each_substrate == 'orn':
        met_name = 'orn'
        
    elif each_substrate == 'abu':
        met_name = 'abu'
        
    elif each_substrate == 'iva':
        met_name = 'iva'
        
    elif each_substrate == 'aad':
        met_name = 'L2aadp'
        
    elif each_substrate == 'hpg':
        met_name = 'hpg'
        
    elif each_substrate == 'dhb':
        met_name = '23dhb'
        
    elif each_substrate == 'dhpg':
        met_name = 'dhpg'
        
    elif each_substrate == 'hty':
        met_name = 'hty'
        
    elif each_substrate == 'cit':
        met_name = 'citr_DASH_L'
        
    elif each_substrate == 'pip':
        met_name = 'Lpipecol'
        
    elif each_substrate == 'b-ala':
        met_name = 'ala_DASH_B'
        
    elif each_substrate == 'dab':
        met_name = '24dab'
        
    elif each_substrate == 'phenylacetate' or each_substrate == 'Pha':
        met_name = 'pac'
        
    elif each_substrate == 'tcl':
        met_name = 'tcl'
        
    elif each_substrate == 'qa':
        met_name = 'qa'
        
# t1pks substreate
    elif each_substrate == 'mal':
        met_name = 'malcoa'
        
    elif each_substrate == 'mmal':
        met_name = 'mmcoa_DASH_S'

    elif each_substrate == '2metbut':
        met_name = '2mbcoa'
        
    elif each_substrate == 'Ethyl_mal' or each_substrate == 'emal':
        met_name = 'emcoa_DASH_S'
            
    elif each_substrate == 'isobut':
        met_name = 'ibcoa'
            
    elif each_substrate == 'ace':
        met_name = 'accoa'
            
    elif each_substrate == 'prop':
        met_name = 'ppcoa'
            
    elif each_substrate == '3metbut':
        met_name = 'ivcoa'

    elif each_substrate == 'mxmal':
        met_name = 'mxmalacp'
        
    elif each_substrate == 'CHC-CoA':
        met_name = 'chccoa'

    elif each_substrate == 'N/A':
        met_name = 'N/A'        
        
    else:
        print each_substrate
        raw_input('substrate_not_defined')
    
    return met_name

def integrated_metabolic_reaction3(dic_integrated_metabolic_reaction, list_of_dismatched_substrate):
    
    list_of_reaction_set = []

    if list_of_dismatched_substrate == []:
        
        list_of_reaction_set.append(dic_integrated_metabolic_reaction)
        
    else:
        
        list_of_reaction_set.append(dic_integrated_metabolic_reaction)
        
        for each_pair_of_substrates in list_of_dismatched_substrate:
            print each_pair_of_substrates
    
            template_list = copy.deepcopy(list_of_reaction_set)
            
            list_of_reaction_set = []
        
            for dic_each_metabolic_reaction in template_list:
                
                each_pair_of_substrates = set(each_pair_of_substrates)
                each_pair_of_substrates = list(each_pair_of_substrates)
                
                print each_pair_of_substrates
                
                substrate_decision_number = distincting_each_substrate_in_list_component(each_pair_of_substrates)
                
                temp_reaction_set = converting_nrps_substrates(each_pair_of_substrates, dic_each_metabolic_reaction, substrate_decision_number) 
                
                for each_dic_reaction_set in temp_reaction_set:
                    
                    list_of_reaction_set.append(each_dic_reaction_set)

    return list_of_reaction_set
                
def distincting_each_substrate_in_list_component(each_pair_of_substrates):    
    
    # this logic of code should be fixed
    if len(each_pair_of_substrates) >= 2:
        substrate_decision_number = 1 
          
    else:
        substrate_decision_number = 0
        
    return substrate_decision_number

def converting_nrps_substrates(each_pair_of_substrates, dic_integrated_metabolic_reaction, substrate_decision_number):
    
    temp_list_of_reaction_set = []
    
    for each_substrate in each_pair_of_substrates:
        
        temp_dic ={}
        temp_dic = copy.deepcopy(dic_integrated_metabolic_reaction)
        
        if each_substrate not in temp_dic:
            continue
        
        if each_substrate == 'ala_DASH_L':
            temp_dic['ala_DASH_L'] -= 1
        
        elif each_substrate == 'arg_DASH_L':
            temp_dic['arg_DASH_L'] -= 1 
        
        elif each_substrate == 'asn_DASH_L':
            temp_dic['asn_DASH_L'] -= 1
            
        elif each_substrate == 'asp_DASH_L':
            temp_dic['asp_DASH_L'] -= 1
            
        elif each_substrate == 'cys_DASH_L':
            temp_dic['cys_DASH_L'] -= 1
            
        elif each_substrate == 'gln_DASH_L':
            temp_dic['gln_DASH_L'] -= 1
            
        elif each_substrate == 'glu_DASH_L':
            temp_dic['glu_DASH_L'] -= 1
            
        elif each_substrate == 'gly':
            temp_dic['gly'] -= 1
            
        elif each_substrate == 'his_DASH_L':
            temp_dic['his_DASH_L'] -= 1
            
        elif each_substrate == 'leu_DASH_L':
            temp_dic['leu_DASH_L'] -= 1
            
        elif each_substrate == 'lys_DASH_L':
            temp_dic['lys_DASH_L'] -= 1
            
        elif each_substrate == 'met_DASH_L':
            temp_dic['met_DASH_L'] -= 1
            
        elif each_substrate == 'phe_DASH_L':
            temp_dic['phe_DASH_L'] -= 1
            
        elif each_substrate == 'pro_DASH_L':
            temp_dic['pro_DASH_L'] -= 1
            
        elif each_substrate == 'ser_DASH_L':
            temp_dic['ser_DASH_L'] -= 1
            
        elif each_substrate == 'thr_DASH_L':
            temp_dic['thr_DASH_L'] -= 1
            
        elif each_substrate == 'ile_DASH_L':
            temp_dic['ile_DASH_L'] -= 1
             
        elif each_substrate == 'trp_DASH_L':
            temp_dic['trp_DASH_L'] -= 1
             
        elif each_substrate == 'tyr_DASH_L':
            temp_dic['tyr_DASH_L'] -= 1
             
        elif each_substrate == 'val_DASH_L':
            temp_dic['val_DASH_L'] -= 1
            
        elif each_substrate == 'phg_DASH_L':
            temp_dic['phg_DASH_L'] -= 1
        
        elif each_substrate == 'bht_DASH_L':
            temp_dic['bht_DASH_L'] -= 1
            
        elif each_substrate == 'orn':
            temp_dic['orn'] -= 1
            
        elif each_substrate == 'abu':
            temp_dic['abu'] -= 1
            
        elif each_substrate == 'iva':
            temp_dic['iva'] -= 1
            
        elif each_substrate == 'L2aadp':
            temp_dic['L2aadp'] -= 1
            
        elif each_substrate == 'hpg':
            temp_dic['hpg'] -= 1
            
        elif each_substrate == '23dhb':
            temp_dic['23dhb'] -= 1
            
        elif each_substrate == 'dhpg':
            temp_dic['dhpg'] -= 1
            
        elif each_substrate == 'hty':
            temp_dic['hty'] -= 1
            
        elif each_substrate == 'citr_DASH_L':
            temp_dic['citr_DASH_L'] -= 1
            
        elif each_substrate == 'Lpipecol':
            temp_dic['Lpipecol'] -= 1
            
        elif each_substrate == 'ala_DASH_B':
            temp_dic['ala_DASH_B'] -= 1
            
        elif each_substrate == '24dab':
            temp_dic['24dab'] -= 1
            
        elif each_substrate == 'pac':
            temp_dic['pac'] -= 1
            
        elif each_substrate == 'tcl':
            temp_dic['tcl'] -= 1
            
        elif each_substrate == 'qa':
            temp_dic['qa'] -= 1
            
        elif each_substrate == 'malcoa':
            temp_dic['malcoa'] -= 1
        
        elif each_substrate == 'mmcoa_DASH_S':
            temp_dic['mmcoa_DASH_S'] -= 1

        elif each_substrate == '2mbcoa':
            temp_dic['2mbcoa'] -= 1
        
        elif each_substrate == 'emcoa_DASH_S':
            temp_dic['emcoa_DASH_S'] -= 1
            
        elif each_substrate == 'ibcoa':
            temp_dic['ibcoa'] -= 1
            
        elif each_substrate == 'accoa':
            temp_dic['accoa'] -= 1
            
        elif each_substrate == 'ppcoa':
            temp_dic['ppcoa'] -= 1
            
        elif each_substrate == 'ivcoa':
            temp_dic['ivcoa'] -= 1

        elif each_substrate == 'mxmalacp':
            temp_dic['mxmalacp'] -= 1
            
        elif each_substrate == 'chccoa':
            temp_dic['chccoa'] -= 1
         
        temp_list_of_reaction_set.append(temp_dic)
                
        if substrate_decision_number == 0:
            break
        
    return temp_list_of_reaction_set

def adding_product_to_the_reaction(product, list_of_reaction_set):
    
    product_count = 1
    list_of_reaction_set_with_product = []
    
    print product, list_of_reaction_set
    
    for each_integrated_reaction in list_of_reaction_set:
        
        print each_integrated_reaction
        
        new_product_name = product + '_' +  str(product_count)
        each_integrated_reaction[new_product_name] = 1
        list_of_reaction_set_with_product.append(each_integrated_reaction)
        product_count += 1
    
    print list_of_reaction_set_with_product
    
    return list_of_reaction_set_with_product

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
        converted_MNXMID = ['MNXM144', 'ala_DASH_B', 'beta-alanine', 'C00099']
    elif MetID == '24dab':
        converted_MNXMID = ['MNXM840', '24dab', 'L-2,4-diazaniumylbutyrate', 'C03283'] 
    elif MetID == 'pac':  
        converted_MNXMID = ['MNXM497', 'pac', 'phenylacetatee', 'C00548']
    elif MetID == 'tcl':
        converted_MNXMID = ['MNXM37380', 'tcl', '4-Chlorothreonine', 'N/A'] 
    elif MetID == 'qa':
        converted_MNXMID = ['MNXM4797', 'qa', 'quinaldate', 'C06325'] 
         
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

def second_metab_reactions_addition(cobra_model, product, dic_pksnrps_gene, list_of_reaction_set_with_product, metab_MNXM_dict):
#     fp1 = open('output_set_of_integrated_metabolic_reaction.txt','w')
    fp2 = open('output_participated_gene_list_of_backbone_biosynthesis.txt','w')
    fp3 = open("output_database_format_file.txt", 'w')
    product_count = 1
    
    dic_novel_nrps_cluster = {}
    list_reaction_name_SM = []
    list_novel_secondary_metabolite_reactions = []
    
#Extracting each reaction dictionary
    for each_integrated_reaction in list_of_reaction_set_with_product:
        
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
        for each_gene in dic_pksnrps_gene:
            if gpr_count == 0:
                gpr_list = each_gene 
                gpr_count += 1
            else:
                gpr_list = gpr_list + ' AND ' + each_gene
        
        print gpr_list
        print >>fp2, gpr_list 
        reaction.add_gene_reaction_rule(gpr_list)

#Adding the new reaction to the model
        cobra_model.add_reaction(reaction)
        
        reaction_name = reaction.id
        strain_name = reaction_name.split("_")
        strain_name = strain_name[0].strip()

        print "\n", "Cluster reaction:", reaction
        print "Cluster genes:", reaction.gene_reaction_rule
        print reaction.reaction
#         print >>fp1, reaction.reaction
        print >>fp3, "%s\t%s\t%s\t%s" % (strain_name, reaction, reaction.gene_reaction_rule, reaction.reaction)

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
        
#     fp1.close()
    fp2.close()
    fp3.close()
    return cobra_model, list_reaction_name_SM, list_novel_secondary_metabolite_reactions

# performing FBA for metabolic reactions of novel secondary metabolties
def performing_FBA_for_each_reaction_of_SMRs(cobra_model, list_novel_secondary_metabolite_reaction_names):
    
#     fp1 = open(outputfile, 'w')
    
    for each_MR in list_novel_secondary_metabolite_reaction_names:
       
        output = each_MR + '_max.txt'
        query_reaction = 'Transport_'+each_MR
        calculated_cobra_model = cobra_model_FBA(cobra_model, query_reaction, 'maximize', output)

#     fp1.close()
    
    return 

    
# running FBA
def cobra_model_FBA(cobra_model, objective, sense, output):

#     fp1 = open(output,"w")

    reaction = cobra_model.reactions.get_by_id(objective)

    # Run the optimization for the objective reaction and medium composition
    cobra_model.optimize(new_objective=reaction, objective_sense=sense, solver='gurobi')

#     print "\n"
    print cobra_model.solution.f
#     print cobra_model.solution.status
#     fp1.write(str(cobra_model.solution.status)+"\n")
#     fp1.write(str(cobra_model.solution.f)+"\n"+"\n")

#     for the_reaction, the_value in cobra_model.solution.x_dict.items():
#         #print '%s: %1.2f' %(the_reaction, the_value)
#         fp1.write(str(the_reaction)+"\t"+str(the_value)+"\n")

#     print "Done"
#     fp1.close()

    return cobra_model
 
       
if __name__ == '__main__':
# making template model in order to (V)
    cobra_model = create_cobra_model_from_sbml_file('SCO_model_snu.xml', print_time=True)
    inputfile = './NC_020990.1.cluster023.gbk'
    
# making the dictionary file of metabolites in template model (V)
    MNXM_dict = pickle.load(open("Pickle_MNXM_dict.p","rb"))
    
# making dictionary file of metabolites in template model (V)
# (e.g. ''MNXM37': 'gln_DASH_L_p'])
    metab_MNXM_dict = COBRA_TemplateModel_checking_MNXref_metabolites(cobra_model, MNXM_dict)
    
# Enrolling participated substrates for backbone biosynthesis in Type I polyketide (V)
    monomer_mnx_dic = monomers_list_for_nrps('Input_monomers_pksnrps.txt')
    
# extracting information of total combined (V)
    second_total_monomers = second_metab_monomers(inputfile, "genbank", monomer_mnx_dic)
    
# To generate the product name of target 'type I pks' genes (V)
    product = second_metab_reaction_product_names(inputfile, "genbank")
    
# To extract genes related to backbone biosynthesis in type I polyketide synthase (VV)
# Extracting information of locus_tag, locus_tag, substrate_inform(aSProdPred) from query genebankfile 
# dic_pksnrps_gene['SCAB_43961'] = ['polyketide synthase', '(pk-mal-nrp)']
    dic_pksnrps_gene = second_metab_genes(inputfile, "genbank")  
    
# Extracting information of domains from dictionary file 'dic_t1pks_domain' (VV)
#  dic_t1pks_domain[SAV_942_DM12] = ['PKS_KR', '(5084-5264)']
    dic_info_of_bundle_set_met = extracting_sub_set_met_info_from_genebank(inputfile, "genbank", dic_pksnrps_gene)
    
# Extracting information of domains from dictionary file 'dic_info_of_bundle_set_met' and 'dic_t1pks_gene' (VV)
# dic_t1pks_domain['SAV_942_DM6'] = ['PKS_AT', '(2425-2719)', 'SAV_942']
# dic_t1pks_gene_domain['SAV_938'] = ['PKS_AT', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_KR', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP', 'PKS_Docking_Cterm']
    dic_pksnrps_domain, dic_pksnrps_gene_domain, dic_pksnrps_gene_substrate, dic_t1pks_PKS_KR_activity = second_metab_domain(dic_info_of_bundle_set_met, dic_pksnrps_gene)
    
# Extracting information of substrates with their substrates from dictionary file 'dic_info_of_bundle_set_met' and 'dic_t1pks_gene' (VV)
# dic_t1pks_domain_substrate['SAV_943_M1'] = ['mmal', 'Ethyl_mal', 'pk']
    dic_pksnrps_domain_substrate = second_metab_substrate(dic_info_of_bundle_set_met, dic_pksnrps_gene)
   
# Extracting information of modules from dictionary file 'dic_t1pks_domain' (VV)
# dic_t1pks_module['SAV_943_M1'] = ['PKS_KS', 'PKS_AT', 'ACP']
# dic_pksnrps_module_KR_activity[''] = []
    dic_pksnrps_module, dic_pksnrps_module_KR_activity = second_metab_module(dic_pksnrps_gene, dic_pksnrps_gene_domain, dic_t1pks_PKS_KR_activity)

###@@ logic is critically changed at this point

# Generating rules for biosynthesis of type I PKS and converting module and its substrate to metabolic reactions
# For example : dic_converted_metabolic_reaction_without_substrate['SAV_943_M0'] = ['coa': 1, 'nadph': -1, 'nadp': 1, 'hco3': 1, 'h': -1]
    dic_converted_metabolic_reaction_without_substrate = generating_each_module_of_backbone_biosynthesis_for_pksnrps(dic_pksnrps_module, dic_pksnrps_module_KR_activity)
    
# generating integrated metabolic reaction without participated substrate such as malonyl-coenzyme A
    dic_integrated_metabolic_reaction_without_cofactors = integrated_metabolic_reaction1(dic_converted_metabolic_reaction_without_substrate) #####
    
# adding matched participated substrate by using 'product' and dictionary 'dic_t1pks_domain_substrate'
# dic_semiintegrated_metabolic_reaction {'coa': 13, 'mmalcoa': -4, 'h': -10, 'malcoa': -7, 'hco3': 13, 'nadph': -10, 'h2o': 5, 'nadp': 10}
# list_of_dismatched_substrate = [['mmal', 'Ethyl_mal'], ['2metbut', '2metbut']]
    dic_semiintegrated_metabolic_reaction, list_of_dismatched_substrate = integrated_metabolic_reaction2(dic_pksnrps_domain_substrate, dic_integrated_metabolic_reaction_without_cofactors)
    
# completing integrated metabolic reaction by adding product and dismatched substrate to the reaction.
# list_of_reaction_set = [{'nadph': -10, 'nadp': 10, 'ahcys': 0, '2mbcoa': -1, 'nad': 0, 'h': -10, 'fadh2': 0, 'malcoa': -7, 'hco3': 13, 'amet': 0, 'coa': 13, 'h2o': 5, 'nadh': 0, '13dpg': 0, 'mmalcoa': -5, 'pi': 0, 'emalcoa': 0, 'fad': 0}, ...]
    list_of_reaction_set = integrated_metabolic_reaction3(dic_semiintegrated_metabolic_reaction, list_of_dismatched_substrate)

# adding product name to the integrated reaction
# list_of_reaction_set_with_product = {'nadph': -10, 'nadp': 10, 'ahcys': 0, '2mbcoa': -1, 'nad': 0, 'h': -10, 'fadh2': 0, 'malcoa': -7, 'hco3': 13, 'Cluster_05_t1pks_1': 1, 'amet': 0, 'coa': 13, 'h2o': 5, 'nadh': 0, '13dpg': 0, 'mmalcoa': -5, 'pi': 0, 'emalcoa': 0, 'fad': 0}] 
    list_of_reaction_set_with_product = adding_product_to_the_reaction(product, list_of_reaction_set)
    
# making all metabolic reactions from each putative type 1 PKS gene clusters.
# dic_reaction_info_set[gene_info] = {reaction_set}
#     dic_converted_metabolic_reactions_with_gene = making_metabolic_reactions_from_gene_cluster_data(dic_t1pks_gene, list_of_reaction_set_with_product)
         
# making total list of metabolites in model (e.g. metab_MNXM_dict[MNX_ID(e.g. MNXM37)] = metabolite_name_in_model(e.g. gln_DASH_L_p))
#     metab_MNXM_dict = COBRA_TemplateModel_checking_MNXref_metabolites(cobra_model, MNXM_dict)
     
# adding metabolic reactions to the model
    modified_cobra_model, list_reaction_name_SM, list_novel_secondary_metabolite_reactions = second_metab_reactions_addition(cobra_model, product, dic_pksnrps_gene, list_of_reaction_set_with_product, metab_MNXM_dict)
          
## simulating FBA 
    performing_FBA_for_each_reaction_of_SMRs(modified_cobra_model, list_reaction_name_SM)   
    


