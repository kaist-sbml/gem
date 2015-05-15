'''
Created by Kyu-Sang Hwang, Sept 2014

This file generates metabolic reactions for the genes newly annotated 
to be present in the secondary metabolite-biosynthetic gene cluster from antiSMASH.
'''

from time import time
import cobra
# from phenotype_phase_plane import calculate_phenotype_phase_plane
from Bio import SeqIO
from sets import Set
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file
from MNX_checker2 import COBRA_TemplateModel_checking_MNXref_metabolites
from MNX_checker2 import fix_legacy_id
from Multi_objective_LP_v2 import cobra_2D_optimization
from Multi_objective_LP_v2 import cobra_3D_optimization
# from run_cobra_Strepto import cobra_model_FBA
# import utils
import pickle
import copy
import COBRA_GeneTargeting_v1_3

def monomers_list_for_typeIPKS(inputFile):
    fp1 = open(inputFile,"r")
    monomer_mnx_dic = {}

    monomer = fp1.readline()

    while monomer:
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

    dic_t1pks_gene = {}
    list_t1pks_gene = []

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
                        
            if 'Type: t1pks' in qualifier_np_type:
                dic_t1pks_gene[locus_tag] = [product, order_of_substrates]
                print >>fp1, "%s\t%s\t%s\t%s" % (whole_cluster_name, locus_tag, product, order_of_substrates)
                print whole_cluster_name, locus_tag, product, order_of_substrates
                      
    fp1.close()
    return dic_t1pks_gene


def extracting_sub_set_met_info_from_genebank(gbkFile, FileType, dic_t1pks_gene):

    dic_info_of_bundle_set_met = {}
     
#Reads GenBank file
    record = SeqIO.read(gbkFile, FileType)
    
    for dic_gene_key in dic_t1pks_gene:
        
        list_t1pks_domain = []
        
        for feature in record.features: 
        
            if feature.type == 'CDS':
        
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')
                        
                if qualifier_locus_tag[0] == dic_gene_key:    

                    qualifier_set_met = feature.qualifiers.get('sec_met')
                    list_t1pks_domain.append(qualifier_set_met)
        
        dic_info_of_bundle_set_met[dic_gene_key] = list_t1pks_domain
                    
    return dic_info_of_bundle_set_met

def second_metab_domain(dic_info_of_bundle_set_met, dic_t1pks_gene):
    fp1 = open('Output_second_metab_gene_domain.txt','w')
    fp2 = open('Output_second_metab_gene_substrate.txt','w')
    fp3 = open('Output_second_metab_gene_KR_activity.txt','w')
    
    dic_t1pks_domain = {}
    dic_t1pks_gene_domain = {}
    dic_t1pks_gene_substrate = {}
    dic_t1pks_PKS_KR_activity = {}
        
    for each_gene in dic_t1pks_gene:
        
        list_set_met = dic_info_of_bundle_set_met[each_gene][0]

        domain_count = 0
        list_t1pks_domain = []
        list_t1pks_gene_substrate = []
        list_t1pks_PKS_KR_activity = []
    
        for each_sub_set in list_set_met:

            if "NRPS/PKS Domain" in each_sub_set:
                            
                sptline1 = each_sub_set.split('. ')
                crude_domain_info = sptline1[0] 
                
                # to extract the information of KR activity and AT substrate specificity
                sptline2 = each_sub_set.split('; ')
                
                if "PKS_AT" in crude_domain_info:
                    spt_PKS_AT = sptline2[1].split(':')
                    spt_PKS_AT_met = spt_PKS_AT[1].split(', ')
                    spt_PKS_AT_met = spt_PKS_AT_met[2].split(' (')
                    spt_PKS_AT_met = spt_PKS_AT_met[0].strip()
                    
                    list_t1pks_gene_substrate.append(spt_PKS_AT_met)
                
                if "PKS_KR" in crude_domain_info:
                    spt_PKS_KR = sptline2[1].split(': ')
                    spt_PKS_KR = spt_PKS_KR[1].strip()
                    
                    list_t1pks_PKS_KR_activity.append(spt_PKS_KR)
                    
                spt_domain_info = crude_domain_info.split(':')
                whole_domain_info = spt_domain_info[1]
                            
                spt_list_domain_info = whole_domain_info.split()
                spt_list_domain_info.append(each_gene)
                            
                each_t1pks_domain = spt_list_domain_info[0]
                list_t1pks_domain.append(each_t1pks_domain)
                
                # domain information
                            
                domain_number = each_gene + '_DM' + str(domain_count)
                domain_count = domain_count + 1
                            
                dic_t1pks_domain[domain_number] = spt_list_domain_info
                print domain_number, spt_list_domain_info
                    
        dic_t1pks_gene_domain[each_gene] = list_t1pks_domain
        dic_t1pks_gene_substrate[each_gene] = list_t1pks_gene_substrate
        dic_t1pks_PKS_KR_activity[each_gene] = list_t1pks_PKS_KR_activity
        
        print each_gene, list_t1pks_domain
        print each_gene, list_t1pks_gene_substrate
        print each_gene, list_t1pks_PKS_KR_activity

        print >>fp1, "%s\t%s" % (each_gene, list_t1pks_domain)
        print >>fp2, "%s\t%s" % (each_gene, list_t1pks_gene_substrate)
        print >>fp3, "%s\t%s" % (each_gene, list_t1pks_PKS_KR_activity)

    fp1.close()
    fp2.close()
    fp3.close()      
    return dic_t1pks_domain, dic_t1pks_gene_domain, dic_t1pks_gene_substrate, dic_t1pks_PKS_KR_activity

def second_metab_substrate(dic_info_of_bundle_set_met, dic_t1pks_gene):
    fp1 = open('Output_second_metab_substrate_with_domain.txt','w')

    dic_t1pks_domain_substrate = {}
    
    for each_gene in dic_t1pks_gene:

        module_count = 0
        list_set_met =  dic_info_of_bundle_set_met[each_gene][0]

        for each_sub_set in list_set_met:
     
            if "Substrate specificity predictions" in each_sub_set:
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
                                
                module_number = each_gene + '_M' + str(module_count)
                                
                dic_t1pks_domain_substrate[module_number] = list_participated_sustrate
                
                print module_number, list_participated_sustrate
                print >>fp1, "%s\t%s" % (module_number, list_participated_sustrate)
                    
                module_count = module_count + 1
    
    for each_key in dic_t1pks_domain_substrate:          
        print each_key, dic_t1pks_domain_substrate[each_key]
    
    fp1.close()
    return dic_t1pks_domain_substrate


def second_metab_module(dic_t1pks_gene ,dic_t1pks_gene_domain, dic_t1pks_PKS_KR_activity):
    fp1 = open('Output_second_metab_module.txt','w')
    
    dic_t1pks_module = {}
    
    for t1pks_gene in dic_t1pks_gene:
        
        count = 0
  
        if t1pks_gene in dic_t1pks_gene_domain:
            
            list_each_t1pks_domain = dic_t1pks_gene_domain[t1pks_gene]

            list_KR_activity = dic_t1pks_PKS_KR_activity[t1pks_gene]
        
            list_module_info = []
            
            number_of_list = len(list_each_t1pks_domain)
            
            KR_number = 0
            
            for each_domain in list_each_t1pks_domain:
                
                if each_domain == 'PKS_Docking_Nterm' or each_domain == 'PKS_Docking_Cterm':
                    number_of_list = number_of_list - 1
                    continue
                
                list_module_info.append(each_domain)
                number_of_list = number_of_list - 1
                
                if each_domain == 'ACP' or each_domain == 'PCP':
                    
                    module_number = t1pks_gene + '_M' + str(count)
                    
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
                                        
                    dic_t1pks_module[module_number] = [list_module_info, each_module_KR_activity]
                    
                    print >>fp1, "%s\t%s\t%s\t%s" % (t1pks_gene, module_number, list_module_info, each_module_KR_activity)
                    print t1pks_gene, module_number, dic_t1pks_module[module_number]

                    list_module_info = []
                    
                    count = count + 1
                    
                elif each_domain == 'Thioesterase':
                    terminated_domain = 'TE'
                    
                    count = count - 1
                    
                    module_number = t1pks_gene + '_M' + str(count)
                    
                    
                    list_module_info = dic_t1pks_module[module_number][0]
                    KR_activity_in_TE_module = dic_t1pks_module[module_number][1]
                    print list_module_info, KR_activity_in_TE_module
                    
                    list_module_info.append(terminated_domain)
                    print list_module_info
                    
                    A = dic_t1pks_module.pop(module_number)
                    
                    print A

                    dic_t1pks_module[module_number] = [list_module_info, KR_activity_in_TE_module]
                    print module_number, dic_t1pks_module[module_number]
                    
#                   dic_t1pks_module[module_number].append(each_domain)
                    
                    print >>fp1, "%s\t%s\t%s\t%s" % (t1pks_gene, module_number, list_module_info, each_module_KR_activity)
                    
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
                    
                    dic_t1pks_module[module_number] = [list_module_info, each_module_KR_activity]
                    print >>fp1, "%s\t%s\t%s\t%s" % (t1pks_gene, module_number, list_module_info, each_module_KR_activity)

                    list_module_info = []

                    list_module_info.append(poped_domain)   
                    
                    count = count + 1
                    
                elif float(number_of_list) == 0:
                    
                    module_number = t1pks_gene + '_M' + str(count)
                    
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
                    
                    dic_t1pks_module[module_number] = [list_module_info, each_module_KR_activity]
                    print >>fp1, "%s\t%s\t%s\t%s" % (t1pks_gene, module_number, list_module_info, each_module_KR_activity)
                    print module_number, dic_t1pks_module[module_number]

                    list_module_info = []
          
                    count = count + 1
    
    for each_key in dic_t1pks_module:
        print "module number : %s, module_composition : %s" % (each_key, dic_t1pks_module[each_key][0])
                
    fp1.close()
    
    return dic_t1pks_module
#######################################################################

def generating_each_module_of_backbone_biosynthesis_for_t1pks(dic_t1pks_module):
    fp1 = open('Output_backbone_biosynthesis_of_each_module.txt','w')
    
    dic_converted_metabolic_reaction_without_substrate = {}

    for each_module in dic_t1pks_module:
        
        domain_comb = dic_t1pks_module[each_module][0]
        module_KR_activity = dic_t1pks_module[each_module][1]
        
        print domain_comb, module_KR_activity

        each_module_substrates = {}
        
        discriminant = module_discriminator(domain_comb)
        discriminant_with_KRact = Identifier_KR_activity(discriminant, module_KR_activity)
        
        if discriminant_with_KRact == 'None':
            print "this discriminant is not defined : %s" % (domain_comb)
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, 'None')
            continue
 
        if discriminant_with_KRact == 'AT_ACP' or discriminant_with_KRact == 'AT_KR(inactive)_ACP' or discriminant_with_KRact == 'AT_DH_KR(inactive)_ACP' or discriminant_with_KRact == 'AT_DH_ER_KR(inactive)_ACP': 
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates     
            print 'reaction 1: AT-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant_with_KRact == 'AT_KS_ACP' or discriminant_with_KRact == 'AT_KS' or discriminant_with_KRact == 'AT_KS_KR(inactive)_ACP' or discriminant_with_KRact == 'AT_KS_KR(inactive)' or discriminant_with_KRact == 'AT_KS_DH_KR(inactive)_ACP' or discriminant_with_KRact == 'AT_KS_DH_KR(inactive)' or discriminant_with_KRact == 'AT_KS_DH_KR(inactive)_ACP' or discriminant_with_KRact == 'AT_KS_DH_ER_KR(inactive)': 
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates     
            print 'reaction 2: AT-KS-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
             
        elif discriminant_with_KRact == 'AT_KS_KR_ACP' or discriminant_with_KRact == 'AT_KR_ACP' or discriminant_with_KRact == 'AT_KS_KR' or discriminant_with_KRact == 'AT_ER_KR_ACP' or discriminant == 'AT_KS_ER_KR' or discriminant == 'AT_KS_ER_KR_ACP':
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 3: KS-AT-KR-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
             
        elif discriminant_with_KRact == 'AT_KS_DH_KR_ACP' or discriminant_with_KRact == 'AT_DH_KR_ACP' or discriminant_with_KRact == 'AT_KS_DH_KR':
            each_module_substrates['coa'] = 1 #'Coenzyme A', (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['h2o'] = 1 #H2O (C00001)
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 4: KS-AT-DH-KR-ACP', each_module_substrates
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates) 
             
        elif discriminant_with_KRact == 'AT_KS_DH_ER_KR_ACP' or discriminant_with_KRact == 'AT_DH_ER_KR_ACP' or discriminant_with_KRact == 'AT_KS_DH_ER_KR':
            each_module_substrates['coa'] = 1 #'Coenzyme A', (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['h2o'] = 1 #'H2O', (C00001), 'MNXM2'
            each_module_substrates['nadp'] = 2 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -2 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -2 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates 
            print 'reaction 5: KS-AT-DH-ER-KR-ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant_with_KRact == 'ACP':
            continue
            
        elif discriminant_with_KRact == 'AT_KS_cMT_ACP' or discriminant_with_KRact == 'AT_KS_cMT':
            each_module_substrates['coa'] = 1 #'Coenzyme A' (C00010), 'MNXM12'
            each_module_substrates['hco3'] = 1 #'bicarbonate', (C00288), 'MNXM60'
            each_module_substrates['amet'] = -1 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
            each_module_substrates['ahcys'] = 1 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 6: AT_KS_cMT_ACP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant_with_KRact == 'AT_KS_KR_cMT_ACP' or discriminant_with_KRact == 'AT_KS_KR_cMT':
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
            
        elif discriminant_with_KRact == 'AT_KS_DH_KR_cMT_ACP' or discriminant_with_KRact == 'AT_KS_DH_KR_cMT':
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
            
        elif discriminant_with_KRact == 'AT_KS_DH_KR_cMT_ER_ACP' or discriminant_with_KRact == 'AT_KS_DH_KR_cMT_ER':
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
            
        elif discriminant_with_KRact == 'TD':
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 8: TD', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
              
    fp1.close()
                
    return dic_converted_metabolic_reaction_without_substrate

def module_discriminator(domain_comb):
    
    print domain_comb
# Starter units 
    if 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
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
        
# terminal region
    elif 'TD' in domain_comb and 'PKS_AT' not in domain_comb and 'PKS_KS' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'TD'
        
    else:
        discriminant = 'None'
    
    return discriminant

def Identifier_KR_activity(discriminant, each_module_KR_activity):
    
    print discriminant, each_module_KR_activity

    if each_module_KR_activity == 1 or each_module_KR_activity == 0:
        discriminant_with_KRact = discriminant
        return discriminant_with_KRact
    
    elif each_module_KR_activity == 2:
        discriminant_with_KRact = discriminant.replace('KR','KR(inactive)')      
        return discriminant_with_KRact

def integrated_metabolic_reaction1(participated_cofactor_info):
    fp1 = open('Output_participated_cofactors.txt','w')
    
    dic_integrated_metabolic_reaction = {}
    dic_integrated_metabolic_reaction['h'] = 0 
    dic_integrated_metabolic_reaction['h2o'] = 0
    dic_integrated_metabolic_reaction['nadp'] = 0
    dic_integrated_metabolic_reaction['nadph'] = 0
    dic_integrated_metabolic_reaction['coa'] = 0
    dic_integrated_metabolic_reaction['hco3'] = 0
    dic_integrated_metabolic_reaction['ahcys'] = 0  #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
    dic_integrated_metabolic_reaction['amet'] = 0   #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
    
    for each_module in participated_cofactor_info:
        print each_module
        for each_metabolite in participated_cofactor_info[each_module]:
            met_coff = participated_cofactor_info[each_module][each_metabolite]       
            if participated_cofactor_info[each_module][each_metabolite] > 0:
                dic_integrated_metabolic_reaction[each_metabolite] += met_coff
                print dic_integrated_metabolic_reaction
                print >>fp1, "participated_cofactors:\t%s\t%s\t%s" % (each_module, each_metabolite, met_coff)
            else:
                dic_integrated_metabolic_reaction[each_metabolite] += met_coff
                print >>fp1, "participated_cofactors:\t%s\t%s\t%s" % (each_module, each_metabolite, met_coff)
    
  
    print >>fp1, "#####"
    for each_participated_cofactor in dic_integrated_metabolic_reaction:
        print >>fp1, "total_participated_cofactor:\t%s\t%s" % (each_participated_cofactor, dic_integrated_metabolic_reaction[each_participated_cofactor])
       
    fp1.close()
    
    return dic_integrated_metabolic_reaction

def integrated_metabolic_reaction2(dic_t1pks_domain_substrate, dic_integrated_metabolic_reaction_without_cofactors):
    fp1 = open('Output_participated_substrates_from_uniformed_prediction.txt','w')
    
    list_of_dismatched_substrate = []
    
    dic_integrated_metabolic_reaction_without_cofactors['malcoa'] = 0 #'malonyl-CoA', 'C00083', 'MNXM40'
    dic_integrated_metabolic_reaction_without_cofactors['mmcoa_DASH_S'] = 0 #'(S)-methylmalonyl-CoA(5-)','C00683', 'MNXM190', 'not detected in bigg database'
    
    for each_module in dic_t1pks_domain_substrate:
        # CAL_domain : NH2 (Minowa)
        if len(dic_t1pks_domain_substrate[each_module]) < 3:
            continue
        # seperating between unified prediction and not
        if dic_t1pks_domain_substrate[each_module][2] == 'pk':
            temp_list = []
            temp_list.append(dic_t1pks_domain_substrate[each_module][0])
            temp_list.append(dic_t1pks_domain_substrate[each_module][1])
            list_of_dismatched_substrate.append(temp_list)
        
        elif dic_t1pks_domain_substrate[each_module][2] == 'mal':
            dic_integrated_metabolic_reaction_without_cofactors['malcoa'] -= 1
        
        elif dic_t1pks_domain_substrate[each_module][2] == 'mmal':
            dic_integrated_metabolic_reaction_without_cofactors['mmcoa_DASH_S'] -= 1 
    
    for each_participated_substrate in dic_integrated_metabolic_reaction_without_cofactors:
        print >>fp1, "total_participated_cofactor_with_substrates1:\t%s\t%s" % (each_participated_substrate, dic_integrated_metabolic_reaction_without_cofactors[each_participated_substrate])
    
    fp1.close()
    
    return dic_integrated_metabolic_reaction_without_cofactors, list_of_dismatched_substrate

def integrated_metabolic_reaction3(dic_integrated_metabolic_reaction, list_of_dismatched_substrate):
    dic_integrated_metabolic_reaction['2mbcoa'] = 0 #'2-methylbutanoyl-CoA', C01033,'MNXM569'
    dic_integrated_metabolic_reaction['emcoa_DASH_S'] = 0 #'ethylmalonyl-CoA','C18026', 'MNXM2043', 'not detected in bigg database'
    
    dic_integrated_metabolic_reaction['ibcoa'] = 0 #'2-Methylpropanoyl-CoA', 'C00630', 'MNXM470'
    dic_integrated_metabolic_reaction['accoa'] = 0 #'Acetyl-CoA', 'C00024', 'MNXM21'
    dic_integrated_metabolic_reaction['ppcoa'] = 0 #'Propionyl-CoA', 'C00100', 'MNXM86'
    dic_integrated_metabolic_reaction['ivcoa'] = 0 #'3-Methylbutanoyl-CoA', 'C02939', 'MNXM471'
    dic_integrated_metabolic_reaction['chccoa'] = 0 #'cyclohexane-1-carboxyl-CoA', 'C09823', 'MNXM5111'
    
    dic_integrated_metabolic_reaction['mxmalacp'] = 0 #'Methoxymalonyl-[acp]', 'C18616', 'MNXM61686'
    dic_integrated_metabolic_reaction['13dpg'] = 0 #'3-phosphonato-D-glyceroyl phosphate', 'C00236', 'MNXM261'
    dic_integrated_metabolic_reaction['pi'] = 0 #'phosphate', 'C00009', 'MNXM9'
#     dic_integrated_metabolic_reaction['amet'] = 0 #'S-adenosyl-L-methionine', 'C00019', 'MNXM16'
#     dic_integrated_metabolic_reaction['ahcys'] = 0 #'S-adenosyl-L-homocysteine', 'C00021', 'MNXM19'
    dic_integrated_metabolic_reaction['nad'] = 0 #'NAD(+)', 'C00003', 'MNXM8'
    dic_integrated_metabolic_reaction['nadh'] = 0 #'NADH', 'C00004', 'MNXM10'
    dic_integrated_metabolic_reaction['fad'] = 0 #'FAD', 'C00016', 'MNXM96415'
    dic_integrated_metabolic_reaction['fadh2'] = 0 #'FADH2', 'C00016', 'MNXM96415'
    
    list_of_reaction_set = []

    if list_of_dismatched_substrate == []:
        
        list_of_reaction_set.append(dic_integrated_metabolic_reaction)
        
    else:
        
        list_of_reaction_set.append(dic_integrated_metabolic_reaction)
        
        for each_pair_of_substrates in list_of_dismatched_substrate:
                
            template_list = copy.deepcopy(list_of_reaction_set)
            list_of_reaction_set = []
        
            for dic_each_metabolic_reaction in template_list:
                
                substrate_decision_number = distincting_each_substrate_in_list_component(each_pair_of_substrates)
                temp_reaction_set = converting_pk_substrates(each_pair_of_substrates, dic_each_metabolic_reaction, substrate_decision_number) 
                
                for each_dic_reaction_set in temp_reaction_set:
                    
                    list_of_reaction_set.append(each_dic_reaction_set)  
                    
    return list_of_reaction_set
                
def distincting_each_substrate_in_list_component(each_pair_of_substrates):  
    
    print each_pair_of_substrates
## unifying met_name as emal        
    if each_pair_of_substrates[1] == 'Ethyl_mal':
        each_pair_of_substrates[1] = 'emal'

## checking_similarity of pair of substrates
    if each_pair_of_substrates[0] != each_pair_of_substrates[1]:
        substrate_decision_number = 1
#     elif each_pair_of_substrates[0] == ('emal' or 'Ethyl_mal') and  each_pair_of_substrates[1] == ('emal' or 'Ethyl_mal'):
#         substrate_decision_number = 0
    else:
        substrate_decision_number = 0
        
    return substrate_decision_number

def converting_pk_substrates(each_pair_of_substrates, dic_integrated_metabolic_reaction, substrate_decision_number):
    
    temp_list_of_reaction_set = []
    
    for each_substrate in each_pair_of_substrates:
        temp_dic ={}
        temp_dic = copy.deepcopy(dic_integrated_metabolic_reaction)
        
        if each_substrate == 'mal':
            temp_dic['malcoa'] -= 1
            
        elif each_substrate == 'mmal':
            temp_dic['mmcoa_DASH_S'] -= 1
            
        elif each_substrate == '2metbut':
            temp_dic['2mbcoa'] -= 1
            
        elif each_substrate == 'Ethyl_mal' or each_substrate == 'emal':
            temp_dic['emcoa_DASH_S'] -= 1
            
        elif each_substrate == 'isobut':
            temp_dic['ibcoa'] -= 1
            
        elif each_substrate == 'ace':
            temp_dic['accoa'] -= 1
            
        elif each_substrate == 'pro':
            temp_dic['ppcoa'] -= 1
            
        elif each_substrate == '3metbut':
            temp_dic['ivcoa'] -= 1
            
        elif each_substrate == 'CHC-CoA':
            temp_dic['chccoa'] -= 1      
            
            ## mxmal = methoxymalonyl-[ACP] (C18616, MNXM61686)??    
        elif each_substrate == 'mxmal':
            temp_dic['mxmalacp'] -= 1
#         elif each_substrate == 'mxmal':
#             temp_dic['13dpg'] -= 1 #C00236, 1,3-bisphospho-D=glycerate (V)
#             temp_dic['pi'] +=2 #C00009, phosphate (V)
#             temp_dic['amet'] -= 1 #C00019, S-adenosyl-L-methionine (V)
#             temp_dic['ahcys'] += 1 #C00021, S-Adenosyl-L-homocysteine, 'MNXM19' (V)
#             temp_dic['nad'] -= 1 #C00003, NAD+ (V), 'MNXM8'
#             temp_dic['nadh'] += 1 #C00004, NADH (V), 'MNXM10'
#             temp_dic['h'] += 1 #C00080, H+, 'MNXM1' 
#             temp_dic['fad'] -= 1 #C00016, FAD
#             temp_dic['fadh2'] += 1 #C01352, FADH2
            
        elif each_substrate == 'N/A':
            continue
            
        else:
            print each_substrate
            raw_input('substrate_not_defined')
        
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
    
    return list_of_reaction_set_with_product

def converting_MNXMID_to_biggid(MetID):
    converted_MNXMID = []
    
    if MetID == 'coa':
        converted_MNXMID = ['MNXM12', 'coa', 'Coenzyme A', 'C00010']
    elif MetID == 'hco3':
        converted_MNXMID = ['MNXM60', 'hco3', 'bicarbonate', 'C00288']
    elif MetID == 'nadp':
        converted_MNXMID = ['MNXM5', 'nadp', 'NADP+', 'C00006'] 
    elif MetID == 'nadph':
        converted_MNXMID = ['MNXM6', 'nadph', 'NADPH' ,'C00005']
    elif MetID == 'h':
        converted_MNXMID = ['MNXM1', 'h', 'H+', 'C00080']
    elif MetID == 'h2o': 
        converted_MNXMID = ['MNXM2', 'h2o', 'H2O', 'C00001']
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
    elif MetID == 'chccoa':
        converted_MNXMID = ['MNXM5111','chccoa', 'cyclohexane-1-carboxyl-CoA', 'C09823']
    elif MetID == 'mxmalacp':
        converted_MNXMID = ['MNXM61686','mxmalacp', 'Methoxymalonyl-[acp]', 'C18616']        

    elif MetID == '13dpg':
        converted_MNXMID = ['MNXM261', '13dpg', '3-phosphonato-D-glyceroyl phosphate', 'C00236']
    elif MetID == 'pi':
        converted_MNXMID = ['MNXM9', 'pi', 'phosphate', 'C00009']
    elif MetID == 'amet':
        converted_MNXMID = ['MNXM16', 'amet', 'S-adenosyl-L-methionine', 'C00019']
    elif MetID == 'ahcys':
        converted_MNXMID = ['MNXM19', 'ahcys', 'S-adenosyl-L-homocysteine', 'C00021']
    elif MetID == 'nad':
        converted_MNXMID = ['MNXM8', 'nad', 'NAD(+)', 'C00003']
    elif MetID == 'nadh':
        converted_MNXMID = ['MNXM10', 'nadh', 'NADH', 'C00004']
    elif MetID == 'fad':
        converted_MNXMID = ['MNXM96415', 'fad', 'FAD', 'C00016']
    elif MetID == 'fadh2':
        converted_MNXMID = ['MNXM38', 'fadh2', 'FADH2', 'C00016']
    else:
        converted_MNXMID = ['', MetID, '', '']
    
    return converted_MNXMID

def second_metab_reactions_addition(cobra_model, product, dic_t1pks_gene, list_of_reaction_set_with_product, outputfile1, outputfile2):
#     fp1 = open('output_set_of_integrated_metabolic_reaction.txt','w')
    fp2 = open(outputfile1,'w')
    fp3 = open(outputfile2, 'w')
    product_count = 1
    
    dic_novel_t1pks_cluster = {}
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
            
            each_metabolite = each_metabolite + "_c"
            
            if each_metabolite in cobra_model.metabolites:
                reaction.add_metabolites({cobra_model.metabolites.get_by_id(each_metabolite):each_integrated_reaction[abbr_MetID]})
            elif each_integrated_reaction[abbr_MetID] == 0:
                continue
            else:
                obj_each_metabolite = Metabolite(abbr_MetID+'_c', name=Met_Name, compartment='c')
                reaction.add_metabolites({obj_each_metabolite:each_integrated_reaction[abbr_MetID]})

#Setting GPR association
        gpr_count = 0
        for each_gene in dic_t1pks_gene:
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
    return cobra_model, list_novel_secondary_metabolite_reactions, list_reaction_name_SM

def exports_reactions_addition(cobra_model, list_export_met_name, outputfile1, outputfile2):
    
    fp1 = open(outputfile1, 'w')
    fp2 = open(outputfile2, 'w')
    
    for each_export_met_name in list_export_met_name:

#Creating a transport reaction
#Creating reaction ID
        TransportR = "Transport_" + each_export_met_name
        Cofactor_transport_reaction = TransportR
        print >>fp1, "%s" % TransportR
        print TransportR
        reaction = Reaction(TransportR)

#Setting bounds
        reaction.reversibility = 0 # 1: reversible
        reaction.lower_bound = 0
        reaction.upper_bound = 999999

#Adding a substrate metabolite
#    print cobra_model.metabolites.get_by_id(str(product_c))
        reaction.add_metabolites({cobra_model.metabolites.get_by_id(str(each_export_met_name+'_c')):-1})

#Adding product metabolite(s)
        new_product_name_e = Metabolite(each_export_met_name+"_e", name='', compartment='e')
        reaction.add_metabolites({new_product_name_e:1})

#Adding the new reaction to the model
        cobra_model.add_reaction(reaction)

        print "\n", "Transport reaction:", reaction
        print reaction.reaction


#Creating an exchange reaction
#Creating reaction ID
        ExportR = "Ex_" + each_export_met_name
        Cofactor_export_reaction = ExportR
        print >>fp2, "%s" % ExportR
        reaction = Reaction(ExportR)

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
    
    fp1.close()
    fp2.close()
    
    return cobra_model, Cofactor_transport_reaction, Cofactor_export_reaction

def generating_control_reactions_against_NovelSMR(list_reaction_name_SM, list_novel_secondary_metabolite_reactions, nadph_coef_number='None', amet_coef_number='None', fadh2_coef_number='None', nadh2_coef_number='None'):
    
    list_of_completed_reaction_set_without_cofactor = []
    list_reaction_name_SM_without_cofactor = []
    count = 0
    
    each_reaction_name = list_reaction_name_SM[count]
    modified_product = each_reaction_name[:-1] + "x_cof"
    
    for each_reaction_set in list_novel_secondary_metabolite_reactions:
        
        each_reaction_name = list_reaction_name_SM[count]
        count = count + 1
        
        temp_dic = {}
        temp_dic = copy.deepcopy(each_reaction_set)
    
        if nadph_coef_number != 'None':
            temp_dic['nadph'] = -(nadph_coef_number)
            temp_dic['nadp'] = +(nadph_coef_number)
            temp_dic['pi'] = +(nadph_coef_number)
        
        if amet_coef_number != 'None':
            temp_dic['amet'] = -(amet_coef_number)
            temp_dic['ahcys'] = +(amet_coef_number)
        
        if fadh2_coef_number != 'None':
            temp_dic['fadh2'] = +(fadh2_coef_number)
            temp_dic['fad'] = -(fadh2_coef_number)
        
        if nadh2_coef_number != 'None':
            temp_dic['nadh'] = +(nadh2_coef_number)
            temp_dic['nad'] = -(nadh2_coef_number)
            
        if nadph_coef_number != 'None' or amet_coef_number != 'None' or fadh2_coef_number != 'None' or nadh2_coef_number != 'None':
            temp_dic['h'] = 0
#             temp_dic['hco3'] = 0
            temp_dic['h2o'] = 0
            temp_dic.pop(each_reaction_name)
            if float(count) < 10:
                modified_met_name = each_reaction_name[:-1] + "x_cof_" + each_reaction_name[-1:]
            else:
                modified_met_name = each_reaction_name[:-2] + "x_cof_" + each_reaction_name[-2:]
            temp_dic[modified_met_name] = 1
#             list_reaction_name_SM_without_cofactor.append(modified_met_name)
            list_of_completed_reaction_set_without_cofactor.append(temp_dic)
    
    return modified_product, list_of_completed_reaction_set_without_cofactor

def removing_any_metabolic_reaction_from_model(cobra_model, reaction_list):
    
    for each_reaction in reaction_list:
        
        removed_MR = cobra_model.reactions.get_by_id(each_reaction)
        cobra_model.remove_reactions(removed_MR)
    
    return cobra_model
    
# performing FBA for metabolic reactions of novel secondary metabolties
def performing_FBA_for_each_reaction_of_SMRs(cobra_model, list_novel_secondary_metabolite_reaction_names):
    
#     fp1 = open(outputfile, 'w')
    
    for each_MR in list_novel_secondary_metabolite_reaction_names:
       
        output = each_MR + '_max.txt'
        query_reaction = 'Transport_'+each_MR
        cobra_model = cobra_model_FBA(cobra_model, query_reaction, 'maximize', output)
        
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
# 
#     for the_reaction, the_value in cobra_model.solution.x_dict.items():
#         #print '%s: %1.2f' %(the_reaction, the_value)
#         fp1.write(str(the_reaction)+"\t"+str(the_value)+"\n")
# 
#     print "Done"
#     fp1.close()

    return cobra_model

# 'Biomass_SCO'
# 

def performing_2D_optimization_between_product_biomass(cobra_model, biomassR_name ,list_novel_SMR_names, start_number, spt_max_number):
    
    for each_reaction_name in list_novel_SMR_names:
        
        output = each_reaction_name + '_2D_opt.txt'
        
        x_axis_reaction = each_reaction_name
        y_axis_reaction = biomassR_name
        
        cobra_model2 = cobra_2D_optimization(cobra_model, y_axis_reaction, x_axis_reaction, start_number, spt_max_number, output)
        
    return 

def performing_3D_optimization_between_product_cofactor_and_biomass(cobra_model, biomassR_name , cofactorR_name ,list_novel_SMR_names, start_number, spt_max_number):
    
    for each_reaction_name in list_novel_SMR_names:
        
        output = each_reaction_name + '_3D_opt.txt'
        
        x_axis_reaction = each_reaction_name
        y_axis_reaction = cofactorR_name
        z_axis_reaction = biomassR_name
        
        cobra_model2 = cobra_3D_optimization(cobra_model, z_axis_reaction, y_axis_reaction, x_axis_reaction, start_number, spt_max_number, output)
        
    return cobra_model2

       
if __name__ == '__main__':
# making template model in order to (V)
    cobra_model = create_cobra_model_from_sbml_file('SCO_model_snu.xml', print_time=True)
    inputfile = 'NC_013929.1.cluster031.gbk'
    
# making the dictionary file of metabolites in template model (V)
    MNXM_dict = pickle.load(open("Pickle_MNXM_dict.p","rb"))
    
# making dictionary file of metabolites in template model (V)
# (e.g. ''MNXM37': 'gln_DASH_L_p'])
    metab_MNXM_dict = COBRA_TemplateModel_checking_MNXref_metabolites(cobra_model, MNXM_dict)
    
# Enrolling participated substrates for backbone biosynthesis in Type I polyketide (V)
    monomer_mnx_dic = monomers_list_for_typeIPKS('Input_monomers_pks.txt')
    
# extracting information of total combined (V)
    second_total_monomers = second_metab_monomers(inputfile, "genbank", monomer_mnx_dic)
    
# To generate the product name of target 'type I pks' genes (V)
    product = second_metab_reaction_product_names(inputfile, "genbank")
    
# To extract genes related to backbone biosynthesis in type I polyketide synthase (VV)
# Extracting information of locus_tag, locus_tag, substrate_inform(aSProdPred) from query genebankfile 
# dic_t1pks_gene['SAV_938'] = ['type I polyketide synthase AVES 1', 'pk-mmal-mal']
    dic_t1pks_gene = second_metab_genes(inputfile, "genbank")  
    
# Extracting information of domains from dictionary file 'dic_t1pks_domain' (VV)
#  dic_t1pks_domain[SAV_942_DM12] = ['PKS_KR', '(5084-5264)']
    dic_info_of_bundle_set_met = extracting_sub_set_met_info_from_genebank(inputfile, "genbank", dic_t1pks_gene)
    
# Extracting information of domains from dictionary file 'dic_info_of_bundle_set_met' and 'dic_t1pks_gene' (VV)
# dic_t1pks_domain['SAV_942_DM6'] = ['PKS_AT', '(2425-2719)', 'SAV_942']
# dic_t1pks_gene_domain['SAV_938'] = ['PKS_AT', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_KR', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP', 'PKS_Docking_Cterm']
    dic_t1pks_domain, dic_t1pks_gene_domain, dic_t1pks_gene_substrate, dic_t1pks_PKS_KR_activity = second_metab_domain(dic_info_of_bundle_set_met, dic_t1pks_gene)
    
# Extracting information of substrates with their substrates from dictionary file 'dic_info_of_bundle_set_met' and 'dic_t1pks_gene' (VV)
# dic_t1pks_domain_substrate['SAV_943_M1'] = ['mmal', 'Ethyl_mal', 'pk']
    dic_t1pks_domain_substrate = second_metab_substrate(dic_info_of_bundle_set_met, dic_t1pks_gene)
   
# Extracting information of modules from dictionary file 'dic_t1pks_domain' (VV)
# dic_t1pks_module['SAV_943_M1'] = ['PKS_KS', 'PKS_AT', 'ACP']
    dic_t1pks_module = second_metab_module(dic_t1pks_gene, dic_t1pks_gene_domain, dic_t1pks_PKS_KR_activity)

# Generating rules for biosynthesis of type I PKS and converting module and its substrate to metabolic reactions
# For example : dic_converted_metabolic_reaction_without_substrate['SAV_943_M0'] = ['coa': 1, 'nadph': -1, 'nadp': 1, 'hco3': 1, 'h': -1]
    dic_converted_metabolic_reaction_without_substrate = generating_each_module_of_backbone_biosynthesis_for_t1pks(dic_t1pks_module)
    
# generating integrated metabolic reaction without participated substrate such as malonyl-coenzyme A
    dic_integrated_metabolic_reaction_without_cofactors = integrated_metabolic_reaction1(dic_converted_metabolic_reaction_without_substrate) #####
    
# adding matched participated substrate by using 'product' and dictionary 'dic_t1pks_domain_substrate'
# dic_semiintegrated_metabolic_reaction {'coa': 13, 'mmalcoa': -4, 'h': -10, 'malcoa': -7, 'hco3': 13, 'nadph': -10, 'h2o': 5, 'nadp': 10}
# list_of_dismatched_substrate = [['mmal', 'Ethyl_mal'], ['2metbut', '2metbut']]
    dic_semiintegrated_metabolic_reaction, list_of_dismatched_substrate = integrated_metabolic_reaction2(dic_t1pks_domain_substrate, dic_integrated_metabolic_reaction_without_cofactors)
    
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
     
# adding metabolic reactions (with cofactor) to the model
    modified_cobra_model1, list_novel_SMRs, list_reaction_name_SM = second_metab_reactions_addition(cobra_model, product, dic_t1pks_gene, list_of_reaction_set_with_product, 'output_participated_gene_list_of_backbone_biosynthesis.txt', 'output_database_format_file.txt')
    
# generating_control_reactions_against_newly_generated_biosynthetic_reactions_for_secondary_metabolites
    modified_product, list_of_completed_reaction_set_without_cofactor= generating_control_reactions_against_NovelSMR(list_reaction_name_SM, list_novel_SMRs, nadph_coef_number=0, amet_coef_number='None', fadh2_coef_number='None', nadh2_coef_number='None') 
    
# adding metabolic reactions (without cofactor) to the model
    modified_cobra_model2, list_novel_SMRs_without_cofactor, list_reaction_name_SM_without_cofactor = second_metab_reactions_addition(cobra_model = modified_cobra_model1, product = modified_product, dic_t1pks_gene = dic_t1pks_gene, list_of_reaction_set_with_product = list_of_completed_reaction_set_without_cofactor, outputfile1 = 'output_participated_gene_list_of_backbone_biosynthesis_x_cof.txt', outputfile2 = 'output_database_format_file_x_cof.txt')

# ## optional: removing uninterpretated reactiosn from model
#     reaction_list = ['AKGDH2']
#     modified_cobra_model2 = removing_any_metabolic_reaction_from_model(modified_cobra_model2, reaction_list) 
# 
# 
# # saving modified SBML model as xml form
#     sptcluster = inputfile.split('.')
#     name1 = sptcluster[1][1:]
#     name2 = sptcluster[3]
#     modified_model_name = name1 + '_' + name2 + '_' + 'model_optRMrevmoved.xml'
#     write_cobra_model_to_sbml_file(modified_cobra_model2, modified_model_name)
#     print "modified_model_is_generated"
#     
# # temporary file : writing textfile for newly generated reactions for secondary metabolites
#     outputfile1 = name1 + '_' + name2 + '_novel_SMRs.txt'
#     textfile1 = open(outputfile1,'w')
#     count = 0
#     for each_reaction1 in list_reaction_name_SM:
#         each_reaction2 = list_reaction_name_SM_without_cofactor[count]
#         print >>textfile1, "%s\t%s" % (each_reaction1, each_reaction2)
#         count += 1
#     print "generating_the_list_of_novel_SMRs_name"
        
# # creat_model
#     cobra_model2 = create_cobra_model_from_sbml_file(modified_model_name, print_time=True)
          
## simulating FBA 
    performing_FBA_for_each_reaction_of_SMRs(modified_cobra_model2, list_reaction_name_SM)

#     performing_FBA_for_each_reaction_of_SMRs(cobra_model2, list_reaction_name_SM)
#     performing_FBA_for_each_reaction_of_SMRs(cobra_model2, list_reaction_name_SM_without_cofactor)

## calculating phenotype_phase_plane (2D: version problem)
#     performing_2D_optimization_between_product_biomass(cobra_model = cobra_model2, biomassR_name = 'Biomass_SCO', list_novel_SMR_names = list_reaction_name_SM, start_number = 1, spt_max_number = 10)
#     performing_2D_optimization_between_product_biomass(cobra_model = cobra_model2, biomassR_name = 'Biomass_SCO', list_novel_SMR_names = list_reaction_name_SM_without_cofactor, start_number = 1, spt_max_number = 10)

## MOLP simulation tools
# generating reactions related to transport, export in order to 3D optimization
#     cobra_model3, Cofactor_transport_reaction, Cofactor_export_reaction = exports_reactions_addition(cobra_model = cobra_model2, list_export_met_name = ['nadph'], outputfile1 = './list_of_cofactor_transport_reaction.txt', outputfile2 = './list_of_cofactor_export_reaction.txt')
    
# simulating MOLP by 3D (product, cofactor, biomass)
#     performing_3D_optimization_between_product_cofactor_and_biomass(cobra_model = cobra_model3, biomassR_name = 'Biomass_SCO', cofactorR_name = Cofactor_export_reaction ,list_novel_SMR_names = list_reaction_name_SM, start_number = 1, spt_max_number = 10)
#     performing_3D_optimization_between_product_cofactor_and_biomass(cobra_model = cobra_model3, biomassR_name = 'Biomass_SCO', cofactorR_name = Cofactor_export_reaction ,list_novel_SMR_names = list_reaction_name_SM_without_cofactor, start_number = 1, spt_max_number = 10)

    
    


