'''
2015 Kyu-Sang Hwang
2014-2015 Hyun Uk Kim

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

def get_defined_sec_metab_monomers(inputFile):
    fp1 = open(inputFile,"r")
    monomer_mnx_dict = {}

    monomer = fp1.readline()

    while monomer:
        print monomer
        monomer = monomer.split("\t")
        monomer[0] = monomer[0].strip()
        monomer[1] = monomer[1].strip()
        monomer[2] = monomer[2].strip()
        monomer_mnx_dict[monomer[0]] = [monomer[1],monomer[2]]
        monomer = fp1.readline()

    print "\n", "List of secondary metabolic monomers:"
    print monomer_mnx_dict
    fp1.close()
    return monomer_mnx_dict


#Output: e.g., ['nrp', 'val-pro']
def get_monomers_from_cluster_gbk(gbkFile, FileType, monomer_mnx_dic):
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
            #e.g., ['pk-mmal-mal', 'mal-mal-mal-mmal', 'mmal-mal-mmal', 'mal-pk-mal']
            count = 0
            
            for each_module_substrate in second_total_monomers:
                each_module_substrate = each_module_substrate.replace('(','')
                each_module_substrate = each_module_substrate.replace(')','')
                each_module_substrate = each_module_substrate.strip()
                second_total_monomers[count] = each_module_substrate
                
                count += 1
            count = 1
    print second_total_monomers, "\n"
    return second_total_monomers

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
            #print "\n", "Cluster number:", clusterNo

            #Retrieving "product"
            product = feature.qualifiers.get('product')
            product = product[0]
            
            gene_strain = record.id
            gene_strain = gene_strain.split('.')
            gene_strain = gene_strain[0].strip()
            gene_strain = gene_strain.replace('_','') 
            #if gene_strain != None:
                #gene_strain = gene_strain[0].split(':')
                #gene_strain = gene_strain[1].split('.')
                #gene_strain = gene_strain[0].strip()
                #print "gene strain:", gene_strain
             #else:
                #gene_strain = 'unknown'
                #print "gene strain:", gene_strain
                  
    if float(clusterNo) < 10:
        product = gene_strain+"_"+"Cluster_0"+clusterNo+"_"+product
    else:
        product = gene_strain+"_"+"Cluster_"+clusterNo+"_"+product

    print product, "\n"
    return product


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
                #print qualifier_sec_et
                sec_met_info_list.append(qualifier_sec_met)

                cluster_info_dict[qualifier_locus_tag] = sec_met_info_list

    #print 'cluster_info_dict'
    #print cluster_info_dict, '\n'
    for i in cluster_info_dict.keys():
        print i
        print cluster_info_dict[i], '\n'
    return cluster_info_dict

def get_cluster_domain(cluster_info_dict):
    fp1 = open('Output_second_metab_gene_domain.txt','w')
    fp2 = open('Output_second_metab_gene_substrate.txt','w')
    #fp3 = open('Output_second_metab_gene_KR_activity.txt','w')
    
    dic_nrps_gene_domain = {}
    #dic_t1pks_PKS_KR_activity = {}
        
    for each_gene in cluster_info_dict.keys():
        
        list_set_met = cluster_info_dict[each_gene][0]

        domain_count = 0
        list_nrps_domain = []
        list_nrps_gene_substrate = []
        #list_t1pks_PKS_KR_activity = []
    
        for each_sub_set in list_set_met:

            if "NRPS/PKS Domain" in each_sub_set:
                            
                sptline1 = each_sub_set.split('. ')
                crude_domain_info = sptline1[0] 

                #Extract the information of KR activity and AT substrate specificity
                sptline2 = each_sub_set.split('; ')                
                
                if "AMP-binding" in crude_domain_info:
                    spt_NRPS_AMP = sptline2[1].split(':')
                    spt_NRPS_AMP_met = spt_NRPS_AMP[1].split(', ')
                    spt_NRPS_AMP_met = spt_NRPS_AMP_met[3].split(' (')
                    spt_NRPS_AMP_met = spt_NRPS_AMP_met[0].strip()
                    list_nrps_gene_substrate.append(spt_NRPS_AMP_met)
                
                #if "PKS_KR" in crude_domain_info:
                    #spt_PKS_KR = sptline2[1].split(': ')
                    #spt_PKS_KR = spt_PKS_KR[1].strip()
                    #list_t1pks_PKS_KR_activity.append(spt_PKS_KR)
                    
                spt_domain_info = crude_domain_info.split(':')
                whole_domain_info = spt_domain_info[1]

                spt_list_domain_info = whole_domain_info.split()
                spt_list_domain_info.append(each_gene)
                            
                each_nrps_domain = spt_list_domain_info[0]
                list_nrps_domain.append(each_nrps_domain)
                
                #Domain information
                domain_number = each_gene + '_DM' + str(domain_count)
                domain_count = domain_count + 1
                            
        dic_nrps_gene_domain[each_gene] = list_nrps_domain
        #dic_t1pks_PKS_KR_activity[each_gene] = list_t1pks_PKS_KR_activity
        
        #print each_gene, list_nrps_domain
        #print each_gene, list_nrps_gene_substrate
        #print each_gene, list_t1pks_PKS_KR_activity

        print >>fp1, "%s\t%s" % (each_gene, list_nrps_domain)
        print >>fp2, "%s\t%s" % (each_gene, list_nrps_gene_substrate), "\n"
        #print >>fp3, "%s\t%s" % (each_gene, list_t1pks_PKS_KR_activity)

    fp1.close()
    fp2.close()
    #fp3.close()

    print 'dic_nrps_gene_domain'
    print dic_nrps_gene_domain, '\n'
    return dic_nrps_gene_domain

def second_metab_substrate(cluster_info_dict):
    fp1 = open('Output_second_metab_substrate_with_domain.txt','w')

    dic_nrps_domain_substrate = {}
    for each_gene in cluster_info_dict.keys():
        module_count = 0
        list_set_met =  cluster_info_dict[each_gene][0]
        
        for each_sub_set in list_set_met:
     
            if "Substrate specificity predictions" in each_sub_set:
                sptline2 = each_sub_set.split(';')
                whole_substrate_info = sptline2[1]
                                
                participated_substrates = whole_substrate_info.split(':')
                sptSubstrates = participated_substrates[1]
                                
                if ', ' not in sptSubstrates:
                    print "insufficient substrate_info"
                    continue
                
                substrates = sptSubstrates.split(', ')

                list_participated_sustrate = []
                                
                for each_substrate in substrates:
                    sptSubstrate = each_substrate.split('(')
                    participated_substrate = sptSubstrate[0].strip()
                    list_participated_sustrate.append(participated_substrate)
                                
                module_number = each_gene + '_M' + str(module_count)
                dic_nrps_domain_substrate[module_number] = list_participated_sustrate
                
                #print module_number, list_participated_sustrate
                print >>fp1, "%s\t%s" % (module_number, list_participated_sustrate)
                    
                module_count = module_count + 1
    
    fp1.close()
    print 'dic_nrps_domain_substrate'
    print dic_nrps_domain_substrate, '\n'
    return dic_nrps_domain_substrate


def second_metab_module(locustag_product_monomer_dict, dic_nrps_gene_domain):
    fp1 = open('Output_second_metab_module.txt','w')
    
    dic_nrps_module = {}
    
    for nrps_gene in locustag_product_monomer_dict:
        
        count = 0
        if nrps_gene in dic_nrps_gene_domain:
            
            list_each_nrps_domain = dic_nrps_gene_domain[nrps_gene]
            #print list_each_nrps_domain
#             list_KR_activity = dic_t1pks_PKS_KR_activity[t1pks_gene]
        
            list_module_info = []
            number_of_list = len(list_each_nrps_domain)

#             KR_number = 0
            
            for each_domain in list_each_nrps_domain:    
#                 if each_domain == 'PKS_Docking_Nterm' or each_domain == 'PKS_Docking_Cterm':
#                     number_of_list = number_of_list - 1
#                     continue
                
                list_module_info.append(each_domain)
                number_of_list = number_of_list - 1
                
                if each_domain == 'PCP' or each_domain == 'ACP':
                    module_number = nrps_gene + '_M' + str(count)
                    dic_nrps_module[module_number] = list_module_info
                    print >>fp1, "%s\t%s\t%s" %(nrps_gene, module_number, list_module_info)
                    list_module_info = []
                    count = count + 1

                elif each_domain == 'Epimerization':
                    count = count - 1
                    module_number = nrps_gene + '_M' + str(count)
                    list_module_info = dic_nrps_module[module_number]
                    list_module_info.append('Epimerization')
                    A = dic_nrps_module.pop(module_number)
                    dic_nrps_module[module_number] = list_module_info
                    print >>fp1, "%s\t%s\t%s" %(nrps_gene, module_number, list_module_info)
                    list_module_info = []
                    count = count + 1

                elif each_domain == 'Thioesterase':
                    count = count - 1
                    module_number = nrps_gene + '_M' + str(count)
                    list_module_info = dic_nrps_module[module_number]
                    #print list_module_info
                    
                    list_module_info.append('Thioesterase')
                    #print list_module_info
                    
                    A = dic_nrps_module.pop(module_number)
                    #print A

                    dic_nrps_module[module_number] = list_module_info
                    print >>fp1, "%s\t%s\t%s" %(nrps_gene, module_number, list_module_info)
                    
                elif list_module_info.count('Condensation_DCL') == 2 or list_module_info.count('Condensation_LCL') == 2 or list_module_info.count('Condensation_LCL') + list_module_info.count('Condensation_DCL') == 2:
                    module_number = nrps_gene + '_M' + str(count)
                    list_module_info.pop()
                    dic_nrps_module[module_number] = list_module_info
                    list_module_info = []
                    count = count + 1
                    
                elif float(number_of_list) == 0:
                    module_number = nrps_gene + '_M' + str(count)
                    dic_nrps_module[module_number] = list_module_info
                    print >>fp1, "%s\t%s\t%s" % (nrps_gene, module_number, list_module_info)
                    list_module_info = []
                    count = count + 1
                
    fp1.close()
    print 'dic_nrps_module'
    print dic_nrps_module, '\n'
    return dic_nrps_module


def generating_each_module_of_backbone_biosynthesis_for_t1pks(dic_nrps_module):
    fp1 = open('Output_backbone_biosynthesis_of_each_module.txt','w')
    
    dic_converted_metabolic_reaction_without_substrate = {}

    for each_module in dic_nrps_module:
        domain_comb = dic_nrps_module[each_module]
        
        print domain_comb

        each_module_substrates = {}
        discriminant = module_discriminator(domain_comb)
        
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
            
        elif discriminant == 'TD':
            each_module_substrates['nadp'] = 1 #'NADP+', (C00006), 'MNXM5'
            each_module_substrates['nadph'] = -1 #'NADPH' ,(C00005), 'MNXM6'
            each_module_substrates['h'] = -1 #'H+', (C00080), 'MNXM1'
            dic_converted_metabolic_reaction_without_substrate[each_module] = each_module_substrates                
            print 'reaction 8: HC_Aox_MT_PCP', each_module_substrates 
            print >>fp1, "module_type:\t%s\t%s\t%s" % (each_module, domain_comb, each_module_substrates)
            
        elif discriminant == 'PCP':
            continue
    
    fp1.close()
                
    return dic_converted_metabolic_reaction_without_substrate

def module_discriminator(domain_comb):
## domain information :
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

# Exceptionsal cases
    if ('Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb) and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'A'
    
# Starter unit_AMP-binding+PCP
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'A_PCP' #
        
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'MT_A_PCP' ##
        
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'MT_A_PCP'
              

# Condensation starter domain : C-domain
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_PCP' #   
       
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_PCP' ##
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_PCP' ##
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_PCP_E' ####
        
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_E_PCP' ###

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_E_PCP'
        
# Condensation domain that links D-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_PCP' #
         
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_PCP' ##
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_PCP' ##
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_PCP_E' ####
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_E_PCP' ###
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_E_PCP' ###        

# Condensation domain that links D-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_PCP' #
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_PCP_E' ####
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_E_PCP' ###
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_E_PCP' ###


# Condensation domain that links L-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_PCP' #
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_PCP' ##
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_PCP' ##
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_PCP_E' ####
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_E_PCP' ###
        
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_E_PCP' ###


# Condensation_dual : C-domain
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cd_A_PCP' #### + epimerization
        
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cd_A_MT_PCP' ###
        
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cd_A_MT_PCP' ###
        

# Glycopeptide condensation domain (O) : C-domain
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_PCP' #
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_PCP' ##
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_PCP' ##
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_PCP_E' ####
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_E_PCP' ###
        
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_E_PCP' ###

# Glycopeptide condensation domain (X) : C-domain
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_PCP' #
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_PCP' ##
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_PCP' ##
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_PCP_E' ####
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_E_PCP' ###
        
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_E_PCP' ###

# Heterocyclization : C-domain
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_E_PCP' ###
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_E_PCP'  
        
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_E_PCP'  

# Terminal reductase domain : terminal domain       
    elif 'TD' in domain_comb and ('Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb) and 'AMP-binding' not in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'TD'
    
    else:
        discriminant = 'None'
        
    return discriminant

# def Identifier_KR_activity(discriminant, each_module_KR_activity):
#     
#     print discriminant, each_module_KR_activity
# 
#     if each_module_KR_activity == 1 or each_module_KR_activity == 0:
#         discriminant_with_KRact = discriminant
#         return discriminant_with_KRact
#     
#     elif each_module_KR_activity == 2:
#         discriminant_with_KRact = discriminant.replace('KR','KR(inactive)')      
#         return discriminant_with_KRact

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

def integrated_metabolic_reaction2(dic_nrps_domain_substrate, dic_integrated_metabolic_reaction_without_cofactors):
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
    dic_integrated_metabolic_reaction_without_cofactors['hpg'] = 0 #'D-4-Hydroxyphenylglycine', 'C03493', 'MNXM4544'
    dic_integrated_metabolic_reaction_without_cofactors['23dhb'] = 0 #'2,3-Dihydroxybenzoic acid', 'C00196', 'MNXM455'
    dic_integrated_metabolic_reaction_without_cofactors['dhpg'] = 0 #'3,5-Dihydroxy-phenylglycine', 'C12026', 'MNXM9962'
    dic_integrated_metabolic_reaction_without_cofactors['hty'] = 0 #'L-Homotyrosine', 'C18622', 'MNXM59438'
    dic_integrated_metabolic_reaction_without_cofactors['citr_DASH_L'] = 0 #'L-citruline', 'C00327', 'MNXM211' 
    dic_integrated_metabolic_reaction_without_cofactors['Lpipecol'] = 0 #'L-pipecolate', 'C00408', 'MNXM684'
    dic_integrated_metabolic_reaction_without_cofactors['24dab'] = 0 #'L-2,4-diazaniumylbutyrate', 'C03283', 'MNXM840'
    dic_integrated_metabolic_reaction_without_cofactors['ala_DASH_B'] = 0 #'beta-alanine zwitterion', 'C00099', 'MNXM144'
    dic_integrated_metabolic_reaction_without_cofactors['tcl'] = 0 #'4-Chlorothreonine', 'N/A', 'MNXM37380'
    dic_integrated_metabolic_reaction_without_cofactors['qa'] = 0 #'quinoxaline', 'C18575','MNXM80501' VV
    
    for each_module in dic_nrps_domain_substrate:

        if dic_nrps_domain_substrate[each_module][3] == 'nrp':
            
            sptlist1 = dic_nrps_domain_substrate[each_module][0].split(',')
            temp_list = []
            
            if len(sptlist1) >= 2:
                for each_substrate in sptlist1:
                    converted_met1 = converting_substrate_to_met_coeff(each_substrate)
                    temp_list.append(converted_met1)
            elif dic_nrps_domain_substrate[each_module][0] == 'hydrophobic-aliphatic' or dic_nrps_domain_substrate[each_module][0] == 'hydrophilic':
                each_substrate = 'N/A'
                temp_list.append(each_substrate)    
            else:
                query_met2 = dic_nrps_domain_substrate[each_module][0]
                converted_met2 = converting_substrate_to_met_coeff(query_met2)
                temp_list.append(converted_met2)
                
            query_met3 = dic_nrps_domain_substrate[each_module][1]
            converted_met3 = converting_substrate_to_met_coeff(query_met3)
            temp_list.append(converted_met3)
            
            query_met4 = dic_nrps_domain_substrate[each_module][2]
            converted_met4 = converting_substrate_to_met_coeff(query_met4)
            temp_list.append(converted_met4)
            
            list_of_dismatched_substrate.append(temp_list)
        
        elif dic_nrps_domain_substrate[each_module][3] != 'nrp':
            query_met5 = dic_nrps_domain_substrate[each_module][3]
            converted_met5 = converting_substrate_to_met_coeff(query_met5)
            dic_integrated_metabolic_reaction_without_cofactors[converted_met5] -= 1
             
    for each_participated_substrate in dic_integrated_metabolic_reaction_without_cofactors:
        print >>fp1, "total_participated_cofactor_with_substrates1:\t%s\t%s" % (each_participated_substrate, dic_integrated_metabolic_reaction_without_cofactors[each_participated_substrate])
    
    fp1.close()
    
    return dic_integrated_metabolic_reaction_without_cofactors, list_of_dismatched_substrate

def converting_substrate_to_met_coeff(each_substrate):
    
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
        
    elif each_substrate == 'tcl':
        met_name = 'tcl'
        
    elif each_substrate == 'b-ala':
        met_name = 'ala_DASH_B'
        
    elif each_substrate == 'phenylacetate' or each_substrate == 'Pha':
        met_name = 'pac'
        
    elif each_substrate == 'qa':
        met_name = 'qa'
        
    elif each_substrate == 'pip':
        met_name = 'Lpipecol'
        
    elif each_substrate == 'dab':
        met_name = '24dab'
        
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
            template_list = copy.deepcopy(list_of_reaction_set)
            list_of_reaction_set = []
        
            for dic_each_metabolic_reaction in template_list:
                each_pair_of_substrates = set(each_pair_of_substrates)
                each_pair_of_substrates = list(each_pair_of_substrates)
                substrate_decision_number = distincting_each_substrate_in_list_component(each_pair_of_substrates)
                temp_reaction_set = converting_nrps_substrates(each_pair_of_substrates, dic_each_metabolic_reaction, substrate_decision_number) 
                
                for each_dic_reaction_set in temp_reaction_set:
                    list_of_reaction_set.append(each_dic_reaction_set)
    
    return list_of_reaction_set
                
def distincting_each_substrate_in_list_component(each_pair_of_substrates):    
    
    # this logic of code should be fixed
    if each_pair_of_substrates[0] != each_pair_of_substrates[1] or each_pair_of_substrates[1] != each_pair_of_substrates[2]:
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
            
        elif each_substrate == 'tcl':
             met_name = 'tcl'
        
        elif each_substrate == 'ala_DASH_B':
            met_name = 'ala_DASH_B'
        
        elif each_substrate == 'pac':
            met_name = 'pac'
        
        elif each_substrate == 'qa':
            met_name = 'qa'    
            
        elif each_substrate == 'Lpipecol':
            met_name = 'Lpipecol'   
            
        elif each_substrate == '24dab':
            met_name = '24dab'                       
         
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
    elif MetID == 'tcl':
        converted_MNXMID = ['MNXM37380', 'tcl', '4-Chlorothreonine', 'N/A'] 
    elif MetID == 'ala_DASH_B':
        converted_MNXMID = ['MNXM144', 'ala_DASH_B', 'beta-Alanine', 'C00099']  
    elif MetID == 'pac':
        converted_MNXMID = ['NXM497', 'pac', 'phenylacetate', 'C00548']  
    elif MetID == 'qa':
        converted_MNXMID = ['MNXM4797', 'qa', 'Quinaldinic acid', 'C06325'] 
    elif MetID == 'Lpipecol':
        converted_MNXMID = ['MNXM684', 'Lpipecol', 'L-pipecolate', 'C00408']       
    elif MetID == '24dab':
        converted_MNXMID = ['MNXM840', '24dab', 'L-2,4-diazaniumylbutyrate', 'C03283']
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

def second_metab_reactions_addition(cobra_model, product, locustag_product_monomer_dict, list_of_reaction_set_with_product, metab_MNXM_dict):
#     fp1 = open('output_set_of_integrated_metabolic_reaction.txt','w')
    fp2 = open('output_participated_gene_list_of_backbone_biosynthesis.txt','w')
    fp3 = open("output_database_format_file.txt", 'w')
    product_count = 1
    
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
        for each_gene in locustag_product_monomer_dict:
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
#     number = 1
    
    for each_MR in list_novel_secondary_metabolite_reaction_names:
       
        output = each_MR + '_max.txt'
        query_reaction = 'Transport_'+each_MR
#         print number
        calculated_cobra_model = cobra_model_FBA(cobra_model, query_reaction, 'maximize', output)
#         number += 1
#     fp1.close()
    
    return 

    
# running FBA
def cobra_model_FBA(cobra_model, objective, sense, output):

#     fp1 = open(output,"w")
    reaction = cobra_model.reactions.get_by_id(objective)

    # Run the optimization for the objective reaction and medium composition
    cobra_model.optimize(new_objective=reaction, objective_sense=sense, solver='gurobi')
    print cobra_model.solution.f
#     fp1.write(str(cobra_model.solution.status)+"\n")
#     fp1.write(str(cobra_model.solution.f)+"\n"+"\n")

#     for the_reaction, the_value in cobra_model.solution.x_dict.items():
#         #print '%s: %1.2f' %(the_reaction, the_value)
#         fp1.write(str(the_reaction)+"\t"+str(the_value)+"\n")

#     fp1.close()

    return cobra_model
       
