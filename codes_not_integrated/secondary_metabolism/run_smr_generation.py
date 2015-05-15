'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
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
from secondary_metabolite_biosynthesis_nrps import *

print "Generating NRP biosynthesis reactions.."

# making template model in order to (V)
cobra_model = create_cobra_model_from_sbml_file('SCO_model_snu.xml', print_time=True)
inputfile = './NC_021055.1.cluster002.gbk'

# making the dictionary file of metabolites in template model (V)
MNXM_dict = pickle.load(open("Pickle_MNXM_dict.p","rb"))

# making dictionary file of metabolites in template model (V)
# (e.g. ''MNXM37': 'gln_DASH_L_p'])
metab_MNXM_dict = COBRA_TemplateModel_checking_MNXref_metabolites(cobra_model, MNXM_dict)
print metab_MNXM_dict

# Enrolling participated substrates for backbone biosynthesis in Type I polyketide (V)
#monomer_mnx_dic = monomers_list_for_nrps('Input_monomers_nrps.txt')

# extracting information of total combined (V)
#second_total_monomers = second_metab_monomers(inputfile, "genbank", monomer_mnx_dic)

# To generate the product name of target 'type I pks' genes (V)
#product = second_metab_reaction_product_names(inputfile, "genbank")

# To extract genes related to backbone biosynthesis in type I polyketide synthase (VV)
# Extracting information of locus_tag, locus_tag, substrate_inform(aSProdPred) from query genebankfile
# dic_t1pks_gene['SAV_938'] = ['type I polyketide synthase AVES 1', 'pk-mmal-mal']
#dic_nrps_gene = second_metab_genes(inputfile, "genbank")

# Extracting information of domains from dictionary file 'dic_t1pks_domain' (VV)
#  dic_t1pks_domain[SAV_942_DM12] = ['PKS_KR', '(5084-5264)']
#dic_info_of_bundle_set_met = extracting_sub_set_met_info_from_genebank(inputfile, "genbank", dic_nrps_gene)

# Extracting information of domains from dictionary file 'dic_info_of_bundle_set_met' and 'dic_t1pks_gene' (VV)
# dic_t1pks_domain['SAV_942_DM6'] = ['PKS_AT', '(2425-2719)', 'SAV_942']
# dic_t1pks_gene_domain['SAV_938'] = ['PKS_AT', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_KR', 'ACP', 'PKS_KS', 'PKS_AT', 'PKS_DH', 'PKS_KR', 'ACP', 'PKS_Docking_Cterm']
#dic_nrps_domain, dic_nrps_gene_domain, dic_nrps_gene_substrate = second_metab_domain(dic_info_of_bundle_set_met, dic_nrps_gene)

# Extracting information of substrates with their substrates from dictionary file 'dic_info_of_bundle_set_met' and 'dic_t1pks_gene' (VV)
# dic_t1pks_domain_substrate['SAV_943_M1'] = ['mmal', 'Ethyl_mal', 'pk']
#dic_nrps_domain_substrate = second_metab_substrate(dic_info_of_bundle_set_met, dic_nrps_gene)

# Extracting information of modules from dictionary file 'dic_t1pks_domain' (VV)
# dic_t1pks_module['SAV_943_M1'] = ['PKS_KS', 'PKS_AT', 'ACP']
#dic_nrps_module = second_metab_module(dic_nrps_gene, dic_nrps_gene_domain)

# Generating rules for biosynthesis of type I PKS and converting module and its substrate to metabolic reactions
# For example : dic_converted_metabolic_reaction_without_substrate['SAV_943_M0'] = ['coa': 1, 'nadph': -1, 'nadp': 1, 'hco3': 1, 'h': -1]
#dic_converted_metabolic_reaction_without_substrate = generating_each_module_of_backbone_biosynthesis_for_t1pks(dic_nrps_module)

# generating integrated metabolic reaction without participated substrate such as malonyl-coenzyme A
#dic_integrated_metabolic_reaction_without_cofactors = integrated_metabolic_reaction1(dic_converted_metabolic_reaction_without_substrate) #####

# adding matched participated substrate by using 'product' and dictionary 'dic_t1pks_domain_substrate'
# dic_semiintegrated_metabolic_reaction {'coa': 13, 'mmalcoa': -4, 'h': -10, 'malcoa': -7, 'hco3': 13, 'nadph': -10, 'h2o': 5, 'nadp': 10}
# list_of_dismatched_substrate = [['mmal', 'Ethyl_mal'], ['2metbut', '2metbut']]
#dic_semiintegrated_metabolic_reaction, list_of_dismatched_substrate = integrated_metabolic_reaction2(dic_nrps_domain_substrate, dic_integrated_metabolic_reaction_without_cofactors)

# completing integrated metabolic reaction by adding product and dismatched substrate to the reaction.
# list_of_reaction_set = [{'nadph': -10, 'nadp': 10, 'ahcys': 0, '2mbcoa': -1, 'nad': 0, 'h': -10, 'fadh2': 0, 'malcoa': -7, 'hco3': 13, 'amet': 0, 'coa': 13, 'h2o': 5, 'nadh': 0, '13dpg': 0, 'mmalcoa': -5, 'pi': 0, 'emalcoa': 0, 'fad': 0}, ...]
#list_of_reaction_set = integrated_metabolic_reaction3(dic_semiintegrated_metabolic_reaction, list_of_dismatched_substrate)

# adding product name to the integrated reaction
# list_of_reaction_set_with_product = {'nadph': -10, 'nadp': 10, 'ahcys': 0, '2mbcoa': -1, 'nad': 0, 'h': -10, 'fadh2': 0, 'malcoa': -7, 'hco3': 13, 'Cluster_05_t1pks_1': 1, 'amet': 0, 'coa': 13, 'h2o': 5, 'nadh': 0, '13dpg': 0, 'mmalcoa': -5, 'pi': 0, 'emalcoa': 0, 'fad': 0}]
#list_of_reaction_set_with_product = adding_product_to_the_reaction(product, list_of_reaction_set)

# making all metabolic reactions from each putative type 1 PKS gene clusters.
# dic_reaction_info_set[gene_info] = {reaction_set}
#     dic_converted_metabolic_reactions_with_gene = making_metabolic_reactions_from_gene_cluster_data(dic_t1pks_gene, list_of_reaction_set_with_product)

# making total list of metabolites in model (e.g. metab_MNXM_dict[MNX_ID(e.g. MNXM37)] = metabolite_name_in_model(e.g. gln_DASH_L_p))
#     metab_MNXM_dict = COBRA_TemplateModel_checking_MNXref_metabolites(cobra_model, MNXM_dict)

# adding metabolic reactions to the model
#modified_cobra_model, list_reaction_name_SM, list_novel_secondary_metabolite_reactions = second_metab_reactions_addition(cobra_model, product, dic_nrps_gene, list_of_reaction_set_with_product, metab_MNXM_dict)

## simulating FBA
#performing_FBA_for_each_reaction_of_SMRs(modified_cobra_model, list_reaction_name_SM)

