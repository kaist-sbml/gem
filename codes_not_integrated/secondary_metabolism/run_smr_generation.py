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
from sec_met_rxn_generation import *

print "Generating NRP biosynthesis reactions.."

#WHY DO WE NEED THIS?
# making template model in order to (V)
#cobra_model = create_cobra_model_from_sbml_file('SCO_model_snu.xml', print_time=True)
#inputfile = './NC_021055.1.cluster002.gbk' #NRPS
#inputfile = './NC_013929.1.cluster031.gbk' #PKS
inputfile = './NC_020990.1.cluster023.gbk' #Hybrid

mnxm_bigg_compound_dict = pickle.load(open('mnxm_bigg_compound_dict.p','rb'))

#WHO DO WE NEED THIS?
# making dictionary file of metabolites in template model (V)
# (e.g. ''MNXM37': 'gln_DASH_L_p'])
#metab_MNXM_dic = COBRA_TemplateModel_checking_MNXref_metabolites(cobra_model, mnxm_bigg_compound_dict)
#print metab_MNXM_dict

#TO CREATE A INPUT PICKLE FILE; NO NEED TO REPEAT EVERY TIME
monomer_mnx_dict = get_defined_sec_metab_monomers('Input_monomers_nrps.txt')

#MIGHT BE REMOVED
#second_total_monomers = get_monomers_from_cluster_gbk(inputfile, "genbank", monomer_mnx_dict)

#MIGHT BE UNNECESSARY
product = get_product_from_cluster_gbk(inputfile, "genbank")

cluster_info_dict = get_cluster_info_from_cluster_gbk(inputfile, "genbank")

#locustag_product_monomer_dict: MIGHT NOT BE NECESSARY
locustag_domain_dict = get_cluster_domain(cluster_info_dict)

#locustag_monomer_dict = get_cluster_monomers(cluster_info_dict)

locustag_module_domain_dict = get_cluster_module(locustag_domain_dict)

# Generating rules for biosynthesis of type I PKS and converting module and its substrate to metabolic reactions
# For example : dic_converted_metabolic_reaction_without_substrate['SAV_943_M0'] = ['coa': 1, 'nadph': -1, 'nadp': 1, 'hco3': 1, 'h': -1]
dic_converted_metabolic_reaction_without_substrate = generate_currency_metabolites(locustag_module_domain_dict)

# generating integrated metabolic reaction without participated substrate such as malonyl-coenzyme A
#dic_integrated_metabolic_reaction_without_cofactors = integrated_metabolic_reaction1(dic_converted_metabolic_reaction_without_substrate) #####

# adding matched participated substrate by using 'product' and dictionary 'dic_t1pks_domain_substrate'
# dic_semiintegrated_metabolic_reaction {'coa': 13, 'mmalcoa': -4, 'h': -10, 'malcoa': -7, 'hco3': 13, 'nadph': -10, 'h2o': 5, 'nadp': 10}
# list_of_dismatched_substrate = [['mmal', 'Ethyl_mal'], ['2metbut', '2metbut']]
#dic_semiintegrated_metabolic_reaction, list_of_dismatched_substrate = integrated_metabolic_reaction2(locustag_monomer_dict, dic_integrated_metabolic_reaction_without_cofactors)

# completing integrated metabolic reaction by adding product and dismatched substrate to the reaction.
# list_of_reaction_set = [{'nadph': -10, 'nadp': 10, 'ahcys': 0, '2mbcoa': -1, 'nad': 0, 'h': -10, 'fadh2': 0, 'malcoa': -7, 'hco3': 13, 'amet': 0, 'coa': 13, 'h2o': 5, 'nadh': 0, '13dpg': 0, 'mmalcoa': -5, 'pi': 0, 'emalcoa': 0, 'fad': 0}, ...]
#list_of_reaction_set = integrated_metabolic_reaction3(dic_semiintegrated_metabolic_reaction, list_of_dismatched_substrate)

# adding product name to the integrated reaction
# list_of_reaction_set_with_product = {'nadph': -10, 'nadp': 10, 'ahcys': 0, '2mbcoa': -1, 'nad': 0, 'h': -10, 'fadh2': 0, 'malcoa': -7, 'hco3': 13, 'Cluster_05_t1pks_1': 1, 'amet': 0, 'coa': 13, 'h2o': 5, 'nadh': 0, '13dpg': 0, 'mmalcoa': -5, 'pi': 0, 'emalcoa': 0, 'fad': 0}]
#list_of_reaction_set_with_product = adding_product_to_the_reaction(product, list_of_reaction_set)

# adding metabolic reactions to the model
#modified_cobra_model, list_reaction_name_SM, list_novel_secondary_metabolite_reactions = second_metab_reactions_addition(cobra_model, product, locustag_product_monomer_dict, list_of_reaction_set_with_product, metab_MNXM_dict)

## simulating FBA
#performing_FBA_for_each_reaction_of_SMRs(modified_cobra_model, list_reaction_name_SM)

