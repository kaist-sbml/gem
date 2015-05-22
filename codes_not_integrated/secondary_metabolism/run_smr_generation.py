'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

from Bio import SeqIO
from sets import Set
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from MNX_checker2 import COBRA_TemplateModel_checking_MNXref_metabolites, fix_legacy_id
import pickle
import copy
from sec_met_rxn_generation import get_defined_sec_metab_monomers, get_product_from_cluster_gbk, get_cluster_info_from_cluster_gbk, get_cluster_domain, get_cluster_monomers, get_cluster_module, get_currency_metabolites, get_total_currency_metab_coeff, get_all_metab_coeff, integrated_metabolic_reaction3, adding_product_to_the_reaction 

print "Generating NRP biosynthesis reactions.."

#WHY DO WE NEED THIS?
# making template model in order to (V)
#cobra_model = create_cobra_model_from_sbml_file('SCO_model_snu.xml', print_time=True)
inputfile = './NC_021055.1.cluster002.gbk' #NRPS
#inputfile = './NC_013929.1.cluster031.gbk' #PKS
#inputfile = './NC_020990.1.cluster023.gbk' #Hybrid

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

#MIGHT BE REMOVED
locustag_monomer_dict = get_cluster_monomers(cluster_info_dict)

locustag_module_domain_dict = get_cluster_module(locustag_domain_dict)

module_currency_metab_dict = get_currency_metabolites(locustag_module_domain_dict)

currency_metab_coeff_dict = get_total_currency_metab_coeff(module_currency_metab_dict)

#Working on this
metab_coeff_dict, dismatched_substrate_list = get_all_metab_coeff(locustag_monomer_dict, currency_metab_coeff_dict)

#Working on this
#Current version of this creates different reactions for non-consensus monomers: need to be modified
#Two nested functions
list_of_reaction_set = integrated_metabolic_reaction3(metab_coeff_dict, dismatched_substrate_list)

# adding product name to the integrated reaction
# list_of_reaction_set_with_product = {'nadph': -10, 'nadp': 10, 'ahcys': 0, '2mbcoa': -1, 'nad': 0, 'h': -10, 'fadh2': 0, 'malcoa': -7, 'hco3': 13, 'Cluster_05_t1pks_1': 1, 'amet': 0, 'coa': 13, 'h2o': 5, 'nadh': 0, '13dpg': 0, 'mmalcoa': -5, 'pi': 0, 'emalcoa': 0, 'fad': 0}]
#list_of_reaction_set_with_product = adding_product_to_the_reaction(product, list_of_reaction_set)

# adding metabolic reactions to the model
#modified_cobra_model, list_reaction_name_SM, list_novel_secondary_metabolite_reactions = second_metab_reactions_addition(cobra_model, product, locustag_product_monomer_dict, list_of_reaction_set_with_product, metab_MNXM_dict)

## simulating FBA
#performing_FBA_for_each_reaction_of_SMRs(modified_cobra_model, list_reaction_name_SM)

