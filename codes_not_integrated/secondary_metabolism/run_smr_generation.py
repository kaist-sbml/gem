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
from sec_met_rxn_generation import get_product_from_cluster_gbk, get_cluster_info_from_cluster_gbk, get_cluster_domain, get_cluster_monomers, get_cluster_module, get_currency_metabolites, get_total_currency_metab_coeff, get_all_metab_coeff, add_sec_met_rxn 

print "Generating secondary metabolite biosynthesizing reactions.."

#WHY DO WE NEED THIS?
cobra_model = create_cobra_model_from_sbml_file('sma_target_model_sco.xml')
#inputfile = './NC_021055.1.cluster002.gbk' #NRPS
inputfile = './NC_013929.1.cluster031.gbk' #PKS
#inputfile = './NC_020990.1.cluster023.gbk' #Hybrid

#TO MODIFY LATER
mnxm_bigg_compound_dict = pickle.load(open('bigg_mnxm_compound_dict.p','rb'))
mnxm_bigg_compound_dict = pickle.load(open('mnxm_compoundInfo_dict.p','rb'))

#WHO DO WE NEED THIS?
# making dictionary file of metabolites in template model (V)
# (e.g. ''MNXM37': 'gln_DASH_L_p'])
#metab_MNXM_dic = COBRA_TemplateModel_checking_MNXref_metabolites(cobra_model, mnxm_bigg_compound_dict)
##########################################################################


cluster_info_dict = get_cluster_info_from_cluster_gbk(inputfile, "genbank")

#MIGHT BE UNNECESSARY
product = get_product_from_cluster_gbk(inputfile, "genbank")

#locustag_product_monomer_dict: MIGHT NOT BE NECESSARY
locustag_domain_dict = get_cluster_domain(cluster_info_dict)

#MIGHT BE REMOVED
locustag_monomer_dict = get_cluster_monomers(cluster_info_dict)

locustag_module_domain_dict = get_cluster_module(locustag_domain_dict)

module_currency_metab_dict = get_currency_metabolites(locustag_module_domain_dict)

currency_metab_coeff_dict = get_total_currency_metab_coeff(module_currency_metab_dict)

metab_coeff_dict = get_all_metab_coeff(locustag_monomer_dict, currency_metab_coeff_dict, product)

#IN PROGRESS
#Metabolit parts to check
add_sec_met_rxn(cobra_model, metab_coeff_dict, product, bigg_mnxm_compound_dict, mnxm_compoundInfo_dict, cluster_info_dict)
