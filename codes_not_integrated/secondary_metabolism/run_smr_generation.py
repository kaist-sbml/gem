'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
import pickle
import os
from sec_met_rxn_generation import get_product_from_cluster_gbk, get_cluster_info_from_cluster_gbk, get_cluster_domain, get_cluster_monomers, get_cluster_module, get_currency_metabolites, get_total_currency_metab_coeff, get_all_metab_coeff, add_sec_met_rxn 

print "Generating secondary metabolite biosynthesizing reactions.."

target_model = create_cobra_model_from_sbml_file('sma_target_model_sco.xml')
#inputfile = './NC_021055.1.cluster002.gbk' #NRPS
#inputfile = './NC_013929.1.cluster031.gbk' #PKS
#inputfile = './NC_020990.1.cluster023.gbk' #Hybrid

bigg_mnxm_compound_dict = pickle.load(open('bigg_mnxm_compound_dict.p','rb'))
mnxm_compoundInfo_dict = pickle.load(open('mnxm_compoundInfo_dict.p','rb'))

dirname = './target_genome/'
#dirname = './'
cluster_files = []
for inputfile in os.listdir(dirname):
    if inputfile.endswith('.gbk') and 'cluster' in inputfile:
        cluster_files.append(inputfile)

cluster_files.sort()

for cluster_f in cluster_files:
    print '\n', cluster_f

    cluster_info_dict = get_cluster_info_from_cluster_gbk(dirname+cluster_f, "genbank")

    product = get_product_from_cluster_gbk(dirname+cluster_f, "genbank")

    if 't1pks' in product or 'nrps' in product:

        locustag_domain_dict = get_cluster_domain(cluster_info_dict)

        locustag_monomer_dict = get_cluster_monomers(cluster_info_dict)

        locustag_module_domain_dict = get_cluster_module(locustag_domain_dict)

        module_currency_metab_dict = get_currency_metabolites(locustag_module_domain_dict)

        currency_metab_coeff_dict = get_total_currency_metab_coeff(module_currency_metab_dict)

        metab_coeff_dict = get_all_metab_coeff(locustag_monomer_dict, currency_metab_coeff_dict, product)

        target_model = add_sec_met_rxn(target_model, metab_coeff_dict, product, bigg_mnxm_compound_dict, mnxm_compoundInfo_dict, cluster_info_dict)

write_cobra_model_to_sbml_file(target_model, 'sma_target_model_sco_complete.xml')
