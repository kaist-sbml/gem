'''
2015 Hyun Uk Kim
2015 Kyu-Sang Hwang
'''

import logging
import os
import pickle
import sys
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from cobra.manipulation.delete import prune_unused_metabolites
from sec_met_rxn_generation import (
    get_product_from_cluster_gbk,
    get_cluster_info_from_cluster_gbk,
    get_cluster_domain, get_cluster_monomers,
    get_cluster_module, get_currency_metabolites,
    get_total_currency_metab_coeff,
    get_all_metab_coeff,
    add_sec_met_rxn,
    check_producibility_sec_met,
    get_monomers_nonprod_sec_met,
    get_monomers_prod_sec_met,
)
from gapfill_network_manipulation import (
    get_mnxr_bigg_in_target_model,
    get_mnxr_unique_to_universal_model,
    integrate_target_universal_models,
    add_transport_exchange_rxn_nonprod_monomer,
    check_producibility_nonprod_monomer,
    get_unique_nonprod_monomers_list,
    execute_gapfill,
    check_gapfill_rxn_biomass_effects,
    add_gapfill_rxn_target_model,
)

logging.debug("Generating secondary metabolite biosynthesizing reactions..")

orgname = sys.argv[1]

bigg_mnxm_compound_dict = pickle.load(open('./input2/bigg_mnxm_compound_dict.p','rb'))
mnxm_compoundInfo_dict = pickle.load(open('./input2/mnxm_compoundInfo_dict.p','rb'))
template_exrxnid_flux_dict = pickle.load(open('./input1/%s/tempModel_exrxnid_flux_dict.p' %'sco','rb'))

if '/' in orgname:
    orgname = orgname[:-1]

dirname = './%s/' %orgname
    
cluster_files = []
for f in os.listdir(dirname):
    if f.endswith('.gbk') and 'cluster' in f:
        cluster_files.append(f)

cluster_files.sort()

for f in os.listdir(dirname+'2_primary_metabolic_model'):
    if f.endswith('.xml'):
       model_sbml = f

target_model = create_cobra_model_from_sbml_file(dirname+'2_primary_metabolic_model/'+model_sbml)

#if __name__ == '__main__':
#    cluster_f = 'NC_018750.1.cluster003.gbk'

prod_sec_met_dict = {}
nonprod_sec_met_dict = {}

for cluster_f in cluster_files:
    logging.debug(cluster_f)

    cluster_info_dict, record = get_cluster_info_from_cluster_gbk(dirname+cluster_f, "genbank")

    product = get_product_from_cluster_gbk(record)

    if 't1pks' in product or 'nrps' in product:

        locustag_domain_dict, locustag_kr_dict = get_cluster_domain(cluster_info_dict)

        locustag_monomer_dict = get_cluster_monomers(cluster_info_dict)

        locustag_module_domain_dict = get_cluster_module(locustag_domain_dict)

        module_currency_metab_dict = get_currency_metabolites(locustag_module_domain_dict, locustag_kr_dict)

        currency_metab_coeff_dict = get_total_currency_metab_coeff(module_currency_metab_dict)

        metab_coeff_dict = get_all_metab_coeff(locustag_monomer_dict, currency_metab_coeff_dict, product)

        target_model = add_sec_met_rxn(target_model, metab_coeff_dict, product, bigg_mnxm_compound_dict, mnxm_compoundInfo_dict, cluster_info_dict)
        
        target_model, product = check_producibility_sec_met(target_model, product, dirname)

        if target_model.solution.f < 0.0001:
            nonprod_sec_met_metab_list = get_monomers_nonprod_sec_met(metab_coeff_dict)
            nonprod_sec_met_dict[product] = nonprod_sec_met_metab_list
        else:
            prod_sec_met_metab_list = get_monomers_prod_sec_met(metab_coeff_dict)
            prod_sec_met_dict[product] = prod_sec_met_metab_list

logging.debug("Producible secondary metabolites:")
logging.debug(prod_sec_met_dict)

logging.debug("Nonproducible secondary metabolites:")
logging.debug(nonprod_sec_met_dict)

logging.debug("Gap-filling for the production of secondary metabolites..")
logging.debug("Step 1: Network manipulation for gap-filling process..")


universal_model = pickle.load(open("./input2/universal_model.p","rb"))

#From gapfill_network_manipulation.py
bigg_mnxr_dict = pickle.load(open("./input2/bigg_mnxr_dict.p","rb"))


logging.debug("Retrieving reaction information from target_model and universal_model..")
mnxr_bigg_target_model_dict = get_mnxr_bigg_in_target_model(target_model, bigg_mnxr_dict)

mnxr_unique_to_universal_model_list = get_mnxr_unique_to_universal_model(mnxr_bigg_target_model_dict, universal_model)


logging.debug("Merging target_model and universal_model..")
target_model2 = integrate_target_universal_models(mnxr_unique_to_universal_model_list, target_model, universal_model)


logging.debug("Step 2: Optimization-based gap-filling process..")

unique_nonprod_monomers_list = get_unique_nonprod_monomers_list(nonprod_sec_met_dict, prod_sec_met_dict)

#Some monomers used for nonproducible secondary metabolites can be produced from an initial target_model
#They need to be excluded from the list for gap-filling targets
adj_unique_nonprod_monomers_list = []

for nonprod_monomer in unique_nonprod_monomers_list:
    target_model_monomer = add_transport_exchange_rxn_nonprod_monomer(target_model, nonprod_monomer, dirname)
    target_model_monomer = check_producibility_nonprod_monomer(target_model_monomer, nonprod_monomer)
    if target_model_monomer.solution.f < 0.0001:
        adj_unique_nonprod_monomers_list.append(nonprod_monomer)
    else:
        continue

logging.debug("Adjusted unique_nonprod_monomers_list: %s" %adj_unique_nonprod_monomers_list)

for nonprod_monomer in adj_unique_nonprod_monomers_list:

    target_model_temp = add_transport_exchange_rxn_nonprod_monomer(target_model2, nonprod_monomer, dirname)
    target_model_temp = check_producibility_nonprod_monomer(target_model_temp, nonprod_monomer)
    target_model_temp.optimize()

    #Run gap-filling procedure only for monomers producible from target_model with reactions from universal_model
    if target_model_temp.solution.f > 0:

        added_reaction = execute_gapfill(target_model_temp, nonprod_monomer, mnxr_unique_to_universal_model_list)

        added_reaction2  = check_gapfill_rxn_biomass_effects(target_model, universal_model, added_reaction, template_exrxnid_flux_dict, dirname)

        target_model = add_gapfill_rxn_target_model(target_model, universal_model, added_reaction2)

    else:
        logging.debug("Gap-filling not possible: target_model with reactions from universal_model does not produce this monomer: %s" %nonprod_monomer)


#Cleanup of the final version of the target model
prune_unused_metabolites(target_model)

#Output
write_cobra_model_to_sbml_file(target_model, dirname+'4_complete_model/'+'%s_target_model_complete.xml' %orgname)

fp1 = open(dirname+'4_complete_model/'+'%s_target_model_complete_reactions.txt' %orgname, "w")
fp2 = open(dirname+'4_complete_model/'+'%s_target_model_complete_metabolites.txt' %orgname, "w")
fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Lower bound"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
fp2.write("Metabolite ID"+"\t"+"Metabolite name"+"\t"+"Formula"+"\t"+"Compartment"+"\n")

for j in range(len(target_model.reactions)):
    rxn = target_model.reactions[j]
    print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.lower_bound, rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)

for i in range(len(target_model.metabolites)):
    metab = target_model.metabolites[i]
    print >>fp2, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula, metab.compartment)

fp1.close()
fp2.close()

