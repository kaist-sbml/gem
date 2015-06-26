'''
2015 Hyun Uk Kim
2015 Kyu-Sang Hwang
'''

from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
from cobra import Model, Reaction, Metabolite
import pickle
import os
import sys
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
    add_gapfill_rxn_target_model,
)
from gapfill_core import gapfilling_precursor


print "Generating secondary metabolite biosynthesizing reactions.."

orgname = sys.argv[1]

bigg_mnxm_compound_dict = pickle.load(open('./input/bigg_mnxm_compound_dict.p','rb'))
mnxm_compoundInfo_dict = pickle.load(open('./input/mnxm_compoundInfo_dict.p','rb'))

dirname = './%s/' %orgname
cluster_files = []
for f in os.listdir(dirname):
    if f.endswith('.gbk') and 'cluster' in f:
        cluster_files.append(f)

cluster_files.sort()

for f in os.listdir(dirname):
    if f.endswith('.xml'):
       model_sbml = f

target_model = create_cobra_model_from_sbml_file(dirname+model_sbml)

#if __name__ == '__main__':
#    cluster_f = 'NC_018750.1.cluster003.gbk'

prod_sec_met_dict = {}
nonprod_sec_met_dict = {}

for cluster_f in cluster_files:
    print '\n', cluster_f

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
        
        #target_model.solution.f, product = check_producibility_sec_met(target_model, product, dirname)
        target_model, product = check_producibility_sec_met(target_model, product, dirname)

        if target_model.solution.f < 0.0001:
            nonprod_sec_met_metab_list = get_monomers_nonprod_sec_met(metab_coeff_dict)
            nonprod_sec_met_dict[product] = nonprod_sec_met_metab_list
        else:
            prod_sec_met_metab_list = get_monomers_prod_sec_met(metab_coeff_dict)
            prod_sec_met_dict[product] = prod_sec_met_metab_list

print "\n", "Producible secondary metabolites:"
print prod_sec_met_dict, "\n"

print "\n", "Nonproducible secondary metabolites:"
print nonprod_sec_met_dict, "\n"


print "Gap-filling for the production of secondary metabolites.."
print "Step 1: Network manipulation for gap-filling process..", "\n"


universal_model = pickle.load(open("./input/universal_model.p","rb"))

#From gapfill_network_manipulation.py
bigg_mnxr_dict = pickle.load(open("./input/bigg_mnxr_dict.p","rb"))


print "Retrieving reaction information from target_model and universal_model.."
mnxr_bigg_target_model_dict = get_mnxr_bigg_in_target_model(target_model, bigg_mnxr_dict)

mnxr_unique_to_universal_model_list = get_mnxr_unique_to_universal_model(mnxr_bigg_target_model_dict, universal_model)

mnxr_rxn_all_dict = pickle.load(open("./input/mnxr_rxn_all_dict.p","rb"))

print "Merging target_model and universal_model.."
print "Also generating a truncated universal_model with its exclusive reactions.."
print "\n"
target_model2 = integrate_target_universal_models(mnxr_unique_to_universal_model_list, target_model, universal_model)


print "Step 2: Optimization-based gap-filling process..", "\n"

unique_nonprod_monomers_list = get_unique_nonprod_monomers_list(nonprod_sec_met_dict, prod_sec_met_dict)

#Some monomers used for nonproducible secondary metabolites can be produced from an initial target_model
#They need to be excluded from the list for gap-filling targets
for nonprod_monomer in unique_nonprod_monomers_list:
    print nonprod_monomer
    target_model_monomer = add_transport_exchange_rxn_nonprod_monomer(target_model, nonprod_monomer)
    target_model_monomer = check_producibility_nonprod_monomer(target_model_monomer, nonprod_monomer, dirname)
    if target_model_monomer.solution.f > 0:
        print "Optimal value for", nonprod_monomer, ":", target_model_monomer.solution.f
        unique_nonprod_monomers_list.remove(nonprod_monomer)

    else:
        continue

print "Adjusted unique_nonprod_monomers_list", unique_nonprod_monomers_list, "\n"

for nonprod_monomer in unique_nonprod_monomers_list:

    target_model_temp = add_transport_exchange_rxn_nonprod_monomer(target_model2, nonprod_monomer)
    target_model_temp = check_producibility_nonprod_monomer(target_model_temp, nonprod_monomer, dirname)
    target_model_temp.optimize()

    #Run gap-filling procedure only for monomers producible from target_model with reactions from universal_model
    if target_model_temp.solution.f > 0:

        #Run gap-filling algorithm based on MILP in gurobipy
        obj = gapfilling_precursor()

        #Load merged model
        obj.load_cobra_model(target_model_temp)
        #obj.change_reversibility(target_model_temp.reactions.get_by_id('Ex_'+nonprod_monomer), target_model_temp)
        added_reaction = obj.fill_gap("Ex_"+nonprod_monomer, mnxr_unique_to_universal_model_list)
        target_model = add_gapfill_rxn_target_model(target_model, universal_model, added_reaction)
        print "\n"
    else:
        print "Gap-filling not possible: target_model with reactions from universal_model does not produce this monomer", nonprod_monomer, "\n"


#Output
write_cobra_model_to_sbml_file(target_model, dirname+model_sbml[:-4]+'_complete.xml')

fp1 = open(dirname+'%s_target_model_reactions.txt' %orgname, "w")
fp2 = open(dirname+'%s_target_model_metabolites.txt' %orgname, "w")
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

