'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
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
    check_producibility_sec_met
)
from gapfill_network_manipulation import (
    get_mnxr_bigg_in_target_model,
    get_mnxr_unique_to_universal_model,
    get_balanced_rxns_from_mnxr,
    get_manipulated_target_universal_models
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
#inputfile = './NC_021055.1.cluster002.gbk' #NRPS
#inputfile = './NC_013929.1.cluster031.gbk' #PKS
#inputfile = './NC_020990.1.cluster023.gbk' #Hybrid

#if __name__ == '__main__':
#    cluster_f = 'NC_018750.1.cluster003.gbk'

nonproducible_sec_met = []

for cluster_f in cluster_files:
    print '\n', cluster_f

    cluster_info_dict, record = get_cluster_info_from_cluster_gbk(dirname+cluster_f, "genbank")

    product = get_product_from_cluster_gbk(record)

    if 't1pks' in product or 'nrps' in product:

        locustag_domain_dict = get_cluster_domain(cluster_info_dict)

        locustag_monomer_dict = get_cluster_monomers(cluster_info_dict)

        locustag_module_domain_dict = get_cluster_module(locustag_domain_dict)

        module_currency_metab_dict = get_currency_metabolites(locustag_module_domain_dict)

        currency_metab_coeff_dict = get_total_currency_metab_coeff(module_currency_metab_dict)

        metab_coeff_dict = get_all_metab_coeff(locustag_monomer_dict, currency_metab_coeff_dict, product)

        target_model = add_sec_met_rxn(target_model, metab_coeff_dict, product, bigg_mnxm_compound_dict, mnxm_compoundInfo_dict, cluster_info_dict)
        
        product = check_producibility_sec_met(dirname, orgname, target_model, metab_coeff_dict, product)
        if product != None:
            nonproducible_sec_met.append(product)

print "\n", "Nonproducible secondary metabolites:"
print nonproducible_sec_met

'''
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


print "Gap-filling for the production of secondary metabolites..."

target_model = create_cobra_model_from_sbml_file('./E_coli_iAF1260_rewrite.xml')
#cobra_model = create_cobra_model_from_sbml_file('./sma_target_model_sco_complete.xml')

#ivcoa_c, 3-Methylbutanoyl-CoA: one of the precursors for non-ribosomal peptide
if 'ivcoa_c' in target_model.metabolites:
    print "\n"
    print "ivcoa_c (ivcoa_c) is already included in template_model"
else:
    print "\n"
    print "'Metabolite (ivcoa_c) is not included in template_model"

#Load universal network in sbml format
#universal_model = create_cobra_model_from_sbml_file('./universal_network_fiexed_bigg_mnxref.xml')

# fixing the definition of metabolite id (e.g. ala__L_c --> ala_DASH_L_c)
#universal_model = fix_special_characters_compoundid(universal_model)

#pickle.dump(universal_model, open('./input/universal_model.p','wb'))
universal_model = pickle.load(open("./input/universal_model.p","rb"))

#From gapfill_network_manipulation.py
bigg_mnxr_dict = pickle.load(open("./input/bigg_mnxr_dict.p","rb"))
gg_target_model_dict = get_mnxr_bigg_in_target_model(target_model, bigg_mnxr_dict)

mnxr_unique_to_universal_model_list = get_mnxr_unique_to_universal_model(mnxr_bigg_target_model_dict, universal_model)

mnxr_rxn_all_dict = pickle.load(open("./input/mnxr_rxn_all_dict.p","rb"))

balanced_unique_mnxr_list = get_balanced_rxns_from_mnxr(mnxr_unique_to_universal_model_list, mnxr_rxn_all_dict)

target_model, universal_model2 = get_manipulated_target_universal_models(balanced_unique_mnxr_list, target_model, universal_model)


#Run gap-filling algorithm based on MILP in gurobipy
obj = gapfilling_precursor()

#Load merged model
obj.load_cobra_model(target_model)
obj.change_reversibility('Biomass_SCO', target_model)
obj.fill_gap(demand_reaction.id, target_model, universal_model2)


write_cobra_model_to_sbml_file(target_model, dirname+model_sbml[:-4]+'_complete2.xml')

fp1 = open(dirname+'%s_target_model_reactions2.txt' %orgname, "w")
fp2 = open(dirname+'%s_target_model_metabolites2.txt' %orgname, "w")
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
'''
