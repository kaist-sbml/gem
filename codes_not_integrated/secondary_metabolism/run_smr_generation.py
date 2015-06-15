'''
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

from cobra.io.sbml import create_cobra_model_from_sbml_file, write_cobra_model_to_sbml_file
import pickle
import os
from sec_met_rxn_generation import (
    get_product_from_cluster_gbk,
    get_cluster_info_from_cluster_gbk,
    get_cluster_domain, get_cluster_monomers,
    get_cluster_module, get_currency_metabolites,
    get_total_currency_metab_coeff,
    get_all_metab_coeff,
    add_sec_met_rxn
)
import sys

print "Generating secondary metabolite biosynthesizing reactions.."

orgname = sys.argv[1]

bigg_mnxm_compound_dict = pickle.load(open('bigg_mnxm_compound_dict.p','rb'))
mnxm_compoundInfo_dict = pickle.load(open('mnxm_compoundInfo_dict.p','rb'))

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
