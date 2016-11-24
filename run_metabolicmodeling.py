'''
2014-2015
Hyun Uk Kim, Tilmann Weber, Kyu-Sang Hwang and Jae Yong Ryu
'''

#Wildcard imports should never be used in production code.
import argparse
import logging
import multiprocessing
import os
import pickle
import sys
import time

from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
from cobra.manipulation.delete import prune_unused_metabolites
from argparse import Namespace
from modeling import prunPhase
from modeling import augPhase
from modeling import sec_met_rxn_generation
from modeling.gapfilling import gapfill_network_manipulation

def main():
    start = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--model', dest='orgName', default='sco', choices=['eco','sco'], help="Specify a template model for the target modeling")
    parser.add_argument('--disable-modeling', dest='pmr_generation', default=True, action='store_false', help='Disable primary metabolic modeling')
    parser.add_argument('-s', '--smr', dest='smr_generation', default=False, action=('store_true'), help="Specify whether to run secondary metabolic modeling")
    parser.add_argument('-i', '--input', dest='input', default='input', help="Specify input directory")
    parser.add_argument('-o', '--output', dest='output', default='output', help="Specify output directory")
    parser.add_argument('-e', '--ec', dest='eficaz', action='store_true', default=False, help="Run EC number prediction using EFICAz")
    parser.add_argument('-c', '--cpu', dest='cpus', default=multiprocessing.cpu_count(), type=int, help="How many CPUs to use in parallel. (default: %(default)s)")
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help="Print debugging information to stderr")

    options = parser.parse_args()

    if options.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    #Create output folders
    folders = ['0_EFICAz_results', '1_blastp_results', '2_primary_metabolic_model', '3_temp_models', '4_complete_model']

    for folder in folders:
        if not os.path.isdir(options.output+'/'+folder):
            os.makedirs(options.output+'/'+folder)

    options.outputfoldername = options.output+'/'+folders[0]

    get_genome_files(options)

    get_pickles_add_rxn(options)

    if options.pmr_generation:
        get_homolgs(options)

        model = get_pickles_prunPhase(options)

        modelPrunedGPR = run_prunPhase(model, options)

        if options.targetGenome_locusTag_ec_dict:

            target_model = run_augPhase(modelPrunedGPR, options)
        else:
            logging.debug("No EC_number found in the submitted gbk file")

        generate_outputs_primary_model(model, modelPrunedGPR, target_model, options)

    #Secondary metabolic modeling
    if options.smr_generation:
        for model_file in os.listdir(options.output+'/'+'2_primary_metabolic_model'):
            if model_file.endswith('.xml'):
                target_model = create_cobra_model_from_sbml_file(options.output+'/'+'2_primary_metabolic_model/'+model_file)

                logging.debug("Generating secondary metabolite biosynthesizing reactions..")
                cluster_nr = 1
                logging.debug("Total number of clusters: %s" %options.total_cluster)

                prod_sec_met_dict = {}
                nonprod_sec_met_dict = {}

                while cluster_nr <= options.total_cluster:
                    logging.debug("Cluster number: %s" %cluster_nr)
                    target_model = run_sec_met_rxn_generation(cluster_nr, target_model, prod_sec_met_dict, nonprod_sec_met_dict, options)
                    cluster_nr += 1

                target_model2, universal_model = prep_network_for_gapfilling(target_model, options)

                adj_unique_nonprod_monomers_list = get_target_nonprod_monomers_for_gapfilling(target_model, options)

                target_model_complete = run_gapfilling(target_model, target_model2, adj_unique_nonprod_monomers_list, universal_model, options)

                generate_outputs_secondary_model(target_model_complete, options)

            elif not model_file:
                logging.debug("COBRA-compliant SBML file needed")

    if not options.pmr_generation and not options.smr_generation:
        logging.debug("Either primary or secondary metabolic modeling should be performed")

    logging.debug(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


def get_genome_files(options):
    logging.debug("Reading input genome files..")

    logging.debug("Looking for a gbk file of a template model genome..")
    prunPhase.get_temp_fasta(options)

    logging.debug("Looking for a gbk file of a target genome..")
    prunPhase.get_target_gbk(options)

    logging.debug("Reading genbank file of the target genome..")
    prunPhase.get_targetGenomeInfo(options, 'genbank')

    logging.debug("Looking for a fasta file of a target genome..")
    prunPhase.get_target_fasta(options)


#For model augmentation phase in both primary and secondary modeling
def get_pickles_add_rxn(options):
    logging.debug("Loading pickle files necessary for the model augmentation phase..")

    bigg_mnxr_dict = pickle.load(open('./modeling/data/input2/bigg_mnxr_dict.p','rb'))
    options.bigg_mnxr_dict = bigg_mnxr_dict
    kegg_mnxr_dict = pickle.load(open('./modeling/data/input2/kegg_mnxr_dict.p','rb'))
    options.kegg_mnxr_dict = kegg_mnxr_dict
    mnxr_kegg_dict = pickle.load(open('./modeling/data/input2/mnxr_kegg_dict.p','rb'))
    options.mnxr_kegg_dict = mnxr_kegg_dict
    mnxr_rxn_dict = pickle.load(open('./modeling/data/input2/mnxr_rxn_dict.p','rb'))
    options.mnxr_rxn_dict = mnxr_rxn_dict

    bigg_mnxm_compound_dict = pickle.load(open('./modeling/data/input2/bigg_mnxm_compound_dict.p','rb'))
    options.bigg_mnxm_compound_dict = bigg_mnxm_compound_dict
    mnxm_bigg_compound_dict = pickle.load(open('./modeling/data/input2/mnxm_bigg_compound_dict.p','rb'))
    options.mnxm_bigg_compound_dict = mnxm_bigg_compound_dict
    kegg_mnxm_compound_dict = pickle.load(open('./modeling/data/input2/kegg_mnxm_compound_dict.p','rb'))
    options.kegg_mnxm_compound_dict = kegg_mnxm_compound_dict
    mnxm_kegg_compound_dict = pickle.load( open('./modeling/data/input2/mnxm_kegg_compound_dict.p','rb'))
    options.mnxm_kegg_compound_dict = mnxm_kegg_compound_dict

    mnxm_compoundInfo_dict = pickle.load(open('./modeling/data/input2/mnxm_compoundInfo_dict.p','rb'))
    options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    template_exrxnid_flux_dict = pickle.load(open('%s/tempModel_exrxnid_flux_dict.p' %(options.input1),'rb'))
    options.template_exrxnid_flux_dict = template_exrxnid_flux_dict


def get_homolgs(options):
    logging.debug("Searching bidirectional homolgs..")
    logging.debug("Generating a DB for the genes from the target genome..")
    prunPhase.make_blastDB(options)

    logging.debug("Running BLASTP #1: genes in the target genome against genes in the template model..")
    prunPhase.run_blastp(target_fasta='./%s/1_blastp_results/targetGenome_locusTag_aaSeq.fa' %options.output, blastp_result='./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, db_dir = '%s/tempBlastDB' %(options.input1), evalue=1e-30)

    logging.debug("Running BLASTP #2: genes in the template model against genes in the target genome..")
    prunPhase.run_blastp(target_fasta=options.temp_fasta, blastp_result='./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, db_dir = './%s/1_blastp_results/targetBlastDB' %options.output, evalue=1e-30)

    logging.debug("Parsing the results of BLASTP #1..")
    blastpResults_dict1 = prunPhase.parseBlaspResults('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, './%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)

    logging.debug("Parsing the results of BLASTP #2..")
    blastpResults_dict2 = prunPhase.parseBlaspResults('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, './%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

    logging.debug("Selecting the best hits for BLASTP #1..")
    bestHits_dict1 = prunPhase.makeBestHits_dict('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)

    logging.debug("Selecting the best hits for BLASTP #2..")
    bestHits_dict2 = prunPhase.makeBestHits_dict('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

    logging.debug("Selecting the bidirectional best hits..")
    prunPhase.getBBH(bestHits_dict1, bestHits_dict2, options)

    logging.debug("Selecting genes that are not bidirectional best hits..")
    prunPhase.get_nonBBH(options)


#For model pruning phase
#Only model file is not saved in Namespace
def get_pickles_prunPhase(options):
    logging.debug("Loading pickle files of the parsed template model and its relevant genbank data..")
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    tempModel_biggRxnid_locusTag_dict = pickle.load(open('%s/tempModel_biggRxnid_locusTag_dict.p' %(options.input1),'rb'))
    options.tempModel_biggRxnid_locusTag_dict = tempModel_biggRxnid_locusTag_dict

    return model


def run_prunPhase(model, options):
    logging.debug("Pruning phase starting..")
    logging.debug("Labeling reactions with nonhomologous genes to remove from the template model..")
    prunPhase.labelRxnToRemove(model, options)

    logging.debug("Removing reactions with nonhomologous genes from the template model..")
    modelPruned = prunPhase.pruneModel(model, options, 'gurobi')

    logging.debug("Correcting GPR associations in the template model..")
    modelPrunedGPR = prunPhase.swap_locusTag_tempModel(modelPruned, options)

    return modelPrunedGPR


def run_augPhase(modelPrunedGPR, options):
    logging.debug("Augmentation phase starting..")
    logging.debug("Creating various dictionary files for the nonBBH gene-associted reactions...")

    augPhase.get_targetGenome_locusTag_ec_nonBBH_dict(options)

    #Two nested functions
    #def get_rxnid_from_ECNumber(enzymeEC):
    #def get_rxnInfo_from_rxnid(rxnid):
    augPhase.make_all_rxnInfo_fromRefSeq(options)

    augPhase.get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, options)

    logging.debug("Adding the nonBBH gene-associated reactions...")
    augPhase.check_existing_rxns(options)

    augPhase.get_mnxr_using_kegg(options)

    augPhase.extract_rxn_mnxm_coeff(options)

    target_model = augPhase.add_nonBBH_rxn(modelPrunedGPR, options)

    return target_model


def generate_outputs_primary_model(model, modelPrunedGPR, target_model, options):
    #Output files
    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    #Cobrapy IO module seems to have an error for adding new reactions
    write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/target_model_%s.xml' %(options.output, options.orgName))
    target_model = create_cobra_model_from_sbml_file('./%s/2_primary_metabolic_model/target_model_%s.xml' %(options.output, options.orgName))
    write_cobra_model_to_sbml_file(target_model, './%s/2_primary_metabolic_model/target_model_%s.xml' %(options.output, options.orgName))

    #Output on screen
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    logging.debug("Number of genes: %s; %s; %s" %(len(model.genes), len(modelPrunedGPR.genes), len(target_model.genes)))
    logging.debug("Number of reactions: %s; %s; %s" %(len(model.reactions), len(modelPrunedGPR.reactions), len(target_model.reactions)))
    logging.debug("Number of metabolites: %s; %s; %s" %(len(model.metabolites), len(modelPrunedGPR.metabolites), len(target_model.metabolites)))

    fp1 = open('./%s/2_primary_metabolic_model/target_model_reactions.txt' %options.output, "w")
    fp2 = open('./%s/2_primary_metabolic_model/target_model_metabolites.txt' %options.output, "w")
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


def run_sec_met_rxn_generation(cluster_nr, target_model, prod_sec_met_dict, nonprod_sec_met_dict, options):

    sec_met_rxn_generation.get_cluster_location(cluster_nr, options)

    sec_met_rxn_generation.get_cluster_info_from_seq_record(options)

    sec_met_rxn_generation.get_cluster_product(cluster_nr, options)

    if 't1pks' in options.product or 'nrps' in options.product:
        sec_met_rxn_generation.get_cluster_domain(options)

        sec_met_rxn_generation.get_cluster_monomers(options)

        sec_met_rxn_generation.get_cluster_module(options)

        sec_met_rxn_generation.get_currency_metabolites(options)

        sec_met_rxn_generation.get_total_currency_metab_coeff(options)

        sec_met_rxn_generation.get_all_metab_coeff(options)

        target_model = sec_met_rxn_generation.add_sec_met_rxn(target_model, options)

        target_model = sec_met_rxn_generation.check_producibility_sec_met(target_model, options)

        if target_model.solution.f < 0.0001:
            nonprod_sec_met_metab_list = sec_met_rxn_generation.get_monomers_nonprod_sec_met(options)
            nonprod_sec_met_dict[options.product] = nonprod_sec_met_metab_list
        else:
            prod_sec_met_metab_list = sec_met_rxn_generation.get_monomers_prod_sec_met(options)
            prod_sec_met_dict[options.product] = prod_sec_met_metab_list

    else:
        logging.debug("Not type I polyketide synthase ('t1pks'), nonribosomal synthetase ('nrps') or their hybird")

    if cluster_nr == options.total_cluster:
        options.prod_sec_met_dict = prod_sec_met_dict
        options.nonprod_sec_met_dict = nonprod_sec_met_dict

    return target_model


def prep_network_for_gapfilling(target_model, options):

    logging.debug("Gap-filling for the production of secondary metabolites..")
    logging.debug("Step 1: Network manipulation for gap-filling process..")

    universal_model = pickle.load(open("./modeling/data/input2/universal_model.p","rb"))

    logging.debug("Retrieving reaction information from target_model and universal_model..")
    gapfill_network_manipulation.get_mnxr_bigg_in_target_model(target_model, options)

    gapfill_network_manipulation.get_mnxr_unique_to_universal_model(universal_model, options)

    logging.debug("Merging target_model and universal_model..")
    target_model2 = gapfill_network_manipulation.integrate_target_universal_models(target_model, universal_model, options)

    return target_model2, universal_model


def get_target_nonprod_monomers_for_gapfilling(target_model, options):
    logging.debug("Step 2: Optimization-based gap-filling process..")

    unique_nonprod_monomers_list = gapfill_network_manipulation.get_unique_nonprod_monomers_list(options)

    #Some monomers used for nonproducible secondary metabolites can be produced from an initial target_model
    #They need to be excluded from the list for gap-filling targets
    adj_unique_nonprod_monomers_list = []

    for nonprod_monomer in unique_nonprod_monomers_list:
        target_model_monomer = gapfill_network_manipulation.add_transport_exchange_rxn_nonprod_monomer(target_model, nonprod_monomer, options)
        target_model_monomer = gapfill_network_manipulation.check_producibility_nonprod_monomer(target_model_monomer, nonprod_monomer)
        if target_model_monomer.solution.f < 0.0001:
            adj_unique_nonprod_monomers_list.append(nonprod_monomer)
        else:
            continue

    logging.debug("Adjusted unique_nonprod_monomers_list: %s" %adj_unique_nonprod_monomers_list)

    return adj_unique_nonprod_monomers_list


def run_gapfilling(target_model, target_model2, adj_unique_nonprod_monomers_list, universal_model, options):

    for nonprod_monomer in adj_unique_nonprod_monomers_list:
        target_model_temp = gapfill_network_manipulation.add_transport_exchange_rxn_nonprod_monomer(target_model2, nonprod_monomer, options)
        target_model_temp = gapfill_network_manipulation.check_producibility_nonprod_monomer(target_model_temp, nonprod_monomer)
        target_model_temp.optimize()

        #Run gap-filling procedure only for monomers producible from target_model with reactions from universal_model
        if target_model_temp.solution.f > 0:
            added_reaction = gapfill_network_manipulation.execute_gapfill(target_model_temp, nonprod_monomer, options)
            added_reaction2  = gapfill_network_manipulation.check_gapfill_rxn_biomass_effects(target_model, universal_model, added_reaction, options)
            target_model_complete = gapfill_network_manipulation.add_gapfill_rxn_target_model(target_model, universal_model, added_reaction2, options)
        else:
            logging.debug("Gap-filling not possible: target_model with reactions from universal_model does not produce this monomer: %s" %nonprod_monomer)

    #Cleanup of the final version of the target model
    prune_unused_metabolites(target_model_complete)

    return target_model_complete


def generate_outputs_secondary_model(target_model_complete, options):
    #Output files
    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    #Cobrapy IO module seems to have an error for adding new reactions
    write_cobra_model_to_sbml_file(target_model_complete, './%s/4_complete_model/target_model_complete.xml' %options.output)

    #Output on screen
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    logging.debug("Number of genes: %s" %len(target_model_complete.genes))
    logging.debug("Number of reactions: %s" %len(target_model_complete.reactions))
    logging.debug("Number of metabolites: %s" %len(target_model_complete.metabolites))

    fp1 = open('./%s/4_complete_model/target_model_complete_reactions.txt' %options.output, "w")
    fp2 = open('./%s/4_complete_model/target_model_complete_metabolites.txt' %options.output, "w")
    fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Lower bound"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
    fp2.write("Metabolite ID"+"\t"+"Metabolite name"+"\t"+"Formula"+"\t"+"Compartment"+"\n")

    for j in range(len(target_model_complete.reactions)):
        rxn = target_model.reactions[j]
        print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.lower_bound, rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)

    for i in range(len(target_model_complete.metabolites)):
        metab = target_model.metabolites[i]
        print >>fp2, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula, metab.compartment)

    fp1.close()
    fp2.close()


if __name__ == '__main__':
    main()
