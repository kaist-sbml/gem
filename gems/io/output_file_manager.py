
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import collections
import logging
import pickle
import re
from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file


def generate_outputs(folder, cobra_model_no_gapFilled, cobra_model, runtime, options):
    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    write_cobra_model_to_sbml_file(cobra_model,
            './%s/model.xml' %folder, use_fbc_package=False)
    cobra_model = create_cobra_model_from_sbml_file(
            './%s/model.xml' %folder)
    write_cobra_model_to_sbml_file(cobra_model,
            './%s/model.xml' %folder, use_fbc_package=False)

    num_essen_rxn, num_kegg_rxn, num_cluster_rxn = get_model_reactions(
                       folder, cobra_model_no_gapFilled, cobra_model, options)
    get_model_metabolites(folder, cobra_model, options)
    template_model_gene_list = get_model_genes(folder, cobra_model)
    get_summary_report(folder, cobra_model, runtime,
                       num_essen_rxn, num_kegg_rxn, num_cluster_rxn,
                       template_model_gene_list, options)

    if '3_primary_metabolic_model'in folder:
        logging.info("'Primary' metabolic model completed")
    elif '4_complete_model' in folder:
        logging.info("'Secondary' metabolic model completed")

    if options.pmr_generation and options.debug:
        write_data_for_debug(options)


def get_model_reactions(folder, cobra_model_no_gapFilled, cobra_model, options):

    fp1 = open('./%s/model_reactions.txt' %folder, 'w')
    fp2 = open('./%s/remaining_essential_reactions_from_template_model.txt' %folder, 'w')
    fp3 = open('./%s/reactions_added_from_kegg.txt' %folder, 'w')
    fp4 = open('./%s/cluster_fluxes.txt' %folder, 'w')

    fp1.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'
            +'GPR'+'\t'+'pathway'+'\n')
    fp2.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'
            +'GPR'+'\t'+'pathway'+'\n')
    fp3.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'
            +'GPR'+'\t'+'pathway'+'\n')
    fp4.write('reaction_ID'+'\t'+'fluxes without gap-filling reactions'+'\t'
            +'fluxes with gap-filling reactions'+'\n')

    num_essen_rxn = 0
    num_kegg_rxn = 0
    num_cluster_rxn = 0
    for j in range(len(cobra_model.reactions)):
        rxn = cobra_model.reactions[j]
        print >>fp1, '%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem)

        #Remaining essential reactions
        if 'rxnToRemove_dict' in options:
            if rxn.id in options.rxnToRemove_dict.keys():
                if options.rxnToRemove_dict[rxn.id] == False:
                    num_essen_rxn+=1
                    print >>fp2, '%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem)
        else:
            print >>fp2, 'Primary metabolic modeling not performed'
            num_essen_rxn = 'Primary metabolic modeling not performed'

        #Reactions added from KEGG
        if re.search('R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]', rxn.id):
            num_kegg_rxn+=1
            print >>fp3, '%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem)

        #Secondary metabolite biosynthetic reactions
        if re.search('Ex_Cluster', rxn.id):
            num_cluster_rxn+=1

            #Calculated flux values are inaccurate without
            #manual setting of objective_coefficient to zero
            if 'Biomass_SCO' in cobra_model_no_gapFilled.reactions:
                cobra_model_no_gapFilled.reactions.get_by_id('Biomass_SCO').objective_coefficient = 0
            elif 'Ec_biomass_iAF1260_core_59p81M' in cobra_model_no_gapFilled.reactions:
                cobra_model_no_gapFilled.reactions.get_by_id('Ec_biomass_iAF1260_core_59p81M').objective_coefficient = 0

            cobra_model_no_gapFilled.reactions.get_by_id(rxn.id).objective_coefficient = 1
            cobra_model_no_gapFilled.optimize()

            if 'Biomass_SCO' in cobra_model.reactions:
                cobra_model.reactions.get_by_id('Biomass_SCO').objective_coefficient = 0
            elif 'Ec_biomass_iAF1260_core_59p81M' in cobra_model.reactions:
                cobra_model.reactions.get_by_id('Ec_biomass_iAF1260_core_59p81M').objective_coefficient = 0

            cobra_model.reactions.get_by_id(rxn.id).objective_coefficient = 1
            cobra_model.optimize()

            print >>fp4, '%s\t%f\t%f' \
                %(rxn.id, cobra_model_no_gapFilled.solution.f, cobra_model.solution.f)

            cobra_model_no_gapFilled.reactions.get_by_id(rxn.id).objective_coefficient = 0
            cobra_model.reactions.get_by_id(rxn.id).objective_coefficient = 0

    fp1.close()
    fp2.close()
    fp3.close()
    fp4.close()

    return num_essen_rxn, num_kegg_rxn, num_cluster_rxn


def get_model_metabolites(folder, cobra_model, options):

    fp1 = open('./%s/model_metabolites.txt' %folder, "w")
    fp2 = open('./%s/metabolites_gapfilling_needed.txt' %folder, "w")

    fp1.write('metabolite_ID'+'\t'+'metabolite_name'+'\t'
            +'formula'+'\t'+'compartment'+'\n')
    fp2.write('metabolite_ID'+'\t'+'reaction_ID'+'\t'+'reaction_name'+'\t'
            +'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')

    for i in range(len(cobra_model.metabolites)):
        metab = cobra_model.metabolites[i]
        print >>fp1, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula,
                metab.compartment)

        if '4_complete_model' in folder:
            #Remove compartment suffix (e.g., '_c') from 'metab.id'
            if metab.id[:-2] in options.adj_unique_nonprod_monomers_list:
                logging.debug("Metabolite for gap-filling: %s" %metab.id)

                for j in range(len(cobra_model.reactions)):
                    rxn = cobra_model.reactions[j]

                    if metab.id in rxn.reaction:
                        logging.debug("Relevant reactions: %s" %rxn.id)
                        print >>fp2, '%s\t%s\t%s\t%s\t%s\t%s' %(metab.id,
                            rxn.id, rxn.name, rxn.reaction,
                            rxn.gene_reaction_rule, rxn.subsystem)

    fp1.close()
    fp2.close()


def get_model_genes(folder, cobra_model):

    template_model_gene_list = []

    fp1 = open('./%s/remaining_genes_from_template_model.txt' %folder, "w")

    fp1.write('remaining_gene_from_template_model'+'\t'+'reaction_ID'+'\t'
            +'reaction_name'+'\t'+'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')

    for g in range(len(cobra_model.genes)):
        gene = cobra_model.genes[g]

        if re.search('SCO[0-9]+[0-9]+[0-9]+[0-9]', str(gene)) \
            or re.search('b[0-9]+[0-9]+[0-9]+[0-9]', str(gene)):
            template_model_gene_list.append(gene)

            for j in range(len(cobra_model.reactions)):
                rxn = cobra_model.reactions[j]

                if str(gene) in rxn.gene_reaction_rule:
                    print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(str(gene), rxn.id, rxn.name,
                            rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)

    fp1.close()
    return template_model_gene_list


def get_summary_report(folder, cobra_model, runtime,
                       num_essen_rxn, num_kegg_rxn, num_cluster_rxn,
                       template_model_gene_list, options):
    fp1 = open('./%s/summary_report.txt' %folder, "w")

    #log_level
    if options.verbose:
        log_level = 'verbose'
    elif options.debug:
        log_level = 'debug'

    #runtime
    runtime2 = runtime.split()[2]

    model_summary_dict = {}
    model_summary_dict['template_model_organism']= options.orgName
    model_summary_dict['input_file']=options.input
    model_summary_dict['outputfolder']=options.outputfolder
    model_summary_dict['eficaz']=options.eficaz
    model_summary_dict['number_cpu_use']=options.cpus
    model_summary_dict['log_level']=log_level
    model_summary_dict['primary_metabolic_modeling']=options.pmr_generation
    model_summary_dict['secondary_metabolic_modeling']=options.smr_generation
    model_summary_dict['runtime']=runtime2
    model_summary_dict['number_genes']=len(cobra_model.genes)
    model_summary_dict['number_reactions']=len(cobra_model.reactions)
    model_summary_dict['number_metabolites']=len(cobra_model.metabolites)
    model_summary_dict['number_remaining_essential_reactions_from_template_model'] \
        =num_essen_rxn
    model_summary_dict['number_reactions_added_from_kegg']=num_kegg_rxn
    model_summary_dict['number_clusters_for_reactions']=num_cluster_rxn
    model_summary_dict['number_remaining_genes_from_template_model'] \
        =len(template_model_gene_list)

    if '4_complete_model' in folder:
        model_summary_dict['number_metabolites_for_gapfilling'] \
            =len(options.adj_unique_nonprod_monomers_list)
    else:
        model_summary_dict['number_metabolites_for_gapfilling']=0

    #Sort data by keys
    model_summary_dict2 = collections.OrderedDict(sorted(model_summary_dict.items()))

    for key in model_summary_dict2.keys():
        print >>fp1, '%s\t%s' %(key, model_summary_dict2[key])

    fp1.close()


def write_data_for_debug(options):

    with open('./%s/temp_target_BBH_dict.txt' %options.outputfolder2,'w') as f1:
        for locustag in options.temp_target_BBH_dict.keys():
            print >> f1, '%s\t%s' %(locustag, options.temp_target_BBH_dict[locustag])

    pickle.dump(options.temp_target_BBH_dict,
            open('./%s/temp_target_BBH_dict.p' %options.outputfolder2,'wb'))

    with open('./%s/rxnid_to_add_list.txt' %options.outputfolder5,'w') as f2:
        for rxnid in options.rxnid_to_add_list:
            print >>f2, '%s' %rxnid

    with open('./%s/mnxr_to_add_list.txt' %options.outputfolder5,'w') as f3:
        for mnxr in options.mnxr_to_add_list:
            print >>f3, '%s' %mnxr

    with open('./%s/rxnid_info_dict.txt' %options.outputfolder5,'w') as f4:
        for rxnid in options.rxnid_info_dict.keys():
            print >>f4, '%s' %rxnid

    with open('./%s/rxnid_mnxm_coeff_dict.txt' %options.outputfolder5,'w') as f5:
        for rxnid in options.rxnid_mnxm_coeff_dict.keys():
            print >>f5, '%s\t%s' %(rxnid, options.rxnid_mnxm_coeff_dict[rxnid])
