
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import collections
import logging
import re
from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file


def generate_outputs(folder, cobra_model, runtime, options):
    #Output files
    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    #Cobrapy IO module seems to have an error for adding new reactions
    write_cobra_model_to_sbml_file(cobra_model,
            './%s/%s/model.xml' %(options.outputfolder, folder))
    cobra_model = create_cobra_model_from_sbml_file(
            './%s/%s/model.xml' %(options.outputfolder, folder))
    write_cobra_model_to_sbml_file(cobra_model,
            './%s/%s/model.xml' %(options.outputfolder, folder))

    num_essen_rxn, num_kegg_rxn, num_cluster_rxn = get_model_reactions(
                       folder, cobra_model, options)
    get_model_metabolites(folder, cobra_model, options)
    template_model_gene_list = get_model_genes(folder, cobra_model, options)
    get_summary_report(folder, cobra_model, runtime,
                       num_essen_rxn, num_kegg_rxn, num_cluster_rxn,
                       template_model_gene_list, options)

    if folder == '2_primary_metabolic_model':
        logging.info("'Primary' metabolic model completed")
    elif folder == '4_complete_model':
        logging.info("'Secondary' metabolic model completed")


def get_model_reactions(folder, cobra_model, options):

    fp1 = open('./%s/%s/model_reactions.txt'
            %(options.outputfolder, folder), 'w')
    fp2 = open('./%s/%s/remaining_essential_reactions_from_template_model.txt'
            %(options.outputfolder, folder), 'w')
    fp3 = open('./%s/%s/reactions_added_from_kegg.txt'
            %(options.outputfolder, folder), 'w')

    fp1.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')
    fp2.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')
    fp3.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')

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

        #Reactions added from KEGG
        if re.search('R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]', rxn.id):
            num_kegg_rxn+=1
            print >>fp3, '%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem)

        #Secondary metabolite biosynthetic reactions
        if re.search('Cluster', rxn.id):
            num_cluster_rxn+=1

    fp1.close()
    fp2.close()
    fp3.close()

    return num_essen_rxn, num_kegg_rxn, num_cluster_rxn


def get_model_metabolites(folder, cobra_model, options):

    fp1 = open('./%s/%s/model_metabolites.txt' %(options.outputfolder, folder), "w")

    fp1.write('metabolite_ID'+'\t'+'metabolite_name'+'\t'+'formula'+'\t'+'compartment'+'\n')

    for i in range(len(cobra_model.metabolites)):
        metab = cobra_model.metabolites[i]
        print >>fp1, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula,
                metab.compartment)

    fp1.close()


def get_model_genes(folder, cobra_model, options):

    template_model_gene_list = []

    fp1 = open('./%s/%s/remaining_genes_from_template_model.txt'
                %(options.outputfolder, folder), "w")

    fp1.write('remaining_gene_from_template_model'+'\t'+'reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')

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
    fp1 = open('./%s/%s/summary_report.txt' %(options.outputfolder, folder), "w")

    fp1.write(''+'\t'+'primary_metabolic_model'+'\n')

    if options.verbose:
        log_level = 'verbose'
    elif options.debug:
        log_level = 'debug'

    model_summary_dict = {}
    model_summary_dict['template_model_organism']= options.orgName
    model_summary_dict['input_file']=options.input
    model_summary_dict['outputfolder']=options.outputfolder
    model_summary_dict['eficaz']=options.eficaz
    model_summary_dict['number_cpu_use']=options.cpus
    model_summary_dict['log_level']=options.log_level
    model_summary_dict['primary_metabolic_modeling']=options.pmr_generation
    model_summary_dict['secondary_metabolic_modeling']=options.smr_generation
    model_summary_dict['runtime']=runtime
    model_summary_dict['number_genes']=len(cobra_model.genes)
    model_summary_dict['number_reactions']=len(cobra_model.reactions)
    model_summary_dict['number_metabolites']=len(cobra_model.metabolites)
    model_summary_dict['number_remaining_essential_reactions_from_template_model'] \
        =num_essen_rxn
    model_summary_dict['number_reactions_added_from_kegg']=num_kegg_rxn
    model_summary_dict['number_remaining_genes_from_template_model'] \
        =len(template_model_gene_list)
    model_summary_dict['number_clusters_for_reactions']=num_cluster_rxn

    model_summary_dict2 = collections.OrderedDict(sorted(model_summary_dict.items()))

    for key in model_summary_dict.keys():
        print >>fp1, '%s\t%s\t' %(key, model_summary_dict2[key])

    fp1.close()

