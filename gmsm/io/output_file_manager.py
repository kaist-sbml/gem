
import cobra
import collections
import logging
import pickle
import re
import shutil
from cobra.util.solver import linear_reaction_coefficients
from gmsm import utils

def generate_outputs(folder, runtime, options, **kwargs):

    if 'cobra_model' in kwargs:
        cobra_model = kwargs['cobra_model']

    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    cobra_model = utils.stabilize_model(cobra_model, folder, '')

    num_essen_rxn, num_kegg_rxn, num_secondary_rxn = get_model_reactions(
                       folder, options, **kwargs)
    get_model_metabolites(folder, cobra_model, options)
    template_model_gene_list, duplicate_gene_list = \
                       get_model_genes(folder, cobra_model, options)
    get_summary_report(folder, cobra_model, runtime,
                       num_essen_rxn, num_kegg_rxn, num_secondary_rxn,
                       template_model_gene_list, duplicate_gene_list, options)

    if '3_primary_metabolic_model'in folder:
        logging.info("'Primary' metabolic model completed")
    elif '4_complete_model' in folder:
        logging.info("'Secondary' metabolic model completed")

    if options.pmr_generation and options.debug:
        write_data_for_debug(options)


def get_model_reactions(folder, options, **kwargs):

    fp1 = open('./%s/model_reactions.txt' %folder, 'w')
    fp2 = open('./%s/rmc_remaining_essential_reactions_from_template_model.txt' %folder, 'w')
    fp3 = open('./%s/rmc_reactions_added_from_kegg.txt' %folder, 'w')

    fp1.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'
            +'GPR'+'\t'+'pathway'+'\n')
    fp2.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'
            +'GPR'+'\t'+'pathway'+'\n')
    fp3.write('reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'
            +'GPR'+'\t'+'pathway'+'\n')

    if '4_complete_model' in folder:
        if options.anti == 5:
            fp4 = open('./%s/rmc_region_fluxes.txt' %folder, 'w')
        if options.anti == 4:
            fp4 = open('./%s/rmc_cluster_fluxes.txt' %folder, 'w')
        fp4.write('reaction_ID'+'\t'+'fluxes without gap-filling reactions'+'\n')

    if 'cobra_model' in kwargs:
        cobra_model = kwargs['cobra_model']

    num_essen_rxn = 0
    num_kegg_rxn = 0
    num_secondary_rxn = 0
    for j in range(len(cobra_model.reactions)):
        rxn = cobra_model.reactions[j]
        print >>fp1, '%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem)

        #Remaining essential reactions
        if 'rxnToRemove_dict' in options:
            if rxn.id in options.rxnToRemove_dict.keys():
                if options.rxnToRemove_dict[rxn.id] == '0':
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
        if (re.search('Ex_Region', rxn.id) or re.search('Ex_Cluster', rxn.id)) and '4_complete_model' in folder:
            num_secondary_rxn +=1

            #Calculated flux values are inaccurate without
            #manual setting of objective_coefficient to zero
            obj_rxn = linear_reaction_coefficients(cobra_model).keys()[0].id
            cobra_model.reactions.get_by_id(obj_rxn).objective_coefficient = 0
            cobra_model.reactions.get_by_id(rxn.id).objective_coefficient = 1
            flux_dist = cobra_model.optimize()

            print >>fp4, '%s\t%f' %(rxn.id, flux_dist.objective_value)

            # NOTE: Currently disabled due to no gapfilling procedure at the moment
            #if 'cobra_model_no_gapFilled' in kwargs:
            #    cobra_model_no_gapFilled = kwargs['cobra_model_no_gapFilled']

            #    obj_rxn = linear_reaction_coefficients(cobra_model).keys()[0].id
            #    cobra_model_no_gapFilled.reactions.get_by_id(obj_rxn).objective_coefficient = 0
            #    cobra_model_no_gapFilled.reactions.get_by_id(rxn.id).objective_coefficient = 1
            #    flux_dist2 = cobra_model_no_gapFilled.optimize()

            #    print >>fp4, '%s\t%f\t%f' \
            #    %(rxn.id, flux_dist2.objective_value, flux_dist.objective_value)

    fp1.close()
    fp2.close()
    fp3.close()

    if '4_complete_model' in folder:
        fp4.close()

    return num_essen_rxn, num_kegg_rxn, num_secondary_rxn


def get_model_metabolites(folder, cobra_model, options):

    fp1 = open('./%s/model_metabolites.txt' %folder, "w")
    fp1.write('metabolite_ID'+'\t'+'metabolite_name'+'\t'
            +'formula'+'\t'+'compartment'+'\n')

    if '4_complete_model' in folder:
        fp2 = open('./%s/rmc_metabolites_gapfilling_needed.txt' %folder, "w")
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

    if '4_complete_model' in folder:
        fp2.close()


def get_model_genes(folder, cobra_model, options):

    template_model_gene_list = []
    duplicate_gene_list = []

    fp1 = open('./%s/rmc_gpr_associations_from_homology_analysis.txt' %folder, "w")

    fp1.write('gene'+'\t'+'note'+'\t'+'reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')

    if options.orgName == 'bsu':
        locustag_pattern = 'BSU[0-9]+[0-9]+[0-9]+[0-9]'
    elif options.orgName == 'cre':
        locustag_pattern = 'Cre'
    elif options.orgName == 'eco':
        locustag_pattern = 'b[0-9]+[0-9]+[0-9]+[0-9]'
    elif options.orgName == 'mtu':
        locustag_pattern = 'Rv[0-9]+[0-9]+[0-9]+[0-9]'
    elif options.orgName == 'nsal':
        locustag_pattern = 'NSV'
    elif options.orgName == 'ppu':
        locustag_pattern = 'PP_[0-9]+[0-9]+[0-9]+[0-9]'
    elif options.orgName == 'sco':
        locustag_pattern = 'SCO[0-9]+[0-9]+[0-9]+[0-9]'

    for g in range(len(cobra_model.genes)):
        gene = cobra_model.genes[g]

        if re.search(locustag_pattern, gene.id):
            if gene.id not in template_model_gene_list:
                template_model_gene_list.append(gene.id)

            for j in range(len(cobra_model.reactions)):
                rxn = cobra_model.reactions[j]

                if gene.id in rxn.gene_reaction_rule:
                    print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(
                                                    gene.id,
                                                    'remaining_gene_from_template_model',
                                                    rxn.id,
                                                    rxn.name,
                                                    rxn.reaction,
                                                    rxn.gene_reaction_rule,
                                                    rxn.subsystem)

        for j in range(len(cobra_model.reactions)):
            rxn = cobra_model.reactions[j]

            if len(re.findall(gene.id, rxn.gene_reaction_rule)) > 1:
                if gene.id not in duplicate_gene_list:
                    duplicate_gene_list.append(gene.id)
                    print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(
                                                    gene.id,
                                                    'duplicate_gene_from_target_model',
                                                    rxn.id,
                                                    rxn.name,
                                                    rxn.reaction,
                                                    rxn.gene_reaction_rule,
                                                    rxn.subsystem)

    fp1.close()
    return template_model_gene_list, duplicate_gene_list


def get_summary_report(folder, cobra_model, runtime,
                       num_essen_rxn, num_kegg_rxn, num_secondary_rxn,
                       template_model_gene_list, duplicate_gene_list,
                       options):
    fp1 = open('./%s/summary_report.txt' %folder, "w")

    #log_level
    if options.verbose:
        log_level = 'verbose'
    elif options.debug:
        log_level = 'debug'

    #runtime
    runtime2 = runtime.split()[2]

    model_summary_dict = {}
    model_summary_dict['number_cpu_use']=options.cpus
    model_summary_dict['input_file']=options.input
    model_summary_dict['outputfolder']=options.outputfolder
    model_summary_dict['template_model_organism']= options.orgName
    model_summary_dict['eficaz']=options.eficaz
    model_summary_dict['primary_metabolic_modeling']=options.pmr_generation
    model_summary_dict['secondary_metabolic_modeling']=options.smr_generation
    model_summary_dict['eficaz_file'] = options.eficaz_file
    model_summary_dict['compartment_file'] = options.comp
    model_summary_dict['log_level']=log_level
    model_summary_dict['program version']= 'GEMS version %s (%s)'\
                                            %(utils.get_version(), utils.get_git_log())
    model_summary_dict['number_genes']=len(cobra_model.genes)
    model_summary_dict['number_reactions']=len(cobra_model.reactions)
    model_summary_dict['number_metabolites']=len(cobra_model.metabolites)
    model_summary_dict['number_remaining_essential_reactions_from_template_model'] = \
            num_essen_rxn
    model_summary_dict['number_reactions_added_from_kegg']=num_kegg_rxn
    model_summary_dict['number_secondary_reactions']=num_secondary_rxn
    model_summary_dict['number_remaining_genes_from_template_model'] = \
            len(template_model_gene_list)
    model_summary_dict['number_duplicate_genes_in_rxn_from_target_model'] = \
            len(duplicate_gene_list)
    model_summary_dict['runtime']=runtime2

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

    with open('./%s/temp_target_BBH_dict.txt' %options.outputfolder2,'w') as f:
        for locustag in options.temp_target_BBH_dict:
            print >> f, '%s\t%s' %(locustag, options.temp_target_BBH_dict[locustag])

    try:
        with open('./%s/mnxr_to_add_list.txt' %options.outputfolder6,'w') as f:
            for mnxr in options.mnxr_to_add_list:
                print >>f, '%s' %mnxr
    except AttributeError, e:
        logging.warning(e)

    try:
        with open('./%s/targetGenome_locusTag_ec_nonBBH_dict.txt' %options.outputfolder6,'w') as f:
            for rxnid in options.targetGenome_locusTag_ec_nonBBH_dict:
                print >>f, '%s\t%s' %(rxnid, options.targetGenome_locusTag_ec_nonBBH_dict[rxnid])
    except AttributeError, e:
        logging.warning(e)

    try:
        with open('./%s/rxnid_info_dict.txt' %options.outputfolder6,'w') as f:
            for rxnid in options.rxnid_info_dict:
                print >>f, '%s\t%s' %(rxnid, options.rxnid_info_dict[rxnid])
    except AttributeError, e:
        logging.warning(e)

    try:
        with open('./%s/rxnid_locusTag_dict.txt' %options.outputfolder6,'w') as f:
            for rxnid in options.rxnid_locusTag_dict:
                print >>f, '%s\t%s' %(rxnid, options.rxnid_locusTag_dict[rxnid])
    except AttributeError, e:
        logging.warning(e)

    if options.comp:
        try:
            with open('./%s/rxn_newComp_fate_dict.txt' %options.outputfolder6,'w') as f:
                for rxnid in options.rxn_newComp_fate_dict:
                    print >>f, '%s\t%s' %(rxnid, options.rxn_newComp_fate_dict[rxnid])
        except AttributeError, e:
            logging.warning(e)


def remove_tmp_model_files(options):
    shutil.rmtree(options.outputfolder5)
    logging.info("'tmp_model_files' removed")
