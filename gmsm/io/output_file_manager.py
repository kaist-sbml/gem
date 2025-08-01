
import cobra
import collections
import logging
import pickle
import re
import shutil
from cobra.util.solver import linear_reaction_coefficients
from gmsm import utils

def generate_outputs(folder, runtime, run_ns, io_ns, homology_ns, primary_model_ns, secondary_model_ns, **kwargs):
    if 'cobra_model' in kwargs:
        cobra_model = kwargs['cobra_model']

    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    cobra_model = utils.stabilize_model(cobra_model, folder, '')

    num_essen_rxn, num_kegg_rxn, num_bgc_rxn = get_model_reactions(
                       folder, primary_model_ns, **kwargs)
    get_model_metabolites(folder, cobra_model, secondary_model_ns)

    template_model_gene_list, duplicate_gene_list = \
                       get_model_genes(folder, cobra_model, run_ns)
    get_summary_report(folder, cobra_model, runtime,
                       num_essen_rxn, num_kegg_rxn, num_bgc_rxn,
                       template_model_gene_list, duplicate_gene_list, run_ns, secondary_model_ns)

    if '3_primary_metabolic_model'in folder:
        logging.info("'Primary' metabolic model completed")
    elif '4_complete_model' in folder:
        logging.info("'Secondary' metabolic model completed")

    if run_ns.pmr_generation and run_ns.debug:
        write_data_for_debug(run_ns, io_ns, homology_ns, primary_model_ns)


def get_model_reactions(folder, primary_model_ns, **kwargs):

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
        fp4 = open('./%s/rmc_BGCs_fluxes.txt' %folder, 'w')
        fp4.write('reaction_ID'+'\t'+'fluxes without gap-filling reactions'+'\n')

    if 'cobra_model' in kwargs:
        cobra_model = kwargs['cobra_model']

    num_essen_rxn = 0
    num_kegg_rxn = 0
    num_bgc_rxn = 0
    for j in range(len(cobra_model.reactions)):
        rxn = cobra_model.reactions[j]
        print('%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem), file=fp1)

        #Remaining essential reactions
        if 'rxnToRemove_dict' in primary_model_ns:
            if rxn.id in primary_model_ns.rxnToRemove_dict.keys():
                if primary_model_ns.rxnToRemove_dict[rxn.id] == '0':
                    num_essen_rxn+=1
                    print('%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem), file=fp2)
        else:
            print('Primary metabolic modeling not performed', file=fp2)
            num_essen_rxn = 'Primary metabolic modeling not performed'

        #Reactions added from KEGG
        if re.search('R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]', rxn.id):
            num_kegg_rxn+=1
            print('%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name, rxn.reaction,
                                            rxn.gene_reaction_rule, rxn.subsystem), file=fp3)

        #Secondary metabolite biosynthetic reactions
        if (re.search('EX_', rxn.id) and re.search('(Region|Cluster)', rxn.id)) and '4_complete_model' in folder:
            num_bgc_rxn+=1

            #Calculated flux values are inaccurate without
            #manual setting of objective_coefficient to zero
            obj_rxn = list(linear_reaction_coefficients(cobra_model).keys())[0].id
            cobra_model.reactions.get_by_id(obj_rxn).objective_coefficient = 0
            cobra_model.reactions.get_by_id(rxn.id).objective_coefficient = 1
            flux_dist = cobra_model.optimize()

            print('%s\t%f' %(rxn.id, flux_dist.objective_value), file=fp4)

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

    return num_essen_rxn, num_kegg_rxn, num_bgc_rxn


def get_model_metabolites(folder, cobra_model, secondary_model_ns):

    fp1 = open('./%s/model_metabolites.txt' %folder, "w")
    fp1.write('metabolite_ID'+'\t'+'metabolite_name'+'\t'
            +'formula'+'\t'+'compartment'+'\n')

    if '4_complete_model' in folder:
        fp2 = open('./%s/rmc_metabolites_gapfilling_needed.txt' %folder, "w")
        fp2.write('metabolite_ID'+'\t'+'reaction_ID'+'\t'+'reaction_name'+'\t'
                +'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')

    for i in range(len(cobra_model.metabolites)):
        metab = cobra_model.metabolites[i]
        print('%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula,
                metab.compartment), file=fp1)

        if '4_complete_model' in folder:
            #Remove compartment suffix (e.g., '_c') from 'metab.id'
            if metab.id[:-2] in secondary_model_ns.adj_unique_nonprod_monomers_list:
                logging.debug("Metabolite for gap-filling: %s" %metab.id)

                for j in range(len(cobra_model.reactions)):
                    rxn = cobra_model.reactions[j]

                    if metab.id in rxn.reaction:
                        logging.debug("Relevant reactions: %s" %rxn.id)
                        print('%s\t%s\t%s\t%s\t%s\t%s' %(metab.id,
                            rxn.id, rxn.name, rxn.reaction,
                            rxn.gene_reaction_rule, rxn.subsystem), file=fp2)

    fp1.close()

    if '4_complete_model' in folder:
        fp2.close()


def get_model_genes(folder, cobra_model, run_ns):
    template_model_gene_list = []
    duplicate_gene_list = []

    fp1 = open('./%s/rmc_gpr_associations_from_homology_analysis.txt' %folder, "w")

    fp1.write('gene'+'\t'+'note'+'\t'+'reaction_ID'+'\t'+'reaction_name'+'\t'+'reaction_equation'+'\t'+'GPR'+'\t'+'pathway'+'\n')

    if run_ns.orgName == 'bsu':
        locustag_pattern = 'BSU[0-9]+[0-9]+[0-9]+[0-9]'
    elif run_ns.orgName == 'clj':
        locustag_pattern = 'CLJU_RS[0-9]+[0-9]+[0-9]+[0-9]+[0-9]'
    elif run_ns.orgName == 'cre':
        locustag_pattern = 'Cre'
    elif run_ns.orgName == 'eco':
        locustag_pattern = 'b[0-9]+[0-9]+[0-9]+[0-9]'
    elif run_ns.orgName == 'hpy':
        locustag_pattern = 'HP[0-9]+[0-9]+[0-9]+[0-9]'
    elif run_ns.orgName == 'mtu':
        locustag_pattern = 'Rv[0-9]+[0-9]+[0-9]+[0-9]'
    elif run_ns.orgName == 'nsal':
        locustag_pattern = 'NSV'
    elif run_ns.orgName == 'ppu':
        locustag_pattern = 'PP_[0-9]+[0-9]+[0-9]+[0-9]'
    elif run_ns.orgName == 'sce':
        locustag_pattern = 'Y[A-Z]+[A-Z]+[0-9]+[0-9]+[0-9]+[A-Z]'
    elif run_ns.orgName == 'sco':
        locustag_pattern = 'SCO[0-9]+[0-9]+[0-9]+[0-9]'

    for g in range(len(cobra_model.genes)):
        gene = cobra_model.genes[g]

        if re.search(locustag_pattern, gene.id):
            if gene.id not in template_model_gene_list:
                template_model_gene_list.append(gene.id)

            for j in range(len(cobra_model.reactions)):
                rxn = cobra_model.reactions[j]

                if gene.id in rxn.gene_reaction_rule:
                    print('%s\t%s\t%s\t%s\t%s\t%s\t%s' %(
                                                    gene.id,
                                                    'remaining_gene_from_template_model',
                                                    rxn.id,
                                                    rxn.name,
                                                    rxn.reaction,
                                                    rxn.gene_reaction_rule,
                                                    rxn.subsystem), file=fp1)

        for j in range(len(cobra_model.reactions)):
            rxn = cobra_model.reactions[j]

            if len(re.findall(gene.id, rxn.gene_reaction_rule)) > 1:
                if gene.id not in duplicate_gene_list:
                    duplicate_gene_list.append(gene.id)
                    print('%s\t%s\t%s\t%s\t%s\t%s\t%s' %(
                                                    gene.id,
                                                    'duplicate_gene_from_target_model',
                                                    rxn.id,
                                                    rxn.name,
                                                    rxn.reaction,
                                                    rxn.gene_reaction_rule,
                                                    rxn.subsystem), file=fp1)

    fp1.close()
    return template_model_gene_list, duplicate_gene_list


def get_summary_report(folder, cobra_model, runtime,
                       num_essen_rxn, num_kegg_rxn, num_bgc_rxn,
                       template_model_gene_list, duplicate_gene_list, run_ns, secondary_model_ns):

    fp1 = open('./%s/summary_report.txt' %folder, "w")

    #log_level
    if run_ns.verbose:
        log_level = 'verbose'
    elif run_ns.debug:
        log_level = 'debug'

    #runtime
    runtime2 = runtime.split()[2]

    model_summary_dict = {}
    model_summary_dict['number_cpu_use']=run_ns.cpus
    model_summary_dict['input_file']=run_ns.input
    model_summary_dict['outputfolder']=run_ns.outputfolder
    model_summary_dict['template_model_organism']=run_ns.orgName
    model_summary_dict['eficaz']=run_ns.eficaz
    model_summary_dict['primary_metabolic_modeling']=run_ns.pmr_generation
    model_summary_dict['secondary_metabolic_modeling']=run_ns.smr_generation
    model_summary_dict['EC_number_file']=run_ns.ec_file
    model_summary_dict['compartment_file']=run_ns.comp
    model_summary_dict['log_level']=log_level
    model_summary_dict['program version']='GMSM version %s (%s)'\
                                            %(utils.get_version(), utils.get_git_log())
    model_summary_dict['number_genes']=len(cobra_model.genes)
    model_summary_dict['number_reactions']=len(cobra_model.reactions)
    model_summary_dict['number_metabolites']=len(cobra_model.metabolites)
    model_summary_dict['number_remaining_essential_reactions_from_template_model'] = \
            num_essen_rxn
    model_summary_dict['number_reactions_added_from_kegg']=num_kegg_rxn
    model_summary_dict['number_BGCs_for_reactions']=num_bgc_rxn
    model_summary_dict['number_remaining_genes_from_template_model'] = \
            len(template_model_gene_list)
    model_summary_dict['number_duplicate_genes_in_rxn_from_target_model'] = \
            len(duplicate_gene_list)
    model_summary_dict['runtime']=runtime2

    if '4_complete_model' in folder:
        model_summary_dict['number_metabolites_for_gapfilling'] \
            =len(secondary_model_ns.adj_unique_nonprod_monomers_list)
    else:
        model_summary_dict['number_metabolites_for_gapfilling']=0

    #Sort data by keys
    model_summary_dict2 = collections.OrderedDict(sorted(model_summary_dict.items()))

    for key in model_summary_dict2.keys():
        print('%s\t%s' %(key, model_summary_dict2[key]), file=fp1)

    fp1.close()


def write_data_for_debug(run_ns, io_ns, homology_ns, primary_model_ns):

    with open('./%s/temp_target_BBH_dict.txt' %io_ns.outputfolder2,'w') as f:
        for locustag in homology_ns.temp_target_BBH_dict:
            print('%s\t%s' %(locustag, homology_ns.temp_target_BBH_dict[locustag]), file=f)

    try:
        with open('./%s/mnxr_to_add_list.txt' %io_ns.outputfolder6,'w') as f:
            for mnxr in primary_model_ns.mnxr_to_add_list:
                print('%s' %mnxr, file=f)
    except AttributeError as e:
        logging.warning(e)

    try:
        with open('./%s/targetGenome_locusTag_ec_nonBBH_dict.txt' %io_ns.outputfolder6,'w') as f:
            for rxnid in primary_model_ns.targetGenome_locusTag_ec_nonBBH_dict:
                print('%s\t%s' %(rxnid, primary_model_ns.targetGenome_locusTag_ec_nonBBH_dict[rxnid]), file=f)
    except AttributeError as e:
        logging.warning(e)

    try:
        with open('./%s/rxnid_info_dict.txt' %io_ns.outputfolder6,'w') as f:
            for rxnid in primary_model_ns.rxnid_info_dict:
                print('%s\t%s' %(rxnid, primary_model_ns.rxnid_info_dict[rxnid]), file=f)
    except AttributeError as e:
        logging.warning(e)

    try:
        with open('./%s/rxnid_locusTag_dict.txt' %io_ns.outputfolder6,'w') as f:
            for rxnid in primary_model_ns.rxnid_locusTag_dict:
                print('%s\t%s' %(rxnid, primary_model_ns.rxnid_locusTag_dict[rxnid]), file=f)
    except AttributeError as e:
        logging.warning(e)

    if run_ns.comp:
        try:
            with open('./%s/rxn_newComp_fate_dict.txt' %io_ns.outputfolder6,'w') as f:
                for rxnid in primary_model_ns.rxn_newComp_fate_dict:
                    print('%s\t%s' %(rxnid, primary_model_ns.rxn_newComp_fate_dict[rxnid]), file=f)
        except AttributeError as e:
            logging.warning(e)


def remove_tmp_model_files(io_ns):
    shutil.rmtree(io_ns.outputfolder5)
    logging.info("'tmp_model_files' removed")
