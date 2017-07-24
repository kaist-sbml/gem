
import logging
import os
import pickle
from antismash_monomer_info import get_std_id_from_antismash_id
from cobra import Reaction, Metabolite
from gmsm import utils

def get_cluster_location(seq_record, cluster_nr, options):

    for feature in seq_record.features:

        if feature.type == 'cluster':
            cluster_number = 'Cluster number: %s' %cluster_nr
            options.cluster_number = cluster_number

            if options.cluster_number in feature.qualifiers.get('note'):
                options.cluster_loc1 = feature.location.start
                options.cluster_loc2 = feature.location.end


#Exract all the information associated with a particular locus_tag for the selected cluster
def get_cluster_info_from_seq_record(seq_record, options):

    cluster_info_dict = {}

    for feature in seq_record.features:

        if feature.type == 'CDS':
            if feature.location.start >= options.cluster_loc1 \
                    and feature.location.end <= options.cluster_loc2:
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                if feature.qualifiers.get('sec_met'):
                    cluster_info_dict[qualifier_locus_tag] = \
                            feature.qualifiers.get('sec_met')

    #logging.debug('cluster_info_dict: %s' %cluster_info_dict)
    options.cluster_info_dict = cluster_info_dict


def get_cluster_product(seq_record, cluster_nr, options):

    for feature in seq_record.features:

        #Retrieving "Cluster number"
        if feature.type == 'cluster':
            qualifier_cluster = feature.qualifiers.get('note')
            if options.cluster_number in qualifier_cluster:
                product = feature.qualifiers.get('product')[0]

    #Handle legacy problem
    product2 = product.replace("-","_")

    if float(cluster_nr) < 10:
        product3 = "Cluster0"+str(cluster_nr)+"_"+product2
    else:
        product3 = "Cluster"+str(cluster_nr)+"_"+product2

    options.product = product3


#Output: e.g., {'SAV_943_M1':['mmal', 'Ethyl_mal', 'pk']}
def get_cluster_monomers(options):

    locustag_monomer_dict = {}
    for each_gene in options.cluster_info_dict.keys():
        module_count = 0

        for sec_met_info in options.cluster_info_dict[each_gene]:

            if 'Substrate specificity predictions' in sec_met_info:
                #type(sec_met_info) is string
                #Convert 'sec_met_info' into a list
                sec_met_info_list = sec_met_info.split(';')
                for each_sec_met_info in sec_met_info_list:
                    if 'Substrate specificity predictions' in each_sec_met_info:
                        pred_monomer_list = []
                        #Following statement produces a list with 2 elements:
                        #e.g., [' Substrate specificity predictions',
                        #' gly (NRPSPredictor2 SVM), gly (Stachelhaus code),
                        #gly (Minowa), gly (consensus)']
                        monomer_list = each_sec_met_info.split(':')
                        for monomer in monomer_list:
                            #Make sure to include a space ''
                            #for monomers predicted from different engines:
                            #e.g., 'orn,lys,arg (NRPSPredictor2 SVM),
                            #lys (Stachelhaus code), leu (Minowa), nrp (consensus)'
                            if 'Substrate specificity predictions' not in monomer \
                                    and ', ' in monomer:
                                monomer_list2 = monomer.split(', ')
                                for monomer2 in monomer_list2:
                                    monomer2_list = monomer2.split('(')
                                    monomer3 = monomer2_list[0].strip()
                                    pred_monomer_list.append(monomer3)

                module_number = each_gene + '_M' + str(module_count)
                locustag_monomer_dict[module_number] = pred_monomer_list
                module_count += 1

    options.locustag_monomer_dict = locustag_monomer_dict


# Add stoichiometric coeff's of monomers
# Output: e.g., {'mmalcoa': -4, 'malcoa': -7}
def get_all_metab_coeff(options):

    metab_coeff_dict = {}
    for each_module in options.locustag_monomer_dict:

        # For gbk output from antiSMASH 3.0
        #locustag_monomer_dict[each_module] for nrps
        #Position [0]: NRPSPredictor2 SVM
        #Position [1]: Stachelhaus code
        #Position [2]: Minowa
        #Position [3]: consensus
        if len(options.locustag_monomer_dict[each_module]) == 4:

            # In case "consensus" is not reached:
            if options.locustag_monomer_dict[each_module][3] == 'nrp':

                # In case of non-consensus, NRPSPredictor2 SVM is first considered
                # Not considered: e.g., NRPSPredictor2 SVM: val,leu,ile,abu,iva
                aSid_met = options.locustag_monomer_dict[each_module][0]
                biggid_met = get_std_id_from_antismash_id(aSid_met)

                if not biggid_met:
                    # Minowa is next considered
                    aSid_met = options.locustag_monomer_dict[each_module][2]
                    biggid_met = get_std_id_from_antismash_id(aSid_met)

                    if not biggid_met:
                        # Stachelhaus code is next considered
                        aSid_met = options.locustag_monomer_dict[each_module][1]
                        biggid_met = get_std_id_from_antismash_id(aSid_met)

            # In case "consensus" is reached:
            elif options.locustag_monomer_dict[each_module][3] != 'nrp':
                aSid_met = options.locustag_monomer_dict[each_module][3]
                biggid_met = get_std_id_from_antismash_id(aSid_met)

        # For gbk output from antiSMASH 4.0
        #locustag_monomer_dict[each_module] for nrps
        #Position [0]: Stachelhaus code
        #Position [1]: NRPSPredictor3 SVM
        #Position [2]: pHMM
        #Position [3]: PrediCAT
        #Position [4]: SANDPUMA ensemble
        if len(options.locustag_monomer_dict[each_module]) == 5:

            # In case "consensus" (SANDPUMA ensemble) is not reached:
            if options.locustag_monomer_dict[each_module][4] == 'no_call':

                # In case of non-consensus, NRPSPredictor3 SVM is first considered
                # Not considered: e.g., NRPSPredictor2 SVM: val,leu,ile,abu,iva
                aSid_met = options.locustag_monomer_dict[each_module][1]
                biggid_met = get_std_id_from_antismash_id(aSid_met)

                # Stachelhaus code is considered
                if not biggid_met:
                    aSid_met = options.locustag_monomer_dict[each_module][0]
                    biggid_met = get_std_id_from_antismash_id(aSid_met)

                    # If Stachelhaus code has invalid monomer, then TODO is considered
                    if not biggid_met:
                        aSid_met = options.locustag_monomer_dict[each_module][1]
                        biggid_met = get_std_id_from_antismash_id(aSid_met)

            # In case "consensus" (SANDPUMA ensemble) is reached:
            elif options.locustag_monomer_dict[each_module][4] != 'no_call':
                aSid_met = options.locustag_monomer_dict[each_module][4]
                biggid_met = get_std_id_from_antismash_id(aSid_met)

        # For gbk output from antiSMASH 3.0 & 4.0
        #locustag_monomer_dict[each_module] for pks
        #Position [0]: PKS signature
        #Position [1]: Minowa
        #Position [2]: consensus
        elif len(options.locustag_monomer_dict[each_module]) == 3:

            if len(options.locustag_monomer_dict[each_module]) < 3:
                continue

            # In case "consensus" is not reached:
            if options.locustag_monomer_dict[each_module][2] == 'pk':

                # From PKS signature
                # In case of non-consensus, PKS signature is considered
                aSid_met = options.locustag_monomer_dict[each_module][0]
                if aSid_met != 'N/A' and aSid_met != 'mal_or_prop':
                    biggid_met = get_std_id_from_antismash_id(aSid_met)

                # If PKS signature has invalid monomer, then Minowa is considered
                else:
                    aSid_met = options.locustag_monomer_dict[each_module][1]
                    if aSid_met != 'inactive':
                        biggid_met = get_std_id_from_antismash_id(aSid_met)

            # In case "consensus" is reached:
            elif options.locustag_monomer_dict[each_module][2] != 'pk':
                aSid_met = options.locustag_monomer_dict[each_module][2]

                # Original monomers are considered for those reduced by KR, DH and/or ER
                if aSid_met == 'ohmal' or aSid_met == 'ccmal' or aSid_met == 'redmal':
                    aSid_met = 'mal'
                elif aSid_met == 'ohmmal' or aSid_met == 'ccmmal' \
                        or aSid_met == 'redmmal':
                    aSid_met = 'mmal'
                elif aSid_met == 'ohmxmal' or aSid_met == 'ccmxmal' \
                        or aSid_met == 'redmxmal':
                    aSid_met = 'mxmal'
                elif aSid_met == 'ohemal' or aSid_met == 'ccemal' \
                        or aSid_met == 'redemal':
                    aSid_met = 'emal'

                biggid_met = get_std_id_from_antismash_id(aSid_met)

        if options.locustag_monomer_dict[each_module]:
            if biggid_met not in metab_coeff_dict:
                metab_coeff_dict[biggid_met] = -1
            else:
                metab_coeff_dict[biggid_met] -= 1

    # Add secondary metabolite product to the reaction
    metab_coeff_dict[options.product] = 1

    logging.debug('metab_coeff_dict: %s' %metab_coeff_dict)
    options.metab_coeff_dict = metab_coeff_dict


def get_pickles(options):

    if not hasattr(options, 'mnxref'):
        mnxref = pickle.load(open('./gmsm/io/data/input2/MNXref.p','rb'))
        options.mnxref = mnxref

    if not hasattr(options, 'mnxm_compoundInfo_dict'):
        mnxm_compoundInfo_dict = pickle.load(
                open('./gmsm/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
        options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict


def add_sec_met_rxn(target_model, options):

    #ID
    rxn = Reaction(options.product)

    #Reversibility / Lower and upper bounds
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Metabolites and their stoichiometric coeff's
    for metab in options.metab_coeff_dict.keys():

        #Consider only metabolites consumed or produced
        if options.metab_coeff_dict[metab] != 0:
            metab_compt = '_'.join([metab,'c'])

            #Add  metabolites already in the model
            if metab_compt in target_model.metabolites:
                logging.debug('Metabolite %s already present in the model', metab_compt)
                rxn.add_metabolites({target_model.metabolites.get_by_id(
                    metab_compt):options.metab_coeff_dict[metab]})

            #Add metabolites available in the MNXref sbml
            elif metab_compt in options.mnxref.metabolites:
                logging.debug('Metabolite %s available in the MNXref sbml', metab_compt)
                metab_compt = options.mnxref.metabolites.get_by_id(metab_compt)
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

            elif 'Cluster' in metab:
                logging.debug("Secondary metabolite ('Cluster') %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

            #Add metabolites not available in the MNXref sbml
            else:
                logging.debug("Metabolite (MNXM ID) %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt,
                        formula = options.mnxm_compoundInfo_dict[metab][1],
                        name = options.mnxm_compoundInfo_dict[metab][0],
                        compartment='c')
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

    #GPR association
    gpr_count = 0
    for each_gene in options.cluster_info_dict.keys():
        if gpr_count == 0:
            gpr_list = each_gene
            gpr_count += 1
        else:
            gpr_list = ' and '.join([gpr_list, each_gene])

    gpr_list = '( %s )' %(gpr_list)

    rxn.gene_reaction_rule = gpr_list

    #Add a new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("Secondary metabolite reaction: %s" %rxn)
    logging.debug(rxn.reaction)

    ##############################
    #Create a transport reaction
    #Create reaction ID
    rxn = Reaction("Transport_" + options.product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Add a substrate metabolite
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(options.product+'_c')):-1})

    #Add product metabolite(s)
    product_e = Metabolite(options.product+"_e", name='', compartment='e')
    rxn.add_metabolites({product_e:1})

    #Add the new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("Transport reaction: %s" %rxn)
    logging.debug(rxn.reaction)

    ##############################
    #Create an exchange reaction
    #Create reaction ID
    rxn = Reaction("Ex_"+options.product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible 0: irreversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Add a substrate metabolite
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(product_e)):-1})

    #Add a new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("Exchange reaction: %s" %rxn)
    logging.debug(rxn.reaction)

    return target_model


def check_producibility_sec_met(target_model, options):
    for rxn in target_model.reactions:
        rxn.objective_coefficient = 0

    #target_model.reactions.get_by_id('Biomass_SCO').objective_coefficient = 0
    target_model.reactions.get_by_id("Ex_"+options.product).objective_coefficient = 1

    #Model reloading and overwrtting are necessary for model stability
    #Without these, model does not produce an accurate prediction
    target_model = utils.stabilize_model(
            target_model, options.outputfolder5, options.product)

    flux_dist = target_model.optimize()
    logging.debug("Flux: %s" %flux_dist.objective_value)

    target_model.reactions.get_by_id('Biomass_SCO').objective_coefficient = 1
    target_model.reactions.get_by_id("Ex_"+options.product).objective_coefficient = 0

    return target_model, flux_dist


def get_sec_met_monomers(sec_met_list, options):

    for metab in options.metab_coeff_dict.keys():
        if options.metab_coeff_dict[metab] <0:
            if metab not in sec_met_list:
                sec_met_list.append(metab)

    return sec_met_list

