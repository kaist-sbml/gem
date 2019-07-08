
import logging
import os
import pickle
from antismash_monomer_info import get_std_id_from_antismash_id
from cobra import Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from gmsm import utils


def get_region_location(seq_record, secondary_model_ns):

    for feature in seq_record.features:
        if feature.type == 'region':
            if feature.location.start > secondary_model_ns.temp_loc1:
                secondary_model_ns.region_loc1 = feature.location.start
                secondary_model_ns.region_loc2 = feature.location.end
                secondary_model_ns.temp_loc1 = feature.location.start
                secondary_model_ns.temp_loc2 = feature.location.end
                break

def get_cluster_location(seq_record, cluster_nr, secondary_model_ns):

    for feature in seq_record.features:
        if feature.type == 'cluster':
            cluster_number = 'Cluster number: %s' %cluster_nr
            secondary_model_ns.cluster_number = cluster_number

            if secondary_model_ns.cluster_number in feature.qualifiers.get('note'):
                secondary_model_ns.cluster_loc1 = feature.location.start
                secondary_model_ns.cluster_loc2 = feature.location.end


#Exract all the information associated with a particular locus_tag for the selected region
def get_region_info_from_seq_record(seq_record, region_nr, secondary_model_ns):

    region_info_dict = {}

    for feature in seq_record.features:
        if feature.type == 'CDS':
            if feature.location.start >= secondary_model_ns.region_loc1 \
                    and feature.location.end <= secondary_model_ns.region_loc2:
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                if feature.qualifiers.get('sec_met_domain'):
                    region_info_dict[qualifier_locus_tag] = \
                            feature.qualifiers.get('sec_met_domain')

    secondary_model_ns.region_info_dict = region_info_dict


#Exract all the information associated with a particular locus_tag for the selected cluster
def get_cluster_info_from_seq_record(seq_record, secondary_model_ns):

    cluster_info_dict = {}

    for feature in seq_record.features:
        if feature.type == 'CDS':
            if feature.location.start >= secondary_model_ns.cluster_loc1 \
                    and feature.location.end <= secondary_model_ns.cluster_loc2:
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                if feature.qualifiers.get('sec_met'):
                    cluster_info_dict[qualifier_locus_tag] = \
                            feature.qualifiers.get('sec_met')

    secondary_model_ns.cluster_info_dict = cluster_info_dict


def get_region_product(seq_record, region_nr, secondary_model_ns):

    for feature in seq_record.features:
        if feature.type == 'region':

            if feature.location.start == secondary_model_ns.region_loc1 and \
            feature.location.end == secondary_model_ns.region_loc2:
                product_list = feature.qualifiers.get('product')

                # connecting all product of a region
                for i in range(len(product_list)):
                    if i == 0:
                        #Handle legacy problem
                        product = product_list[i].replace('-','_')
                        if float(region_nr) < 10:
                            product2 = seq_record.id.replace('.','_')+"_Region0"+str(region_nr)+"_"+product.lower()
                        else:
                            product2 = seq_record.id.replace('.','_')+"_Region"+str(region_nr)+"_"+product.lower()
                    else:
                        product = product_list[i].replace('-','_')
                        product2 = product2+'_and_'+product.lower()

                    secondary_model_ns.product = product2


def get_cluster_product(seq_record, cluster_nr, secondary_model_ns):

    for feature in seq_record.features:

        #Retrieving "Cluster number"
        if feature.type == 'cluster':
            qualifier_cluster = feature.qualifiers.get('note')
            if secondary_model_ns.cluster_number in qualifier_cluster:
                product = feature.qualifiers.get('product')[0]

    #Handle legacy problem
    product2 = product.replace("-","_")

    if float(cluster_nr) < 10:
        product3 = seq_record.id.replace('.','_')+"_Cluster0"+str(cluster_nr)+"_"+product2
    else:
        product3 = seq_record.id.replace('.','_')+"_Cluster"+str(cluster_nr)+"_"+product2

    secondary_model_ns.product = product3


#Output: e.g., {'SAV_943_M1':['mmal', 'Ethyl_mal', 'pk']}
def get_region_monomers(seq_record, region_nr, secondary_model_ns):

    secondary_model_ns.cds_info_dict = {}
    secondary_model_ns.locustag_monomer_dict = {}

    for feature in seq_record.features: 
        # find location of CDS involved in SM
        if feature.type == 'CDS':
            module_count = 0
            secondary_model_ns.cds_loc1 = feature.location.start
            secondary_model_ns.cds_loc2 = feature.location.end

        # find antiSMASH domain in range of CDS
        if feature.type == 'aSDomain':
            if feature.location.start >= secondary_model_ns.cds_loc1 \
            and feature.location.end <= secondary_model_ns.cds_loc2:
                # collect monomer information in one region
                if feature.location.start >= secondary_model_ns.region_loc1 and \
                feature.location.end <= secondary_model_ns.region_loc2:
                    qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                    module_number = qualifier_locus_tag + '_M' + str(module_count)

                    if feature.qualifiers.get('specificity'):
                        sec_met_info = feature.qualifiers.get('specificity')
                        # exclude KR stochiometry information
                        if 'KR ' not in sec_met_info[0]:
                            secondary_model_ns.cds_info_dict[module_number] = \
                                feature.qualifiers.get('specificity')
                            module_count += 1

    for each_module in secondary_model_ns.cds_info_dict.keys():
        monomer_list = []

        for sec_met_info in secondary_model_ns.cds_info_dict[each_module]:
            #type(sec_met_info) is string
            #Convert 'sec_met_info' into a list
            sec_met = sec_met_info.split(': ')
            monomer = sec_met[1]
            monomer_list.append(monomer)

        secondary_model_ns.locustag_monomer_dict[each_module] = monomer_list


def get_cluster_monomers(secondary_model_ns):

    locustag_monomer_dict = {}
    for each_gene in secondary_model_ns.cluster_info_dict.keys():
        module_count = 0

        for sec_met_info in secondary_model_ns.cluster_info_dict[each_gene]:
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

    secondary_model_ns.locustag_monomer_dict = locustag_monomer_dict


def get_biggid(priority_list, each_module, secondary_model_ns):
    for i in priority_list:
        aSid_met = secondary_model_ns.locustag_monomer_dict[each_module][i]
        biggid_met = get_std_id_from_antismash_id(aSid_met)

        if biggid_met:
            logging.debug("Following monomer prediction taken: '%s'", aSid_met)
            break
        else:
            logging.debug("Following monomer prediction not considered: '%s'", aSid_met)

    if biggid_met:
        return biggid_met
    else:
        return


# Add stoichiometric coeff's of monomers
# Output: e.g., {'mmalcoa': -4, 'malcoa': -7}
def get_all_metab_coeff(io_ns, secondary_model_ns):
    metab_coeff_dict = {}
    for each_module in secondary_model_ns.locustag_monomer_dict:
        logging.debug("Module: %s; monomers: %s", each_module, secondary_model_ns.locustag_monomer_dict[each_module])

        # Paused (antiSMASH 3.0 is not used now)
        # NRPS analyzed with antiSMASH 3.0
        # locustag_monomer_dict[each_module]
        # Index [0]: NRPSPredictor2 SVM
        # Index [1]: Stachelhaus code
        # Index [2]: Minowa
        # Index [3]: consensus
        # Priority: consensus > NRPSPredictor2 SVM > Stachelhaus code > Minowa
        #if len(secondary_model_ns.locustag_monomer_dict[each_module]) == 4:
            #priority_list = [3, 0, 1, 2]
            #biggid_met = get_biggid(priority_list, each_module, secondary_model_ns)

        # NRPS analyzed with antiSMASH 4.0
        # Index [0]: Stachelhaus code
        # Index [1]: NRPSPredictor3 SVM
        # Index [2]: pHMM
        # Index [3]: PrediCAT
        # Index [4]: SANDPUMA ensemble
        # Priority: consensus (SANDPUMA ensemble) > NRPSPredictor3 SVM > PrediCAT >
        #pHMM > Stachelhaus code
        if len(secondary_model_ns.locustag_monomer_dict[each_module]) == 5:
            priority_list = [4, 1, 3, 2, 0]
            biggid_met = get_biggid(priority_list, each_module, secondary_model_ns)

        #NRPS analyzed with antiSMASH 5.0
        elif len(secondary_model_ns.locustag_monomer_dict[each_module]) == 1:
            #index [0]: NRPSpredictor2 SVM
            priority_list = [0]
            biggid_met = get_biggid(priority_list, each_module, secondary_model_ns)

        elif len(secondary_model_ns.locustag_monomer_dict[each_module]) == 3:
            if io_ns.anti_version == 4:
                # PKS_AT analyzed with antiSMASH 4.0 (also available for antiSMASH 3.0 but antiSMASH 3.0 is not used now)
                # Index [0]: PKS signature
                # Index [1]: Minowa
                # Index [2]: consensus
                # Priority: consensus > PKS signature > Minowa
                priority_list = [2, 0, 1]
                biggid_met = get_biggid(priority_list, each_module, secondary_model_ns)
            else:
                # PKS_AT analyzed with antiSMASH 5.0
                # locustag_monomer_dict[each_module] for pks
                # Index [0]: consensus
                # Index [1]: Minowa
                # Index [2]: PKS signature
                # Priority: consensus > PKS signature > Minowa
                priority_list = [0, 2, 1]
                biggid_met = get_biggid(priority_list, each_module, secondary_model_ns)

        # Filter modules without 'Substrate specificity predictions'
        # e.g., 'B446_13275_M0'
        # {'B446_13415_M0': ['gly,ala,val,leu,ile,abu,iva', 'ile', 'val', 'nrp'], 'B446_13275_M0': [], 'B446_13350_M0': ['leu', 'ala', 'sar', 'nrp'], 'B446_13445_M0': ['pro,pip', 'N/A', 'orn', 'nrp']}
        if secondary_model_ns.locustag_monomer_dict[each_module] and biggid_met:
            if biggid_met not in metab_coeff_dict:
                metab_coeff_dict[biggid_met] = -1
            else:
                metab_coeff_dict[biggid_met] -= 1

    # Add secondary metabolite product to the reaction
    metab_coeff_dict[secondary_model_ns.product] = 1
    logging.debug('metab_coeff_dict: %s' %metab_coeff_dict)
    secondary_model_ns.metab_coeff_dict = metab_coeff_dict


def get_pickles(io_ns):

    if not hasattr(io_ns, 'mnxref'):
        mnxref = pickle.load(open('./gmsm/io/data/input2/MNXref.p','rb'))
        io_ns.mnxref = mnxref

    if not hasattr(io_ns, 'mnxm_compoundInfo_dict'):
        mnxm_compoundInfo_dict = pickle.load(
                open('./gmsm/io/data/input2/mnxm_compoundInfo_dict.p','rb'))
        io_ns.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict


def add_sec_met_rxn(target_model, io_ns, secondary_model_ns):

    #ID
    rxn = Reaction(secondary_model_ns.product)

    #Reversibility / Lower and upper bounds
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Metabolites and their stoichiometric coeff's
    for metab in secondary_model_ns.metab_coeff_dict.keys():

        #Consider only metabolites consumed or produced
        if secondary_model_ns.metab_coeff_dict[metab] != 0:
            metab_compt = '_'.join([metab,'c'])

            #Add  metabolites already in the model
            if metab_compt in target_model.metabolites:
                logging.debug('Metabolite %s already present in the model', metab_compt)
                rxn.add_metabolites({target_model.metabolites.get_by_id(
                    metab_compt):secondary_model_ns.metab_coeff_dict[metab]})

            #Add metabolites available in the MNXref sbml
            elif metab_compt in io_ns.mnxref.metabolites:
                logging.debug('Metabolite %s available in the MNXref sbml', metab_compt)
                metab_compt = io_ns.mnxref.metabolites.get_by_id(metab_compt)
                rxn.add_metabolites({metab_compt:secondary_model_ns.metab_coeff_dict[metab]})

            elif 'Region' in metab:
                logging.debug("Secondary metabolite ('Region') %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:secondary_model_ns.metab_coeff_dict[metab]})

            elif 'Cluster' in metab:
                logging.debug("Secondary metabolite ('Cluster') %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:secondary_model_ns.metab_coeff_dict[metab]})

            #Add metabolites not available in the MNXref sbml
            else:
                logging.debug("Metabolite (MNXM ID) %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:secondary_model_ns.metab_coeff_dict[metab]})

    #GPR association
    gpr_count = 0
    if io_ns.anti_version == 5:
        for each_gene in secondary_model_ns.region_info_dict.keys():
            if gpr_count == 0:
                gpr_list = each_gene
                gpr_count += 1
            else:
                gpr_list = ' and '.join([gpr_list, each_gene])

    elif io_ns.anti_version == 4:
        for each_gene in secondary_model_ns.cluster_info_dict.keys():
            if gpr_count == 0:
                gpr_list = each_gene
                gpr_count += 1
            else:
                gpr_list = ' and '.join([gpr_list, each_gene])

    gpr_list = '( %s )' %(gpr_list)

    rxn.gene_reaction_rule = gpr_list

    #Add a new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("%s: %s", rxn, rxn.reaction)

    ##############################
    #Create a transport reaction
    #Create reaction ID
    rxn = Reaction("Transport_" + secondary_model_ns.product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Add a substrate metabolite
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(secondary_model_ns.product+'_c')):-1})

    #Add product metabolite(s)
    product_e = Metabolite(secondary_model_ns.product+"_e", name='', compartment='e')
    rxn.add_metabolites({product_e:1})

    #Add the new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("%s: %s", rxn, rxn.reaction)

    ##############################
    #Create an exchange reaction
    #Create reaction ID
    rxn = Reaction("Ex_"+secondary_model_ns.product)

    #Reversibility / Lower and upper bounds
    rxn.reversibility = 0 # 1: reversible 0: irreversible
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Add a substrate metabolite
    rxn.add_metabolites({target_model.metabolites.get_by_id(str(product_e)):-1})

    #Add a new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("%s: %s", rxn, rxn.reaction)

    return target_model


def check_producibility_sec_met(target_model, io_ns, secondary_model_ns):

    obj_rxn = linear_reaction_coefficients(target_model).keys()[0].id
    target_model.reactions.get_by_id(obj_rxn).objective_coefficient = 0
    target_model.reactions.get_by_id("Ex_"+secondary_model_ns.product).objective_coefficient = 1

    #Model reloading and overwrtting are necessary for model stability
    #Without these, model does not produce an accurate prediction
    target_model = utils.stabilize_model(
            target_model, io_ns.outputfolder5, secondary_model_ns.product)

    flux_dist = target_model.optimize()
    logging.debug("Flux: %s" %flux_dist.objective_value)

    target_model.reactions.get_by_id(obj_rxn).objective_coefficient = 1
    target_model.reactions.get_by_id("Ex_"+secondary_model_ns.product).objective_coefficient = 0

    return target_model, flux_dist


def get_sec_met_monomers(sec_met_list, secondary_model_ns):

    for metab in secondary_model_ns.metab_coeff_dict.keys():
        if secondary_model_ns.metab_coeff_dict[metab] <0:
            if metab not in sec_met_list:
                sec_met_list.append(metab)

    return sec_met_list

