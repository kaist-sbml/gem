
import logging
import os
import pickle
from gmsm.secondary_model.antismash_monomer_info import get_std_id_from_antismash_id
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
                break


#Exract all the information associated with a particular locus_tag for the selected region
def get_region_info_from_seq_record(seq_record, secondary_model_ns):

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


def get_region_product(seq_record, order, region_nr, io_ns, secondary_model_ns):

    for feature in seq_record.features:
        
        if feature.type == 'region':
            if feature.location.start == secondary_model_ns.region_loc1 and \
            feature.location.end == secondary_model_ns.region_loc2:
                product_list = feature.qualifiers.get('product')
                
                #connecting all product of a region
                product = '_'.join(product_list)
                #Handle legacy problem
                product2 = product.replace('-','_')
                
                if len(io_ns.seq_record_BGC_num_lists) == 1:
                    product3 = "Region"+str(region_nr)+"_"+product2.lower()
                else:
                    product3 = "Region"+str(order)+"."+str(region_nr)+"_"+product2.lower()

    secondary_model_ns.product = product3


#Output: e.g., {'SAV_943_M1':['mmal', 'Ethyl_mal', 'pk']}
def get_region_monomers(seq_record, secondary_model_ns):

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
            if 'PKS_AT' in feature.qualifiers.get('aSDomain') or \
            'AMP-binding' in feature.qualifiers.get('aSDomain'):
                
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

        #NRPS analyzed with antiSMASH 8.0
        if len(secondary_model_ns.locustag_monomer_dict[each_module]) == 1:
            #index [0]: substrate consensus
            priority_list = [0]
            biggid_met = get_biggid(priority_list, each_module, secondary_model_ns)

        elif len(secondary_model_ns.locustag_monomer_dict[each_module]) == 3:
            # PKS_AT analyzed with antiSMASH 8.0
            # locustag_monomer_dict[each_module] for pks
            # Index [0]: consensus
            # Index [1]: PKS signature
            # Index [2]: Minowa
            # Priority: consensus > PKS signature > Minowa
            priority_list = [0, 2, 1]
            biggid_met = get_biggid(priority_list, each_module, secondary_model_ns)

        else:
            biggid_met = None
            
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

    #ID and name
    rxn = Reaction(secondary_model_ns.product)
    rxn.name = rxn.id

    #Reversibility / Lower and upper bounds
    rxn.lower_bound = 0
    rxn.upper_bound = 1000

    #Metabolites and their stoichiometric coeff's
    for metab in secondary_model_ns.metab_coeff_dict.keys():

        #Consider only metabolites consumed or produced
        if secondary_model_ns.metab_coeff_dict[metab] != 0:
            metab_compt = '_'.join([metab,'c'])

            #Add metabolites already in the model
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

            #Add metabolites not available in the MNXref sbml
            else:
                logging.debug("Metabolite (MNXM ID) %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:secondary_model_ns.metab_coeff_dict[metab]})

    #GPR association
    gpr_count = 0
    for each_gene in secondary_model_ns.region_info_dict.keys():
        if gpr_count == 0:
            gpr = each_gene
            gpr_count += 1
        else:
            gpr = ' and '.join([gpr, each_gene])

    gpr = '%s' %(gpr)

    rxn.gene_reaction_rule = gpr

    #Add a new reaction to the model
    target_model.add_reaction(rxn)

    logging.debug("%s: %s", rxn, rxn.reaction)

    ##############################
    #Create a transport reaction
    #Create reaction ID and name
    rxn = Reaction("Transport_" + secondary_model_ns.product)
    rxn.name = rxn.id

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
    #Create reaction ID and name
    rxn = Reaction("EX_"+secondary_model_ns.product)
    rxn.name = rxn.id

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

    obj_rxn = list(linear_reaction_coefficients(target_model).keys())[0].id
    target_model.reactions.get_by_id(obj_rxn).objective_coefficient = 0
    target_model.reactions.get_by_id("EX_"+secondary_model_ns.product).objective_coefficient = 1

    #Model reloading and overwrtting are necessary for model stability
    #Without these, model does not produce an accurate prediction
    target_model = utils.stabilize_model(
            target_model, io_ns.outputfolder5, secondary_model_ns.product)

    flux_dist = target_model.optimize()
    logging.debug("Flux: %s" %flux_dist.objective_value)

    target_model.reactions.get_by_id(obj_rxn).objective_coefficient = 1
    target_model.reactions.get_by_id("EX_"+secondary_model_ns.product).objective_coefficient = 0

    return target_model, flux_dist


def get_sec_met_monomers(sec_met_list, secondary_model_ns):

    for metab in secondary_model_ns.metab_coeff_dict.keys():
        if secondary_model_ns.metab_coeff_dict[metab] <0:
            if metab not in sec_met_list:
                sec_met_list.append(metab)

    return sec_met_list

