
import logging
import os
import pickle
from antismash_monomer_info import get_std_id_from_antismash_id
from cobra import Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from gmsm import utils

def get_region_location(seq_record, options):
    options.region_startlist = []
    options.region_endlist = []
    for feature in seq_record.features:
        if feature.type == 'region':
            options.region_startlist.append(feature.location.start)
            options.region_endlist.append(feature.location.end)


#Exract all the information associated with a particular locus_tag for the selected cluster
def get_region_info_from_seq_record(seq_record, region_nr, options):
    region_info_dict = {}
    for feature in seq_record.features:
        if feature.type == 'CDS':
            if feature.location.start >= options.region_startlist[region_nr-1] \
                    and feature.location.end <= options.region_endlist[region_nr-1]:
                qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                if feature.qualifiers.get('sec_met_domain'):
                    region_info_dict[qualifier_locus_tag] = \
                            feature.qualifiers.get('sec_met_domain')

    options.region_info_dict = region_info_dict
             

def get_region_product(seq_record, region_nr, options):
    for feature in seq_record.features:
        if feature.type == 'region':
            if feature.location.start == options.region_startlist[region_nr-1] and \
            feature.location.end == options.region_endlist[region_nr-1]:
                product_list = feature.qualifiers.get('product')
                for i in range(len(product_list)):
                    # connecting all product of a region
                    if i == 0:
                        product = product_list[i].replace('-','_')
                        if float(region_nr) < 10:
                            product2 = "Region0"+str(region_nr)+"_"+product
                        else:
                            product2 = "Region"+str(region_nr)+"_"+product
                    else:
                        product = product_list[i].replace('-','_')
                        product2 = product2+'_'+product
                    options.product = product2      
        
        
        
#Exract all the information associated with a particular locus_tag for the selected region
#Output: e.g., {'SAV_943_M1':['mmal', 'Ethyl_mal', 'pk']}
def get_region_monomers(seq_record, region_nr, options):

    options.cds_info_dict = {}
    options.locustag_monomer_dict = {}

    for feature in seq_record.features:  
        # find location of CDS involved in SM
        if feature.type == 'CDS':       
            module_count = 0
            options.cds_loc1 = feature.location.start
            options.cds_loc2 = feature.location.end
        # find antiSMASH domain in range of CDS    
        if feature.type == 'aSDomain':       
            if feature.location.start >= options.cds_loc1 \
            and feature.location.end <= options.cds_loc2:
                # collect monomer information in one region
                if feature.location.start >= options.region_startlist[region_nr-1] and \
                feature.location.end <= options.region_endlist[region_nr-1]:
                    qualifier_locus_tag = feature.qualifiers.get('locus_tag')[0]
                    module_number = qualifier_locus_tag + '_M' + str(module_count)
                    if feature.qualifiers.get('specificity'):
                        sec_met_info = feature.qualifiers.get('specificity')
                        if 'KR ' not in sec_met_info[0]:     # exclude KR stochiometry information
                            options.cds_info_dict[module_number] = \
                                feature.qualifiers.get('specificity')
                            module_count += 1

    for each_gene in options.cds_info_dict.keys():
        monomer_list = []
        for sec_met_info in options.cds_info_dict[each_gene]:
            sec_met = sec_met_info.split(': ')
            if 'KR ' not in sec_met[0]:
                monomer = sec_met[1]
                monomer_list.append(monomer)
            options.locustag_monomer_dict[each_gene] = monomer_list
    
        

def get_biggid(priority_list, each_module, options):
    for i in priority_list:
        aSid_met = options.locustag_monomer_dict[each_module][i]
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
def get_all_metab_coeff(options):
        metab_coeff_dict = {}
        for each_module in options.locustag_monomer_dict:
            # PKS analyzed with antiSMASH 5.0
            # locustag_monomer_dict[each_module] for pks
            # Index [0]: consensus
            # Index [1]: Minowa
            # Index [2]: PKS signature
            # Priority: consensus > PKS signature > Minowa
            if len(options.locustag_monomer_dict[each_module]) == 3:
                priority_list = [0, 2, 1]
                biggid_met = get_biggid(priority_list, each_module, options)
            elif len(options.locustag_monomer_dict[each_module]) == 1:
                priority_list = [0]
                biggid_met = get_biggid(priority_list, each_module, options)
            # Filter modules without 'Substrate specificity predictions'
            # e.g., 'B446_13275_M0'
            # {'B446_13415_M0': ['gly,ala,val,leu,ile,abu,iva', 'ile', 'val', 'nrp'], 'B446_13275_M0': [], 'B446_13350_M0': ['leu', 'ala', 'sar', 'nrp'], 'B446_13445_M0': ['pro,pip', 'N/A', 'orn', 'nrp']}
            if options.locustag_monomer_dict[each_module] and biggid_met:
                if biggid_met not in metab_coeff_dict:
                    metab_coeff_dict[biggid_met] = -1
                else:
                    metab_coeff_dict[biggid_met] -= 1
        # Add secondary metabolite product to the reaction
        metab_coeff_dict[options.product] = 1
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

            elif 'Region' in metab:
                logging.debug("Secondary metabolite ('Region') %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

            #Add metabolites not available in the MNXref sbml
            else:
                logging.debug("Metabolite (MNXM ID) %s: To be added" %metab)
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:options.metab_coeff_dict[metab]})

    #GPR association
    gpr_count = 0
    for each_gene in options.region_info_dict.keys():
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

    logging.debug("%s: %s", rxn, rxn.reaction)

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

    logging.debug("%s: %s", rxn, rxn.reaction)

    return target_model


def check_producibility_sec_met(target_model, options):

    obj_rxn = linear_reaction_coefficients(target_model).keys()[0].id
    target_model.reactions.get_by_id(obj_rxn).objective_coefficient = 0
    target_model.reactions.get_by_id("Ex_"+options.product).objective_coefficient = 1

    #Model reloading and overwrtting are necessary for model stability
    #Without these, model does not produce an accurate prediction
    target_model = utils.stabilize_model(
            target_model, options.outputfolder5, options.product)

    flux_dist = target_model.optimize()
    logging.debug("Flux: %s" %flux_dist.objective_value)

    target_model.reactions.get_by_id(obj_rxn).objective_coefficient = 1
    target_model.reactions.get_by_id("Ex_"+options.product).objective_coefficient = 0

    return target_model, flux_dist


def get_sec_met_monomers(sec_met_list, options):

    for metab in options.metab_coeff_dict.keys():
        if options.metab_coeff_dict[metab] <0:
            if metab not in sec_met_list:
                sec_met_list.append(metab)

    return sec_met_list

