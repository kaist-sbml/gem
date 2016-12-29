
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import cobra
import logging
import pickle
from sec_met_rxn_generation import(
    get_cluster_location,
    get_cluster_info_from_seq_record,
    get_cluster_product,
    get_cluster_domain,
    get_cluster_monomers,
    get_cluster_module,
    get_currency_metabolites,
    get_total_currency_metab_coeff,
    get_all_metab_coeff,
    get_pickles,
    add_sec_met_rxn,
    check_producibility_sec_met,
    get_monomers_nonprod_sec_met,
    get_monomers_prod_sec_met
)
from gapfilling import(
    get_unique_nonprod_monomers_list,
    add_transport_exchange_rxn_nonprod_monomer,
    check_producibility_nonprod_monomer,
    check_gapfill_rxn_biomass_effects,
    add_gapfill_rxn_target_model
)


def run_sec_met_rxn_generation(cluster_nr, target_model, prod_sec_met_dict,
                                nonprod_sec_met_dict, options):

    get_cluster_location(cluster_nr, options)

    get_cluster_info_from_seq_record(options)

    get_cluster_product(cluster_nr, options)

    if 't1pks' in options.product or 'nrps' in options.product:
        get_cluster_domain(options)

        get_cluster_monomers(options)

        get_cluster_module(options)

        get_currency_metabolites(options)

        get_total_currency_metab_coeff(options)

        get_all_metab_coeff(options)

        get_pickles(options)

        target_model = add_sec_met_rxn(target_model, options)

        target_model = check_producibility_sec_met(target_model, options)

        if target_model.solution.f < 0.0001:
            nonprod_sec_met_metab_list = get_monomers_nonprod_sec_met(options)
            nonprod_sec_met_dict[options.product] = nonprod_sec_met_metab_list
        else:
            prod_sec_met_metab_list = get_monomers_prod_sec_met(options)
            prod_sec_met_dict[options.product] = prod_sec_met_metab_list

    else:
        logging.debug("Not type I polyketide synthase ('t1pks'), nonribosomal synthetase ('nrps') or their hybird")

    if cluster_nr == options.total_cluster:
        options.prod_sec_met_dict = prod_sec_met_dict
        options.nonprod_sec_met_dict = nonprod_sec_met_dict

    return target_model


def get_target_nonprod_monomers_for_gapfilling(target_model, options):
    logging.info("Gap-filling for the production of secondary metabolites..")

    unique_nonprod_monomers_list = get_unique_nonprod_monomers_list(options)

    #Some monomers used for nonproducible secondary metabolites can be produced from an initial target_model
    #They need to be excluded from the list for gap-filling targets
    adj_unique_nonprod_monomers_list = []

    for nonprod_monomer in unique_nonprod_monomers_list:
        target_model_monomer = add_transport_exchange_rxn_nonprod_monomer(target_model,
                               nonprod_monomer, options)
        target_model_monomer = check_producibility_nonprod_monomer(target_model_monomer,
                               nonprod_monomer)
        if target_model_monomer.solution.f < 0.0001:
            adj_unique_nonprod_monomers_list.append(nonprod_monomer)
        else:
            continue

    logging.debug("Adjusted unique_nonprod_monomers_list: %s" %adj_unique_nonprod_monomers_list)

    options.adj_unique_nonprod_monomers_list = adj_unique_nonprod_monomers_list


def get_universal_model(target_model, options):


    universal_model = pickle.load(open("./gems/io/data/input2/universal_model.p","rb"))

    return universal_model


def run_gapfilling(target_model, universal_model, options):

    gapfill_rxns2 = []

    for nonprod_monomer in options.adj_unique_nonprod_monomers_list:
        #Gap-filling via cobrapy
        gapfill_rxns = cobra.flux_analysis.gapfilling.SMILEY(
                target_model, '%s_c' %nonprod_monomer, universal_model)
        logging.debug('gapfill_rxns: %s' %gapfill_rxns)

        #'gapfill_rxns' is a list of list:
        #e.g., [[<Reaction HKSR9 at 0x7f706aa8cad0>,
        #<Reaction MNXR3151_reverse at 0x7f706a86d7d0>]]
        #Element in a list of list is Class.
        if len(gapfill_rxns[0]) > 0:
            for i in gapfill_rxns[0]:
                if '_reverse' in str(i):
                    #Check next gap-filling rxn already covered
                    if str(i)[:-8] not in gapfill_rxns2:
                        gapfill_rxns2.append(str(i)[:-8])
                elif '_reverse' not in str(i):
                    #Check next gap-filling rxn already covered
                    if str(i) not in gapfill_rxns2:
                        gapfill_rxns2.append(str(i))

        else:
            logging.warning("Gap-filling not possible: target_model with reactions from universal_model does not produce this monomer: %s" %nonprod_monomer)

    #Currently this function causes an error;
    #gap-filling reactions are not added to the model being edited
    #gap_rxns3 = check_gapfill_rxn_biomass_effects(target_model,
    #                           universal_model, gapfill_rxns2, options)
    target_model_complete = add_gapfill_rxn_target_model(target_model,
                            universal_model, gapfill_rxns2,options)

    return target_model_complete

