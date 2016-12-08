
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import pickle
rom sec_met_rxn_generation import(
    get_cluster_location,
    get_cluster_info_from_seq_record,
    get_cluster_product,
    get_cluster_domain,
    get_cluster_monomers,
    get_cluster_module,
    get_currency_metabolites,
    get_total_currency_metab_coeff,
    get_all_metab_coeff,
    add_sec_met_rxn,
    check_producibility_sec_met,
    get_monomers_nonprod_sec_met,
    get_monomers_prod_sec_met
)
from gapfilling.gapfill_network_manipulation import(
    get_mnxr_bigg_in_target_model,
    get_mnxr_unique_to_universal_model,
    integrate_target_universal_models,
    get_unique_nonprod_monomers_list,
    add_transport_exchange_rxn_nonprod_monomer,
    check_producibility_nonprod_monomer,
    execute_gapfill,
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

        target_model = add_sec_met_rxn(target_model, options)

        target_model = check_producibility_sec_met(target_model, options)

        if target_model.solution.f < 0.0001:
            nonprod_sec_met_metab_list = get_monomers_nonprod_sec_met(options)
            nonprod_sec_met_dict[options.product] = nonprod_sec_met_metab_list
        else:
            prod_sec_met_metab_list = get_monomers_prod_sec_met(options)
            prod_sec_met_dict[options.product] = prod_sec_met_metab_list

    else:
        logging.warning("Not type I polyketide synthase ('t1pks'), nonribosomal synthetase ('nrps') or their hybird")

    if cluster_nr == options.total_cluster:
        options.prod_sec_met_dict = prod_sec_met_dict
        options.nonprod_sec_met_dict = nonprod_sec_met_dict

    return target_model


def prep_network_for_gapfilling(target_model, options):

    logging.info("Gap-filling for the production of secondary metabolites..")
    logging.debug("Step 1: Network manipulation for gap-filling process..")

    universal_model = pickle.load(open("./modeling/io/data/input2/universal_model.p","rb"))

    logging.debug("Retrieving reaction information from target_model an universal_model..")
    get_mnxr_bigg_in_target_model(target_model, options)

    get_mnxr_unique_to_universal_model(universal_model, options)

    logging.debug("Merging target_model and universal_model..")
    target_model2 = integrate_target_universal_models(target_model,
                    universal_model, options)

    return target_model2, universal_model


def get_target_nonprod_monomers_for_gapfilling(target_model, options):
    logging.debug("Step 2: Optimization-based gap-filling process..")

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

    return adj_unique_nonprod_monomers_list


def run_gapfilling(target_model, target_model2, adj_unique_nonprod_monomers_list, universal_model, options):

    for nonprod_monomer in adj_unique_nonprod_monomers_list:
        target_model_temp = add_transport_exchange_rxn_nonprod_monomer(target_model2,
                            nonprod_monomer, options)
        target_model_temp = check_producibility_nonprod_monomer(target_model_temp,
                            nonprod_monomer)
        target_model_temp.optimize()

        #Run gap-filling procedure only for monomers producible from target_model with reactions from universal_model
        if target_model_temp.solution.f > 0:
            added_reaction = execute_gapfill(target_model_temp, nonprod_monomer, options)
            added_reaction2 = check_gapfill_rxn_biomass_effects(target_model,
                              universal_model, added_reaction, options)
            target_model_complete = add_gapfill_rxn_target_model(target_model,
                                    universal_model, added_reaction2, options)
        else:
            logging.warning("Gap-filling not possible: target_model with reactions from universal_model does not produce this monomer: %s" %nonprod_monomer)

    return target_model_complete

