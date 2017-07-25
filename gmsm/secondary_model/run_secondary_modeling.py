
import cobra
import logging
import pickle
from sec_met_rxn_generation import(
    get_cluster_location,
    get_cluster_info_from_seq_record,
    get_cluster_product,
    get_cluster_monomers,
    get_all_metab_coeff,
    get_pickles,
    add_sec_met_rxn,
    check_producibility_sec_met,
    get_sec_met_monomers
)
from gapfilling import(
    get_unique_nonprod_monomers_list,
    add_transport_exchange_rxn_nonprod_monomer,
    check_producibility_nonprod_monomer,
    check_gapfill_rxn_biomass_effects,
    add_gapfill_rxn_target_model
)


def run_secondary_modeling(seq_record, target_model, options):
    prod_sec_met_dict = {}
    nonprod_sec_met_dict = {}

    cluster_nr = 1

    while cluster_nr <= options.total_cluster:
        logging.info("Generating reactions for Cluster %s.." %cluster_nr)
        target_model = run_sec_met_rxn_generation(
                 seq_record, cluster_nr,
                 target_model,
                 prod_sec_met_dict, nonprod_sec_met_dict,
                 options)

        cluster_nr += 1

    return target_model


def run_sec_met_rxn_generation(seq_record, cluster_nr, target_model, prod_sec_met_dict,
                                nonprod_sec_met_dict, options):

    get_cluster_location(seq_record, cluster_nr, options)

    get_cluster_info_from_seq_record(seq_record, options)

    get_cluster_product(seq_record, cluster_nr, options)

    if 't1pks' in options.product or 'nrps' in options.product:
        get_cluster_monomers(options)

        get_all_metab_coeff(options)

        get_pickles(options)

        target_model = add_sec_met_rxn(target_model, options)

        target_model, flux_dist = check_producibility_sec_met(target_model, options)

        if flux_dist.objective_value < float(options.cobrapy.non_zero_flux_cutoff):
            nonprod_sec_met_list = []
            nonprod_sec_met_metab_list = get_sec_met_monomers(nonprod_sec_met_list, options)
            nonprod_sec_met_dict[options.product] = nonprod_sec_met_metab_list
        else:
            prod_sec_met_list = []
            prod_sec_met_metab_list = get_sec_met_monomers(prod_sec_met_list, options)
            prod_sec_met_dict[options.product] = prod_sec_met_metab_list

    else:
        logging.debug("This BGC does not belong to 't1pks', 'nrps' or their hybird")

    if cluster_nr == options.total_cluster:
        options.prod_sec_met_dict = prod_sec_met_dict
        options.nonprod_sec_met_dict = nonprod_sec_met_dict

    return target_model


def get_target_nonprod_monomers_for_gapfilling(target_model, options):
    logging.info("Producibility of secondary metabolites (gap-filling needed)..")

    unique_nonprod_monomers_list = get_unique_nonprod_monomers_list(options)

    #Some monomers used for nonproducible secondary metabolites can be produced from an initial target_model
    #They need to be excluded from the list for gap-filling targets
    adj_unique_nonprod_monomers_list = []

    for nonprod_monomer in unique_nonprod_monomers_list:
        target_model_monomer = add_transport_exchange_rxn_nonprod_monomer(target_model,
                               nonprod_monomer, options)
        target_model_monomer, flux_dist = \
                check_producibility_nonprod_monomer(target_model_monomer,
                                                    nonprod_monomer)
        if flux_dist.objective_value < float(options.cobrapy.non_zero_flux_cutoff):
            adj_unique_nonprod_monomers_list.append(nonprod_monomer)
        else:
            continue

    logging.debug("Unique set of metabolites for gap-filling: %s",
                    adj_unique_nonprod_monomers_list)

    options.adj_unique_nonprod_monomers_list = adj_unique_nonprod_monomers_list


def run_gapfilling(target_model, options):

    gapfill_rxns2 = []

    for nonprod_monomer in options.adj_unique_nonprod_monomers_list:
        try:
            #Gap-filling via cobrapy
            #TODO: Check downstream functions for 'gapfill_iter' > 1
            gapfill_rxns = cobra.flux_analysis.gapfilling.SMILEY(
                    target_model, '%s_c' %nonprod_monomer,
                    options.mnxref,
                    iterations = int(options.cobrapy.gapfill_iter))
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

        except:
            logging.warning("Gap-filling not possible: target_model with reactions from universal_model does not produce this monomer: %s" %nonprod_monomer)

    #Currently this function causes an error;
    #gap-filling reactions are not added to the model being edited
    #gap_rxns3 = check_gapfill_rxn_biomass_effects(target_model,
    #                           universal_model, gapfill_rxns2, options)
    target_model_complete = add_gapfill_rxn_target_model(target_model,
                            options.mnxref, gapfill_rxns2,options)

    return target_model_complete

