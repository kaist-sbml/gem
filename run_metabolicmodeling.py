#!/usr/bin/env python

#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import argparse
import logging
import multiprocessing
import os
import pickle
import sys
import time

from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
from argparse import Namespace
from modeling.io.input_file_manager import (
    get_genome_files,
    get_pickles_add_rxn,
    get_pickles_prunPhase
)
from modeling.io.input_file_manager import (
    generate_outputs_primary_model,
    generate_outputs_secondary_model
)
from modeling.homology import get_homologs
from modeling.primary_model import run_prunPhase
from modeling.primary_model import run_augPhase
from modeling.run_sec_met_rxn_generation import (
    run_sec_met_rxn_generation,
    prep_network_for_gapfilling,
    get_target_nonprod_monomers_for_gapfilling,
    run_gapfilling
)

def main():
    start = time.time()

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--model',
                        dest='orgName',
                        default='sco',
                        choices=['eco','sco'],
                        help="Specify a template model for the target modeling")
    parser.add_argument('--disable-modeling',
                        dest='pmr_generation',
                        default=True,
                        action='store_false',
                        help='Disable primary metabolic modeling')
    parser.add_argument('-s', '--smr',
                        dest='smr_generation',
                        default=False,
                        action=('store_true'),
                        help="Specify whether to run secondary metabolic modeling")
    parser.add_argument('-i', '--input',
                        dest='input',
                        default='input',
                        help="Specify input directory")
    parser.add_argument('-o', '--output',
                        dest='output',
                        default='output',
                        help="Specify output directory")
    parser.add_argument('-e', '--ec',
                        dest='eficaz',
                        action='store_true',
                        default=False,
                        help="Run EC number prediction using EFICAz")
    parser.add_argument('-c', '--cpu',
                        dest='cpus',
                        default=multiprocessing.cpu_count(),
                        type=int,
                        help="How many CPUs to use in parallel. (default: %(default)s)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        default=False,
                        help="Print verbose status information to stderr")
    parser.add_argument('-d', '--debug',
                        dest='debug',
                        action='store_true',
                        default=False,
                        help="Print debugging information to stderr")

    options = parser.parse_args()

    if options.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    #Create output folders
    folders = ['0_EFICAz_results', '1_blastp_results', '2_primary_metabolic_model', '3_temp_models', '4_complete_model']

    for folder in folders:
        if not os.path.isdir(options.output+'/'+folder):
            os.makedirs(options.output+'/'+folder)

    options.outputfoldername = options.output+'/'+folders[0]

    get_genome_files(options)

    get_pickles_add_rxn(options)

    if options.pmr_generation:
        get_homologs(options)

        model = get_pickles_prunPhase(options)

        modelPrunedGPR = run_prunPhase(model, options)

        if options.targetGenome_locusTag_ec_dict:

            target_model = run_augPhase(modelPrunedGPR, options)
        else:
            logging.warning("No EC_number found in the submitted gbk file")

        generate_outputs_primary_model(model, modelPrunedGPR, target_model, options)

    #Secondary metabolic modeling
    if options.smr_generation:
        for model_file in os.listdir(options.output+'/'+'2_primary_metabolic_model'):
            if model_file.endswith('.xml'):
                target_model = create_cobra_model_from_sbml_file(options.output+'/'+'2_primary_metabolic_model/'+model_file)

                logging.info("Generating secondary metabolite biosynthesizing reactions..")
                cluster_nr = 1
                logging.debug("Total number of clusters: %s" %options.total_cluster)

                prod_sec_met_dict = {}
                nonprod_sec_met_dict = {}

                while cluster_nr <= options.total_cluster:
                    logging.debug("Cluster number: %s" %cluster_nr)
                    target_model = run_sec_met_rxn_generation(cluster_nr, target_model, prod_sec_met_dict, nonprod_sec_met_dict, options)
                    cluster_nr += 1

                target_model2, universal_model = prep_network_for_gapfilling(target_model, options)

                adj_unique_nonprod_monomers_list = get_target_nonprod_monomers_for_gapfilling(target_model, options)

                target_model_complete = run_gapfilling(target_model, target_model2, adj_unique_nonprod_monomers_list, universal_model, options)

                generate_outputs_secondary_model(target_model_complete, options)

            elif not model_file:
                logging.warning("COBRA-compliant SBML file needed")

    if not options.pmr_generation and not options.smr_generation:
        logging.warning("Either primary or secondary metabolic modeling should be performed")

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


if __name__ == '__main__':
    main()
