#!/usr/bin/env python

#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import argparse
import glob
import logging
import multiprocessing
import os
import pickle
import sys
import time
import warnings
from cobra.io.sbml import (
    write_cobra_model_to_sbml_file,
    create_cobra_model_from_sbml_file
)
#cobrapy 0.5.8 appears to have an error in this function.
#Argument 'cobra_model' should be manually removed from
#'the_metabolite.remove_from_model(cobra_model)' in
#'cobra.manipulation.delete.prune_unused_metabolites'.
from cobra.manipulation.delete import prune_unused_metabolites
from argparse import Namespace
from gems import check_prereqs
from gems.io.input_file_manager import (
    get_genome_files,
    get_pickles_prunPhase,
    get_pickles_augPhase
    )
from gems.io.output_file_manager import generate_outputs
from gems.homology.bidirect_blastp_analysis import get_homologs
from gems.primary_model.run_primary_modeling import run_prunPhase, run_augPhase
from gems.secondary_model.run_secondary_modeling import (
    run_sec_met_rxn_generation,
    get_target_nonprod_monomers_for_gapfilling,
    get_universal_model,
    run_gapfilling
    )


def main():
    start = time.time()

    parser = argparse.ArgumentParser()

    #General options
    parser.add_argument('-c', '--cpu',
                        dest='cpus',
                        default=multiprocessing.cpu_count(),
                        type=int,
                        help="How many CPUs to use in parallel. (default: %(default)s)")

    group = parser.add_argument_group('Input and output options')
    group.add_argument('-i', '--input',
                        dest='input',
                        default='input',
                        help="Specify input file")
    group.add_argument('-o', '--outputfolder',
                        dest='outputfolder',
                        default='output',
                        help="Specify output directory")

    group = parser.add_argument_group('GEMS options',
                        "At least one of the three options should be selected:"
                        " '-e', '-p-' and '-s'")
    group.add_argument('-m', '--model',
                        dest='orgName',
                        default='sco',
                        choices=['eco','sco'],
                        help="Specify a template model for the target modeling")
    group.add_argument('-e', '--ec',
                        dest='eficaz',
                        action='store_true',
                        default=False,
                        help="Run EC number prediction using EFICAz")
    group.add_argument('-p', '--primary-modeling',
                        dest='pmr_generation',
                        default=False,
                        action='store_true',
                        help='Run primary metabolic modeling')
    group.add_argument('-s', '--secondary-modeling',
                        dest='smr_generation',
                        default=False,
                        action=('store_true'),
                        help="Run secondary metabolic modeling")

    group = parser.add_argument_group('Debugging and logging options')
    group.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        default=False,
                        help="Print verbose status information to stderr")
    group.add_argument('-d', '--debug',
                        dest='debug',
                        action='store_true',
                        default=False,
                        help="Print debugging information to stderr")
    group.add_argument('-w', '--warning',
                        dest='warning',
                        action='store_true',
                        default=False,
                        help="Print UserWarning messages from cobrapy")

    options = parser.parse_args()

    if options.verbose:
        log_level = logging.INFO
    elif options.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING

    logging.basicConfig(format='%(levelname)s: %(message)s', level=log_level)

    #Warning messages from cobrapy turned off by default
    if not options.warning:
        warnings.filterwarnings("ignore")

    #Get genome files only if one of functional options is selected
    if options.eficaz or options.pmr_generation or options.smr_generation:
        #Check prerequisites of executables and libraries
        check_prereqs()

        #Create output folders
        folders = ['1_EFICAz_results', '2_blastp_results',
                    '3_primary_metabolic_model', '4_complete_model', 'tmp_files']

        if '/' in options.outputfolder:
            options.outputfolder = options.outputfolder[:-1]

        #'1_EFICAz_results'
        options.outputfolder1 = options.outputfolder + os.sep + folders[0]
        #'2_blastp_results'
        options.outputfolder2 = options.outputfolder + os.sep + folders[1]
        #'3_primary_metabolic_model'
        options.outputfolder3 = options.outputfolder + os.sep + folders[2]
        #'4_complete_model'
        options.outputfolder4 = options.outputfolder + os.sep + folders[3]
        #'tmp_files'
        options.outputfolder5 = options.outputfolder + os.sep + folders[4]

        for folder in folders:
            if not os.path.isdir(options.outputfolder + os.sep + folder):
                os.makedirs(options.outputfolder + os.sep + folder)

        get_genome_files(options)

    #Primary metabolic modeling
    if options.pmr_generation:
        get_homologs(options)

        model = get_pickles_prunPhase(options)

        modelPrunedGPR = run_prunPhase(model, options)

        if options.targetGenome_locusTag_ec_dict:
            get_pickles_augPhase(options)
            target_model = run_augPhase(modelPrunedGPR, options)
        else:
            logging.warning("No EC_number found in the submitted gbk file")

        #Cleanup of the model
        prune_unused_metabolites(target_model)

        runtime1 = time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start))

        generate_outputs(options.outputfolder3, target_model, runtime1, options)

    #Secondary metabolic modeling
    if options.smr_generation:
        model_file = []
        files = glob.glob(options.outputfolder3 + os.sep + '*.xml')
        model_file = [each_file for each_file in files if '.xml' in each_file]

        if len(model_file) > 0 and '.xml' in model_file[0]:
            model_file = os.path.basename(model_file[0])
            target_model = create_cobra_model_from_sbml_file(
                    options.outputfolder3 + os.sep + model_file)

            logging.info("Generating secondary metabolite biosynthesizing reactions..")
            logging.debug("Total number of clusters: %s" %options.total_cluster)

            if options.total_cluster > 0:
                prod_sec_met_dict = {}
                nonprod_sec_met_dict = {}

                cluster_nr = 1
                while cluster_nr <= options.total_cluster:
                    logging.info("Generating reactions for Cluster %s.." %cluster_nr)
                    target_model = run_sec_met_rxn_generation(cluster_nr,
                        target_model, prod_sec_met_dict, nonprod_sec_met_dict, options)
                    cluster_nr += 1

                get_target_nonprod_monomers_for_gapfilling(target_model, options)

                universal_model = get_universal_model(target_model, options)

                target_model_complete = run_gapfilling(target_model, universal_model, options)

                #Cleanup of the model
                prune_unused_metabolites(target_model_complete)

                runtime2 = time.strftime("Elapsed time %H:%M:%S",
                        time.gmtime(time.time() - start))

                generate_outputs(options.outputfolder4,
                    target_model_complete, runtime2, options)

            else:
                logging.debug("No cluster information available for secondary metabolic modeling")

        else:
            logging.warning("COBRA-compliant SBML file needed")

    if not options.eficaz and not options.pmr_generation and not options.smr_generation:
        logging.warning("No functional options selected")

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


if __name__ == '__main__':
    main()
