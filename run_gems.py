#!/usr/bin/env python

#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import argparse
import cobra
import copy
import glob
import logging
import multiprocessing
import os
import pickle
import sys
import time
import warnings

# cobrapy == 0.5.11 should be used, which now has a fixed  function:
#'cobra.manipulation.delete.prune_unused_metabolites'.
from cobra.manipulation.delete import prune_unused_metabolites
from gems import check_prereqs, utils
from gems.config import load_config
from gems.eficaz import getECs
from gems.io.input_file_manager import (
    setup_outputfolders,
    check_input_filetype,
    get_target_gbk,
    get_fasta_files,
    get_pickles_prunPhase,
    get_pickles_augPhase
    )
from gems.io.output_file_manager import generate_outputs
from gems.homology.bidirect_blastp_analysis import get_homologs
from gems.primary_model.run_primary_modeling import run_prunPhase, run_augPhase
from gems.secondary_model.run_secondary_modeling import (
    run_sec_met_rxn_generation,
    get_target_nonprod_monomers_for_gapfilling,
    run_gapfilling
    )


def main():
    start = time.time()
    usage = \
            '\rGEMS version {version} ({git_log})\n\nusage: {usage}\n----------------------------------------------------------------------------------'\
            .format(version=utils.get_version(),
                    git_log=utils.get_git_log(),
                    usage='run_gems.py [-h] [Resource management] [Input and output setting] [GEMS modeling options] [Debugging and logging options]')
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            usage=usage
            )

    #General options
    group = parser.add_argument_group('Resource management')
    group.add_argument('-c', '--cpu',
                        dest='cpus',
                        default=multiprocessing.cpu_count(),
                        type=int,
                        help="How many CPUs to use in parallel. (default: %(default)s)")

    group = parser.add_argument_group('Input and output setting')
    group.add_argument('-i', '--input',
                        dest='input',
                        default='input',
                        help="Specify input file")
    group.add_argument('-o', '--outputfolder',
                        dest='outputfolder',
                        default='output',
                        help="Specify output directory (optional)")

    group = parser.add_argument_group('Template model options',
                        "Select a biologically close organism")
    group.add_argument('-m', '--model',
                        dest='orgName',
                        default='sco',
                        choices=['bsu', 'eco','mtu','ppu','sco'],
                        help="Specify a template model for the target modeling\n"
                        "'bsu': iYO844; Bacillus subtilis subsp. subtilis str. 168\n"
                        "'eco': iAF1260; Escherichia coli str. K-12 substr. MG1655\n"
                        "'mtu': iNJ661; Mycobacterium tuberculosis H37Rv\n"
                        "'ppu': iJN746; Pseudomonas putida KT2440\n"
                        "'sco': iMK1208; Streptomyces coelicolor A3(2)"
                        )

    group = parser.add_argument_group('GEMS modeling options',
                        "At least one of the three options should be selected:"
                        " '-e', '-p-' and '-s'")
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
    group.add_argument('-V', '--version',
                        dest='version',
                        action='store_true',
                        default=False,
                        help="Show the program version and exit")

    options = parser.parse_args()

    utils.setup_logging(options)

    #Create output folders
    setup_outputfolders(options)

    utils.setup_logfile_format(options)

    if options.version:
        print 'GEMS version %s (%s)' %(utils.get_version(), utils.get_git_log())
        sys.exit(0)

    #Warning messages from cobrapy turned off by default
    if not options.warning:
        warnings.filterwarnings("ignore")

    logging.info('Starting GEMS ver. %s (%s)', utils.get_version(), utils.get_git_log())

    check_input_filetype(options)

    #Get genome files only if one of functional options is selected
    if options.eficaz or options.pmr_generation or options.smr_generation:

        #Load config data
        load_config(options)

        #Check prerequisites of executables and libraries
        check_prereqs(options)

        #Create output folders
#        setup_outputfolders(options)

#        utils.setup_logfile_format(options)

    # EC number prediction
    if options.eficaz:
        seq_record = get_target_gbk(options)

        if options.eficaz_path and options.targetGenome_locusTag_aaSeq_dict:
            getECs(seq_record, options)
        elif not options.eficaz_path:
            logging.warning("EFICAz not found")
        elif not options.targetGenome_locusTag_aaSeq_dict:
            logging.warning("EFICAz not implemented;")
            logging.warning("No amino acid sequences found in the submitted gbk file")

    # Primary metabolic modeling
    if options.pmr_generation:
        seq_record = get_target_gbk(options)
        get_fasta_files(options)

        if options.targetGenome_locusTag_aaSeq_dict:
            get_homologs(options)
            model = get_pickles_prunPhase(options)
            modelPrunedGPR = run_prunPhase(model, options)

            if options.targetGenome_locusTag_ec_dict:
                get_pickles_augPhase(options)
                target_model = run_augPhase(modelPrunedGPR, options)
            else:
                logging.warning("Primary metabolic modeling not implemented;")
                logging.warning("No EC_numbers found in the submitted gbk file")

            try:
                prune_unused_metabolites(target_model)
            except:
                prune_unused_metabolites(modelPrunedGPR)

            runtime1 = time.strftime("Elapsed time %H:%M:%S",
                    time.gmtime(time.time() - start))

            try:
                generate_outputs(options.outputfolder3, runtime1, options,
                        cobra_model = target_model)
            except:
                generate_outputs(options.outputfolder3, runtime1, options,
                        cobra_model = modelPrunedGPR)
        else:
            logging.warning("Primary metabolic modeling not implemented;")
            logging.warning("No amino acid sequences found in the submitted gbk file")

    # Secondary metabolic modeling
    if options.smr_generation:
        if 'targetGenome_locusTag_aaSeq_dict' not in options:
            seq_record = get_target_gbk(options)

        model_file = []
        files = glob.glob(options.outputfolder3 + os.sep + '*.xml')
        model_file = [each_file for each_file in files if '.xml' in each_file]

        if len(model_file) > 0 and '.xml' in model_file[0]:
            model_file = os.path.basename(model_file[0])
            target_model = cobra.io.read_sbml_model(
                    os.path.join(options.outputfolder3, model_file))

            logging.info("Generating secondary metabolite biosynthesizing reactions..")
            logging.debug("Total number of clusters: %s" %options.total_cluster)

            if options.total_cluster > 0:
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

                target_model_no_gapsFilled = copy.deepcopy(target_model)

                get_target_nonprod_monomers_for_gapfilling(target_model, options)

                target_model_complete = run_gapfilling(target_model, options)

                prune_unused_metabolites(target_model_complete)

                runtime2 = time.strftime("Elapsed time %H:%M:%S",
                        time.gmtime(time.time() - start))

                generate_outputs(options.outputfolder4,
                        runtime2, options,
                        cobra_model_no_gapFilled = target_model_no_gapsFilled,
                        cobra_model = target_model_complete)
            else:
                logging.warning("Secondary metabolic modeling not implemented;")
                logging.warning("No cluster information found in the submitted gbk file")
        else:
            logging.warning("Secondary metabolic modeling not implemented;")
            logging.warning("COBRA-compliant SBML file needed")

    if not options.eficaz and not options.pmr_generation and not options.smr_generation:
        logging.warning("No functional options enabled")

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


if __name__ == '__main__':
    main()
