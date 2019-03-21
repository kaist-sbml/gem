#!/usr/bin/env python

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
from gmsm import check_prereqs, utils
from gmsm.config import load_config
from gmsm.eficaz import getECs1, getECs2
from gmsm.io.input_file_manager import (
    make_folder,
    setup_outputfolders,
    show_input_options,
    check_input_filetype,
    get_target_genome_from_input,
    get_eficaz_file,
    get_fasta_files,
    get_pickles_prunPhase,
    get_pickles_augPhase,
    get_locustag_comp_dict
    )
from gmsm.io.output_file_manager import generate_outputs, remove_tmp_model_files
from gmsm.homology.bidirect_blastp_analysis import get_homologs
from gmsm.primary_model.run_primary_modeling import run_prunPhase, run_augPhase
from gmsm.secondary_model.run_secondary_modeling import (
    run_secondary_modeling,
    get_target_nonprod_monomers_for_gapfilling,
    run_gapfilling
    )


def main():
    start = time.time()
    usage = \
            '\rGEMS version {version} ({git_log})\n\nusage: {usage}\n----------------------------------------------------------------------------------'\
            .format(version=utils.get_version(),
                    git_log=utils.get_git_log(),
                    usage='run_gmsm.py [-h] [Resource management] [Input and output setting] [GEMS modeling options] [Debugging and logging options]')
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
                        default=False,
                        help="Specify input file")
    group.add_argument('-o', '--outputfolder',
                        dest='outputfolder',
                        default='output',
                        help="Specify output directory (optional)")

    group = parser.add_argument_group('Template model options',
                        "Select a biologically close organism\n"
                        "(with corresponding PMID in parenthesis next to the model name)\n"
                        "A default template GEM is 'sco' (24623710)")
    group.add_argument('-m', '--model',
                        dest='orgName',
                        default='sco',
                        choices=['bsu', 'cre', 'eco','mtu','nsal','ppu','sco'],
                        help=
                        "'bsu': iYO844 (17573341); Bacillus subtilis subsp. subtilis str. 168\n"
                        #NOTE: metabolite compartments (other than 'c') NOT standardized
                        "'cre': iCre1355 (26485611); Chlamydomonas reinhardtii\n"
                        "'eco': iAF1260 (17593909); Escherichia coli str. K-12 substr. MG1655\n"
                        "'mtu': iNJ661 (17555602); Mycobacterium tuberculosis H37Rv\n"
                        "'nsal': iNS934 (28676050); Nannochloropsis salina\n"
                        "'ppu': iJN746 (18793442); Pseudomonas putida KT2440\n"
                        "'sco': iMK1208 (24623710); Streptomyces coelicolor A3(2)"
                        )

    group = parser.add_argument_group('GEMS modeling options',
                        "At least one of the three options should be selected:"
                        " '-e', '-p' and '-s'\n"
                        "Primary metabolic modeling option ('-p') should be selected "
                        "when using '-E' and/or '-C' options\n"
                        " - Examples:\n"
                        "   '-e -E': NOT acceptable\n"
                        "   '-p -E': Acceptable\n"
                        "   '-s -p -E': Acceptable\n"
                        "   '-s -E': NOT acceptable"
                        )
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
    group.add_argument('-E', '--EFICAz',
                        dest='eficaz_file',
                        default=False,
                        help="Specify EFICAz output file")
    group.add_argument('-C', '--comp',
                        dest='comp',
                        default=False,
                        help="Specify file on subcellular localizations (compartments)")

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

    # Create an output directory for a log file
    make_folder(options.outputfolder)

    utils.setup_logging(options)

    # Create output folders
    setup_outputfolders(options)

    if options.version:
        print 'GEMS version %s (%s)' %(utils.get_version(), utils.get_git_log())
        sys.exit(0)

    utils.check_input_options(options)

    # Warning messages from cobrapy turned off by default
    if not options.warning:
        warnings.filterwarnings("ignore")

    logging.info('Starting GEMS ver. %s (%s)', utils.get_version(), utils.get_git_log())

    show_input_options(options)

    logging.info("Reading input genome files..")
    filetype = check_input_filetype(options)

    # Load config data
    load_config(options)

    # Check prerequisites of executables and libraries
    check_prereqs(options)

    # EC number prediction
    if options.eficaz:
        seq_records = get_target_genome_from_input(filetype, options)

        if options.eficaz_path and \
                options.targetGenome_locusTag_aaSeq_dict and \
                not options.eficaz_file:

            if filetype == 'fasta' or len(seq_records) > 1:
                logging.info("Input file in FASTA format or with multiple records:")
                logging.info("Raw EFICAz output (.txt)  will be generated, not GenBank")

            if len(seq_records) == 1:
                getECs1(options, seq_record = seq_records[0])
            elif len(seq_records) > 1:
                getECs2(options)
        else:
            logging.warning("EFICAz not implemented;")

            if not options.eficaz_path:
                logging.warning("EFICAz not found")
            elif not options.targetGenome_locusTag_aaSeq_dict:
                logging.warning(
                        "No amino acid sequences found in input genome data")

    # Primary metabolic modeling
    if options.pmr_generation:
        if not options.eficaz:
            seq_records = get_target_genome_from_input(filetype, options)

        if options.eficaz_file:
            get_eficaz_file(options)

        get_fasta_files(options)

        if options.targetGenome_locusTag_aaSeq_dict:
            get_homologs(options)
            model = get_pickles_prunPhase(options)
            modelPrunedGPR = run_prunPhase(model, options)

            if options.targetGenome_locusTag_ec_dict:
                get_pickles_augPhase(options)

                if options.comp:
                    get_locustag_comp_dict(options)

                target_model = run_augPhase(modelPrunedGPR, options)
            else:
                logging.warning("No EC_numbers found in input genome data")
                logging.warning("New reactions will NOT be added")

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
            logging.warning("No amino acid sequences found in input genome data")

    # Secondary metabolic modeling
    if options.smr_generation:
        if not options.eficaz:
            seq_records = get_target_genome_from_input(filetype, options)
        

        model_file = []
        files = glob.glob(options.outputfolder3 + os.sep + '*.xml')
        model_file = [each_file for each_file in files if '.xml' in each_file]

        if len(seq_records) == 1 and len(model_file) > 0 and '.xml' in model_file[0]:
            if options.total_region > 0:
                
                logging.info("Generating secondary metabolite biosynthesizing reactions..")
                logging.debug("Total number of regions: %s" %options.total_region)
                
            elif options.total_cluster > 0:
                
                logging.info("Generating secondary metabolite biosynthesizing reactions..")
                logging.debug("Total number of clusters: %s" %options.total_cluster)
                
            else:
                logging.warning("No cluster/region information found in input genome data")
                
            seq_record = seq_records[0]

            model_file = os.path.basename(model_file[0])
            target_model = cobra.io.read_sbml_model(
                           os.path.join(options.outputfolder3, model_file))

            target_model = run_secondary_modeling(seq_record, target_model, options)

            #target_model_no_gapsFilled = copy.deepcopy(target_model)

            get_target_nonprod_monomers_for_gapfilling(target_model, options)

            # NOTE: Disabled temporarily
            #target_model_complete = run_gapfilling(target_model, options)

            prune_unused_metabolites(target_model)

            runtime2 = time.strftime("Elapsed time %H:%M:%S",
                                     time.gmtime(time.time() - start))

            generate_outputs(options.outputfolder4,
                             runtime2, options,
                             #cobra_model_no_gapFilled = target_model_no_gapsFilled,
                             cobra_model = target_model)
            

        else:
            logging.warning("Secondary metabolic modeling not implemented;")

            if filetype == 'fasta':
                logging.warning("FASTA input file cannot be used for secondary modeling")
            elif len(seq_records) > 1:
                logging.warning(
                    "Input genome data with multiple records is currently not supported")
            elif len(model_file) == 0 or '.xml' not in model_file[0]:
                logging.warning("COBRA-compliant SBML file needed")
                

    remove_tmp_model_files(options)
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))

if __name__ == '__main__':
    main()
