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
    get_ec_file,
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
            '\rGMSM version {version} ({git_log})\n\nusage: {usage}\n----------------------------------------------------------------------------------'\
            .format(version=utils.get_version(),
                    git_log=utils.get_git_log(),
                    usage='run_gmsm.py [-h] [Resource management] [Input and output setting] [GMSM modeling options] [Debugging and logging options]')
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
                        choices=['bsu', 'clj', 'cre', 'eco', 'hpy', 'mtu', 'nsal', 'ppu', 'sce', 'sco'],
                        help=
                        "'bsu': iYO844 (17573341); Bacillus subtilis subsp. subtilis str. 168\n"
                        #NOTE: metabolite compartments (other than 'c') NOT standardized
                        "'clj': iHN637 (24274140); Clostridium ljungdahlii DSM 13528\n"
                        "'cre': iCre1355 (26485611); Chlamydomonas reinhardtii\n"
                        "'eco': iML1515 (29020004); Escherichia coli str. K-12 substr. MG1655\n"
                        "'hpy': iIT341 (16077130); Helicobacter pylori 26695\n"
                        "'mtu': iNJ661 (17555602); Mycobacterium tuberculosis H37Rv\n"
                        "'nsal': iNS934 (28676050); Nannochloropsis salina\n"
                        "'ppu': iJN746 (18793442); Pseudomonas putida KT2440\n"
                        "'sce': iMM904 (19321003); Saccharomyces cerevisiae S288C\n"
                        "'sco': iKS1317 (30525286); Streptomyces coelicolor A3(2)"
                        )

    group = parser.add_argument_group('GMSM modeling options',
                        "At least one of the three options should be selected:"
                        " '-E', '-p' and '-s'\n"
                        "Primary metabolic modeling option ('-p') should be selected "
                        "when using '-e' and/or '-C' options\n"
                        " - Examples:\n"
                        "   '-e -E': NOT acceptable\n"
                        "   '-p -e': Acceptable\n"
                        "   '-s -p -e': Acceptable\n"
                        "   '-s -e': NOT acceptable"
                        )
    group.add_argument('-e', '--ec',
                        dest='ec_file',
                        default=False,
                        help="EC number prediction file")
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
                        dest='eficaz',
                        action='store_true',
                        default=False,
                        help="Run EC number prediction using EFICAz")
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

    # Create a namespace for each modules 
    # *_ns are arranged based on appearing orders in this python file
    run_ns = parser.parse_args()
    
    [io_ns, config_ns, homology_ns, primary_model_ns, \
     secondary_model_ns] = [run_ns for module in range(5)]
  
    # Create an output directory for a log file
    make_folder(run_ns.outputfolder)

    utils.setup_logging(run_ns)

    # Create output folders
    setup_outputfolders(run_ns, io_ns)

    if run_ns.version:
        print('GMSM version %s (%s)' %(utils.get_version(), utils.get_git_log()))
        sys.exit(0)

    utils.check_input_options(run_ns)

    # Warning messages from cobrapy turned off by default
    if not run_ns.warning:
        warnings.filterwarnings("ignore")

    logging.info('Starting GMSM ver. %s (%s)', utils.get_version(), utils.get_git_log())

    show_input_options(run_ns)

    logging.info("Reading input genome files..")
    filetype = check_input_filetype(run_ns)

    # Load config data
    load_config(config_ns)

    # Check prerequisites of executables and libraries
    check_prereqs(run_ns)

    # EC number prediction
    if run_ns.eficaz:
        seq_records = get_target_genome_from_input(filetype, run_ns, io_ns)

        if run_ns.eficaz_path and \
                io_ns.targetGenome_locusTag_aaSeq_dict and \
                not run_ns.ec_file:

            if filetype == 'fasta' or len(seq_records) > 1:
                logging.info("Input file in FASTA format or with multiple records:")
                logging.info("Raw EFICAz output (.txt)  will be generated, not GenBank")

            if len(seq_records) == 1:
                getECs1(run_ns, io_ns, seq_record = seq_records[0])
            elif len(seq_records) > 1:
                getECs2(run_ns, io_ns)
        else:
            logging.warning("EFICAz not implemented;")

            if not run_ns.eficaz_path:
                logging.warning("EFICAz not found")
            elif not io_ns.targetGenome_locusTag_aaSeq_dict:
                logging.warning(
                        "No amino acid sequences found in input genome data")

    # Primary metabolic modeling
    if run_ns.pmr_generation:
        if not run_ns.eficaz:
            get_target_genome_from_input(filetype, run_ns, io_ns)

        if run_ns.ec_file:
            get_ec_file(run_ns, io_ns)

        get_fasta_files(run_ns, io_ns)

        if io_ns.targetGenome_locusTag_aaSeq_dict:
            get_homologs(io_ns, homology_ns)
            model = get_pickles_prunPhase(io_ns)
            modelPrunedGPR = run_prunPhase(model, io_ns, config_ns, homology_ns, primary_model_ns)

            if io_ns.targetGenome_locusTag_ec_dict:
                get_pickles_augPhase(io_ns)

                if run_ns.comp:
                    get_locustag_comp_dict(run_ns, io_ns)

                target_model = run_augPhase(modelPrunedGPR, run_ns, io_ns, config_ns, homology_ns, primary_model_ns)
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
                generate_outputs(io_ns.outputfolder3, runtime1, run_ns, io_ns, homology_ns, primary_model_ns, secondary_model_ns, cobra_model = target_model)
            except:
                generate_outputs(io_ns.outputfolder3, runtime1, run_ns, io_ns, homology_ns, primary_model_ns, secondary_model_ns, cobra_model = modelPrunedGPR)
        else:
            logging.warning("Primary metabolic modeling not implemented;")
            logging.warning("No amino acid sequences found in input genome data")

    # Secondary metabolic modeling
    if run_ns.smr_generation:
        if not run_ns.eficaz:
            get_target_genome_from_input(filetype, run_ns, io_ns)

        model_file = []
        files = glob.glob(io_ns.outputfolder3 + os.sep + '*.xml')
        model_file = [each_file for each_file in files if '.xml' in each_file]

        if len(model_file) > 0 and '.xml' in model_file[0] \
            and (io_ns.total_region > 0 or io_ns.total_cluster > 0):

            if io_ns.total_region > 0:
                logging.info("Generating secondary metabolite biosynthesizing reactions..")
                logging.debug("Total number of regions: %s" %io_ns.total_region)

            elif io_ns.total_cluster > 0:
                logging.info("Generating secondary metabolite biosynthesizing reactions..")
                logging.debug("Total number of clusters: %s" %io_ns.total_cluster)

            model_file = os.path.basename(model_file[0])
            target_model = cobra.io.read_sbml_model(
                           os.path.join(io_ns.outputfolder3, model_file))

            target_model = run_secondary_modeling(target_model, io_ns, config_ns, secondary_model_ns)

            #target_model_no_gapsFilled = copy.deepcopy(target_model)

            get_target_nonprod_monomers_for_gapfilling(target_model, io_ns, config_ns, secondary_model_ns)

            # NOTE: Disabled temporarily
            #target_model_complete = run_gapfilling(target_model, io_ns, config_ns, secondary_model_ns)

            prune_unused_metabolites(target_model)

            runtime2 = time.strftime("Elapsed time %H:%M:%S",
                                     time.gmtime(time.time() - start))

            generate_outputs(io_ns.outputfolder4,
                             runtime2, run_ns, io_ns, homology_ns, primary_model_ns, secondary_model_ns,
                             #cobra_model_no_gapFilled = target_model_no_gapsFilled,
                             cobra_model = target_model)

        else:
            logging.warning("Secondary metabolic modeling not implemented;")

            if filetype == 'fasta':
                logging.warning("FASTA input file cannot be used for secondary modeling")
            elif len(model_file) == 0 or '.xml' not in model_file[0]:
                logging.warning("COBRA-compliant SBML file needed")
            elif io_ns.total_region == 0 and io_ns.total_cluster == 0:
                logging.warning("No cluster/region information found in input genome data")

    remove_tmp_model_files(io_ns)
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))

if __name__ == '__main__':
    main()
