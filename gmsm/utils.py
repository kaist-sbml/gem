
import cobra
import datetime
import logging
import os
import sys
import subprocess
from os.path import getmtime, isfile, join, split


def setup_logging(run_ns):
    if run_ns.verbose:
        log_level = logging.INFO
    elif run_ns.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING

    logging.basicConfig(format='%(levelname)s: %(message)s', level=log_level)

    if run_ns.verbose or run_ns.debug:
        logger = logging.getLogger('')
        fomatter = logging.Formatter(
                '[%(levelname)s|%(filename)s:%(lineno)s] > %(message)s')
        fh = logging.FileHandler(
                os.path.join(run_ns.outputfolder, 'gmsm.log'), mode = 'w')
        fh.setFormatter(fomatter)
        logger.setLevel(log_level)
        logger.addHandler(fh)


def get_version():
    import gmsm
    version = gmsm.__version__

    return version


def get_git_log():
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        out, err, retcode = execute(args)
        return out.strip()
    except OSError:
        pass
    return""


def check_input_options(run_ns):
    if not run_ns.input:
        logging.warning("Provide input file via ('-i')")
        sys.exit(1)

    if not run_ns.ec_file and \
            not run_ns.eficaz and \
            not run_ns.pmr_generation and \
            not run_ns.smr_generation and \
            not run_ns.comp:
                logging.warning("Select one of the options: '-E', '-p' or '-s'")
                sys.exit(1)

    if run_ns.comp:
        if not run_ns.pmr_generation:
            logging.warning(
                    "Primary metabolic modeling option ('-p') should also be selected")
            sys.exit(1)

    if run_ns.ec_file:
        if run_ns.eficaz:
            logging.warning(
                    "EC number file option ('-e') or the EFICAz run option ('-E') should be removed")
            sys.exit(1)

        if not run_ns.pmr_generation:
            logging.warning(
                    "Primary metabolic modeling option ('-p') should also be selected")
            sys.exit(1)


# Adopted from antismash.utils
def locate_executable(name):
    "Find an executable in the path and return the full path"
    # In windows, executables tend to end on .exe
    if sys.platform == 'win32':
        name += ".exe"
    file_path, _ = split(name)
    if file_path != "":
        if isfile(name) and os.access(name, os.X_OK):
            logging.debug("Found executable %r", name)
            return name
    for p in os.environ["PATH"].split(os.pathsep):
        full_name = join(p, name)
        if isfile(full_name) and os.access(full_name, os.X_OK):
            logging.debug("Found executable %r", full_name)
            return full_name

    return None


# Adopted from antismash.utils
# Ignore the pylint warning about input being redifined, as we're just
# following the subprocess names here.
# pylint: disable=redefined-builtin
def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as e:
         logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, e)
         raise


# Adopted from antismash.utils
def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
         # force into a tuple
         types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features


# Adopted from antismash.utils
def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")


# Adopted from antismash.utils
def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"


# For regular update of the cache:
# KEGG updates its EC_number:reaction ID pairs time to time,
#and old caches can cause an error
def time_bomb(cache_file, config_ns):
    today = datetime.datetime.today()
    modified_date = datetime.datetime.fromtimestamp(getmtime(cache_file))
    file_age = today - modified_date

    if int(file_age.days) > int(config_ns.utils.time_bomb_duration):
        logging.debug('File %s is older than %s days (currently %s days)',
                cache_file, config_ns.utils.time_bomb_duration, file_age.days)
        os.remove(cache_file)
        logging.debug('File %s was removed', cache_file)
    else:
        logging.debug('File %s has not reached %s days (currently %s days)',
                cache_file, config_ns.utils.time_bomb_duration, file_age.days)


def get_keggid_from_mnxr(mnxr, io_ns, primary_model_ns):
    if len(io_ns.mnxr_kegg_dict[mnxr]) > 1:
        keggid_list = []

        for keggid in io_ns.mnxr_kegg_dict[mnxr]:
            if keggid in primary_model_ns.rxnid_info_dict:
                keggid_list.append(keggid)

        if len(keggid_list) == 1:
            kegg_id = keggid_list[0]
        # Choose KEGG reaction ID with a greater value for multiple KEGG IDs given to MNXR
        elif len(keggid_list) > 1:
            keggid_list.sort()
            kegg_id = keggid_list[-1]

    elif len(io_ns.mnxr_kegg_dict[mnxr]) == 1:
        kegg_id = io_ns.mnxr_kegg_dict[mnxr][0]

    return kegg_id


def check_duplicate_rxn(model, rxn2):

    duplicate_rxn = []

    for j in range(len(model.reactions)):
        rxn1 = model.reactions[j]

        comparison = compare_rxns(rxn1, rxn2)

        if comparison == 'same':
            duplicate_rxn.append(rxn1.id)

    if len(duplicate_rxn) >= 1:
        return 'duplicate'
    elif len(duplicate_rxn) == 0:
        return 'unique'


def compare_rxns(rxn1, rxn2):

    rxn1_metab_dict = {}
    for i in range(len(rxn1.metabolites)):
        rxn1_metab_dict[str(list(rxn1.metabolites.keys())[i])] = \
                    float(rxn1.metabolites[list(rxn1.metabolites.keys())[i]])

    rxn2_metab_dict = {}
    for i in range(len(rxn2.metabolites)):
        rxn2_metab_dict[str(list(rxn2.metabolites.keys())[i])] = \
                    float(rxn2.metabolites[list(rxn2.metabolites.keys())[i]])

    if rxn1_metab_dict.items() == rxn2_metab_dict.items():

        if rxn1.reversibility == rxn2.reversibility:
            return 'same'
        else:
            return 'different'
    else:
        return 'different'


#'add_reaction' requires writing/reloading of the model
def stabilize_model(model, folder, label, diff_name=False):

    if diff_name:
        model_name = '%s.xml' %label
    else:
        if label:
            model_name = 'model_%s.xml' %label
        else:
            model_name = 'model.xml'

    for rxn in model.reactions:
        if rxn.gene_reaction_rule == "()":
            rxn.gene_reaction_rule = ""
                
    cobra.io.write_sbml_model(model, join('%s' %folder, model_name))
    model = cobra.io.read_sbml_model(join('%s' %folder, model_name))

    return model


#Output: a dictionary file for major Exchange reactions {Exchange reaction ID:flux value}
def get_exrxnid_flux(model, template_exrxnid_flux_dict):

    target_exrxnid_flux_dict = {}
    model.optimize()

    for exrxn_id in template_exrxnid_flux_dict:
        if exrxn_id in model.reactions:
            rxn = model.reactions.get_by_id(exrxn_id)
            target_exrxnid_flux_dict[exrxn_id] = rxn.flux
        else:
            continue
    return target_exrxnid_flux_dict


#Output: a list file having either T or F for major Exchange reactions
def check_exrxn_flux_direction(
        template_exrxnid_flux_dict, target_exrxnid_flux_dict, config_ns):

    exrxn_flux_change_list = []

    for exrxn_id in template_exrxnid_flux_dict:
        if exrxn_id in target_exrxnid_flux_dict:
            template_exrxn_flux = template_exrxnid_flux_dict[exrxn_id]
            target_exrxn_flux = target_exrxnid_flux_dict[exrxn_id]

            if float(template_exrxn_flux) == 0:
                if target_exrxn_flux == 0:
                    ratio_exrxn_flux = 1
                    logging.debug("%s: from zero to zero flux", exrxn_id)

                else:
                    ratio_exrxn_flux = False
                    logging.debug("%s: from zero to non-zero flux", exrxn_id)
                    continue

            else:
                ratio_exrxn_flux = float(target_exrxn_flux)/float(template_exrxn_flux)

            #Similar species are allowed to uptake nutrients within a decent range
            if ratio_exrxn_flux > 0 and ratio_exrxn_flux < float(config_ns.cobrapy.nutrient_uptake_rate) \
                and 1/ratio_exrxn_flux < float(config_ns.cobrapy.nutrient_uptake_rate):
                exrxn_flux_change_list.append('T')

            #Cause drastic changes in Exchange reaction fluxes
            #(direction and/or magnitude)
            else:
                logging.debug("Drastic change occured in %s: %f -> %f" \
                              %(exrxn_id, float(template_exrxn_flux), target_exrxn_flux))
                exrxn_flux_change_list.append('F')

    return exrxn_flux_change_list
