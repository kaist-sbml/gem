
import cobra
import datetime
import logging
import os
import sys
import subprocess
from os.path import getmtime, isfile, join, split


def setup_logging(options):
    if options.verbose:
        log_level = logging.INFO
    elif options.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING

    logging.basicConfig(format='%(levelname)s: %(message)s', level=log_level)

def get_version():
    import gems
    version = gems.__version__

    return version

def get_git_log():
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        out, err, retcode = execute(args)
        return out.strip()
    except OSError:
        pass
    return""

def setup_logfile_format(options):
    if options.debug:
        logger = logging.getLogger('')
        fomatter = logging.Formatter(
                '[%(levelname)s|%(filename)s:%(lineno)s] > %(message)s')
        fh = logging.FileHandler(
                os.path.join(options.outputfolder, 'gems.log'), mode = 'w')
        fh.setFormatter(fomatter)
        logger.addHandler(fh)


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
    except OSError, e:
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
def time_bomb(cache_file, options):
    today = datetime.datetime.today()
    modified_date = datetime.datetime.fromtimestamp(getmtime(cache_file))
    file_age = today - modified_date

    if int(file_age.days) > int(options.utils.time_bomb_duration):
        logging.debug('File %s is older than %s days (currently %s days)',
                cache_file, options.utils.time_bomb_duration, file_age.days)
        os.remove(cache_file)
        logging.debug('File %s was removed', cache_file)
    else:
        logging.debug('File %s has not reached %s days (currently %s days)',
                cache_file, options.utils.time_bomb_duration, file_age.days)


def get_keggid_from_mnxr(mnxr, options):
    if len(options.mnxr_kegg_dict[mnxr]) > 1:
        keggid_list = []

        for keggid in options.mnxr_kegg_dict[mnxr]:
            if keggid in options.rxnid_info_dict:
                keggid_list.append(keggid)

        if len(keggid_list) == 1:
            kegg_id = keggid_list[0]
        # Choose KEGG reaction ID with a greater value for multiple KEGG IDs given to MNXR
        elif len(keggid_list) > 1:
            keggid_list.sort()
            kegg_id = keggid_list[-1]

    elif len(options.mnxr_kegg_dict[mnxr]) == 1:
        kegg_id = options.mnxr_kegg_dict[mnxr][0]

    return kegg_id


#'add_reaction' requires writing/reloading of the model
def stabilize_model(model, folder, label, diff_name=False):

    if diff_name:
        model_name = '%s.xml' %label
    else:
        if label:
            model_name = 'model_%s.xml' %label
        else:
            model_name = 'model.xml'

    cobra.io.write_sbml_model(model, join('%s' %folder, model_name), use_fbc_package=False)
    model = cobra.io.read_sbml_model(join('%s' %folder, model_name))

    return model


#Output: a dictionary file for major Exchange reactions {Exchange reaction ID:flux value}
def get_exrxnid_flux(model, template_exrxnid_flux_dict):

    target_exrxnid_flux_dict = {}
    model.optimize()

    for exrxn_id in template_exrxnid_flux_dict:
        if exrxn_id in model.solution.x_dict:
            target_exrxnid_flux_dict[exrxn_id] = model.solution.x_dict[exrxn_id]
        else:
            continue
    return target_exrxnid_flux_dict


#Output: a list file having either T or F for major Exchange reactions
def check_exrxn_flux_direction(
        template_exrxnid_flux_dict, target_exrxnid_flux_dict, options):

    exrxn_flux_change_list = []

    for exrxn_id in template_exrxnid_flux_dict:
        if exrxn_id in target_exrxnid_flux_dict:
            template_exrxn_flux = template_exrxnid_flux_dict[exrxn_id]
            target_exrxn_flux = target_exrxnid_flux_dict[exrxn_id]

            if float(template_exrxn_flux) != 0:
                ratio_exrxn_flux = float(target_exrxn_flux)/float(template_exrxn_flux)
            else:
                logging.debug("%s has a zero flux", exrxn_id)

            #Similar species are allowed to uptake nutrients within a decent range
            if float(target_exrxn_flux)*float(template_exrxn_flux) > float(options.cobrapy.non_zero_flux_cutoff) \
                    and float(options.cobrapy.nutrient_uptake_rate) < ratio_exrxn_flux \
                    and ratio_exrxn_flux < float(options.cobrapy.nutrient_uptake_rate):
                exrxn_flux_change_list.append('T')

            #Cause drastic changes in Exchange reaction fluxes
            #(direction and/or magnitude)
            else:
                exrxn_flux_change_list.append('F')

    return exrxn_flux_change_list
