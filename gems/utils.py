
import cobra
import datetime
import logging
import os
import sys
import subprocess
from os.path import getmtime, isfile, join, split

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


#'add_reaction' requires writing/reloading of the model
def stabilize_model(model, label, options):
    cobra.io.write_sbml_model(
                model,
                join('%s' %options.outputfolder5, 'model_%s.xml' %label),
                use_fbc_package=False)
    model = cobra.io.read_sbml_model(
                join('%s' %options.outputfolder5, 'model_%s.xml' %label)
                )
    return model
