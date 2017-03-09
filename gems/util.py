
import cobra
import datetime
import logging
import os
from os.path import join

# For regular update of the cache:
# KEGG updates its EC_number:reaction ID pairs time to time,
#and old caches can cause an error
def time_bomb(cache_file, options):
    today = datetime.datetime.today()
    modified_date = datetime.datetime.fromtimestamp(os.path.getmtime(cache_file))
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
