
import datetime
import logging
import os


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
