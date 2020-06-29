
#Copyright 2017 BioInformatics Research Center, KAIST
#Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import configparser
from argparse import Namespace
from os.path import join, abspath, dirname

_cfg_name = 'gmsm.cfg'

def load_config(config_ns):

    _cfg_dir = join(dirname(abspath(__file__)), _cfg_name)

    config = configparser.ConfigParser()
    config.read(_cfg_dir)

    for section in config.sections():
        if section not in config_ns:
            # Create Namespace within 'config_ns' namespace
            config_ns.__dict__[section] = Namespace()

        for key, value in config.items(section):
            config_ns.__dict__[section].__dict__[key] = value

