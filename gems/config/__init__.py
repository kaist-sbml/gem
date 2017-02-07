
#Copyright 2017 BioInformatics Research Center, KAIST
#Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import ConfigParser
from argparse import Namespace
from os.path import join, abspath, dirname

_cfg_name = 'gems.cfg'

def load_config(options):

    _cfg_dir = join(dirname(abspath(__file__)), _cfg_name)

    config = ConfigParser.ConfigParser()
    config.read(_cfg_dir)

    for section in config.sections():
        if section not in options:
            # Create Namespace within 'options' namespace
            options.__dict__[section] = Namespace()

        for key, value in config.items(section):
            options.__dict__[section].__dict__[key] = value

