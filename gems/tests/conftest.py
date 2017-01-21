
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from argparse import Namespace
import pytest
from cobra.io import read_sbml_model, write_sbml_model
from os.path import join, abspath, dirname

data_dir = join(dirname(abspath(__file__)), 'data')

@pytest.fixture(scope="function")
def model():
    model = read_sbml_model(join(data_dir, 'iMK1208Edited4.xml'))

    # Returning model is necessary, otherwise error occurs
    return model

@pytest.fixture(scope="function")
def options():
    options = Namespace()

    return options

