
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import pickle
import pytest
from argparse import Namespace
from cobra.io import read_sbml_model, write_sbml_model
from os.path import join, abspath, dirname

data_model_dir = join(dirname(abspath(__file__)), 'data_model')

@pytest.fixture(scope="function")
def model():
    model = read_sbml_model(join(data_model_dir, 'iMK1208Edited4.xml'))

    # Returning model is necessary, otherwise error occurs
    return model

@pytest.fixture(scope="function")
def options():
    options = Namespace()

    template_exrxnid_flux_dict = pickle.load(
             open(join(data_model_dir, 'sco_tempModel_exrxnid_flux_dict.p'),'rb'))
    options.template_exrxnid_flux_dict = template_exrxnid_flux_dict

    return options

