
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import pickle
import pytest
from argparse import Namespace
from Bio import SeqIO
from cobra.io import read_sbml_model, write_sbml_model
from gems.secondary_model.sec_met_rxn_generation import get_cluster_location
from os.path import join, abspath, dirname

data_model_dir = join(dirname(abspath(__file__)), 'data_model')
data_antismash_dir = join(dirname(abspath(__file__)), 'data_antismash')

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


@pytest.fixture(scope="function")
def seq_record():

    options = Namespace()

    # Streptomyces collinus Tu 365
    seq_record = SeqIO.read(join(data_antismash_dir, 'NC_021985.1.final.gbk'), 'genbank')
    options.seq_record = seq_record

    # Hybrid cluster: nrps-transatpks-t1pks
    # locations: 341017 - 503094
    # Kirromycin biosynthetic gene cluster (81% of genes show similarity)
    cluster_nr = 3
    get_cluster_location(cluster_nr, options)

    return options
