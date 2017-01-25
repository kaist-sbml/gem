
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import pickle
import pytest
from argparse import Namespace
from Bio import SeqIO
from cobra.io import read_sbml_model, write_sbml_model
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

    # Streptomyces collinus Tu 365
    seq_record = SeqIO.read(join(data_antismash_dir, 'NC_021985.1.final.gbk'), 'genbank')
    return seq_record


@pytest.fixture(scope="function")
def locustag_domain_dict():
    # Directly copied from raw data - keys not sorted
    locustag_domain_dict = {
        'B446_01685': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'ACP_DM02', 'PKS_KS_DM03'],
        'B446_01680': ['PKS_DH2_DM00', 'PKS_KR_DM01', 'ACP_DM02', 'PKS_KS_DM03', 'PKS_DH2_DM04', 'PKS_KR_DM05', 'cMT_DM06', 'ACP_DM07', 'PKS_KS_DM08', 'PKS_DHt_DM09', 'ACP_DM10', 'PKS_KS_DM11', 'Trans-AT_docking_DM12', 'PKS_KR_DM13', 'ACP_DM14', 'PKS_KS_DM15', 'Trans-AT_docking_DM16', 'PKS_KR_DM17', 'ACP_DM18', 'PKS_KS_DM19'],
        'B446_01675': ['PKS_DH2_DM00', 'ACP_DM01', 'PKS_KR_DM02', 'PKS_KS_DM03', 'PKS_DH_DM04', 'PKS_KR_DM05', 'ACP_DM06', 'PKS_Docking_Cterm_DM07'],
        'B446_01485': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02', 'Condensation_LCL_DM03', 'AMP-binding_DM04', 'PCP_DM05', 'Condensation_LCL_DM06', 'AMP-binding_DM07', 'PCP_DM08', 'Thioesterase_DM09'],
        'B446_01635': ['PKS_AT_DM00'],
        'B446_01480': ['Condensation_Starter_DM00', 'AMP-binding_DM01', 'PCP_DM02', 'Epimerization_DM03', 'Condensation_DCL_DM04', 'AMP-binding_DM05', 'PCP_DM06'],
        'B446_01525': ['AMP-binding_DM00'],
        'B446_01670': ['PKS_Docking_Nterm_DM00', 'PKS_KS_DM01', 'PKS_AT_DM02', 'PKS_DH_DM03', 'PKS_KR_DM04', 'ACP_DM05', 'PKS_KS_DM06', 'PKS_AT_DM07', 'ACP_DM08'],
        'B446_01730': ['NAD_binding_4_DM00'],
        'B446_01655': ['PKS_AT_DM00', 'PKS_AT_DM01'],
        'B446_01565': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02', 'Condensation_DCL_DM03', 'AMP-binding_DM04', 'PCP_DM05', 'Thioesterase_DM06'],
        'B446_01405': ['Polyketide_cyc_DM00'],
        'B446_01400': ['PKS_ER_DM00'],
        'B446_01695': ['ACP_DM00', 'PKS_KS_DM01', 'PKS_DHt_DM02', 'PKS_KR_DM03', 'ACP_DM04', 'PKS_KS_DM05'],
        'B446_01690': ['PKS_DHt_DM00', 'ACP_DM01', 'PKS_KR_DM02', 'PKS_KS_DM03', 'PKS_DHt_DM04', 'PKS_KR_DM05', 'cMT_DM06', 'ACP_DM07', 'PKS_KS_DM08', 'PKS_DHt_DM09', 'PKS_KR_DM10', 'PKS_KS_DM11', 'Trans-AT_docking_DM12', 'ACP_DM13'],
        'B446_01540': ['PCP_DM00', 'Epimerization_DM01'],
        'B446_01590': ['PKS_ER_DM00'],
        'B446_01700': ['ACPS_DM00'],
        'B446_01530': ['NRPS-COM_Nterm_DM00', 'Cglyc_DM01', 'AMP-binding_DM02', 'PCP_DM03'],
        'B446_01535': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02'],
        'B446_01640': ['oMT_DM00'],
        'B446_01660': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02']}
    return locustag_domain_dict


@pytest.fixture(scope="function")
def locustag_module_domain_dict():
    # Directly copied from raw data - keys not sorted
    custag_module_domain_dict = {
        'B446_01480_M00': ['Condensation_Starter_DM00', 'AMP-binding_DM01', 'PCP_DM02'],
        'B446_01480_M01': ['Epimerization_DM03', 'Condensation_DCL_DM04', 'AMP-binding_DM05', 'PCP_DM06'],
        'B446_01700_M00': ['ACPS_DM00'],
        'B446_01675_M02': ['PKS_Docking_Cterm_DM07'],
        'B446_01675_M00': ['PKS_DH2_DM00', 'ACP_DM01'],
        'B446_01675_M01': ['PKS_KR_DM02', 'PKS_KS_DM03', 'PKS_DH_DM04', 'PKS_KR_DM05', 'ACP_DM06'],
        'B446_01690_M01': ['PKS_KR_DM02', 'PKS_KS_DM03', 'PKS_DHt_DM04', 'PKS_KR_DM05', 'cMT_DM06', 'ACP_DM07'],
        'B446_01730_M00': ['NAD_binding_4_DM00'],
        'B446_01690_M02': ['PKS_KS_DM08', 'PKS_DHt_DM09', 'PKS_KR_DM10', 'PKS_KS_DM11', 'Trans-AT_docking_DM12', 'ACP_DM13'],
        'B446_01695_M00': ['ACP_DM00'],
        'B446_01400_M00': ['PKS_ER_DM00'],
        'B446_01695_M01': ['PKS_KS_DM01', 'PKS_DHt_DM02', 'PKS_KR_DM03', 'ACP_DM04'],
        'B446_01695_M02': ['PKS_KS_DM05'],
        'B446_01690_M00': ['PKS_DHt_DM00', 'ACP_DM01'],
        'B446_01635_M00': ['PKS_AT_DM00'],
        'B446_01540_M01': ['Epimerization_DM01'],
        'B446_01405_M00': ['Polyketide_cyc_DM00'],
        'B446_01655_M00': ['PKS_AT_DM00', 'PKS_AT_DM01'],
        'B446_01565_M02': ['Thioesterase_DM06'],
        'B446_01565_M00': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02'],
        'B446_01565_M01': ['Condensation_DCL_DM03', 'AMP-binding_DM04', 'PCP_DM05'],
        'B446_01535_M00': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02'],
        'B446_01530_M00': ['NRPS-COM_Nterm_DM00', 'Cglyc_DM01', 'AMP-binding_DM02', 'PCP_DM03'],
        'B446_01540_M00': ['PCP_DM00'],
        'B446_01525_M00': ['AMP-binding_DM00'],
        'B446_01485_M01': ['Condensation_LCL_DM03', 'AMP-binding_DM04', 'PCP_DM05'],
        'B446_01485_M00': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02'],
        'B446_01485_M03': ['Thioesterase_DM09'],
        'B446_01485_M02': ['Condensation_LCL_DM06', 'AMP-binding_DM07', 'PCP_DM08'],
        'B446_01640_M00': ['oMT_DM00'],
        'B446_01685_M01': ['PKS_KS_DM03'],
        'B446_01685_M00': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'ACP_DM02'],
        'B446_01590_M00': ['PKS_ER_DM00'],
        'B446_01670_M01': ['PKS_KS_DM06', 'PKS_AT_DM07', 'ACP_DM08'],
        'B446_01670_M00': ['PKS_Docking_Nterm_DM00', 'PKS_KS_DM01', 'PKS_AT_DM02', 'PKS_DH_DM03', 'PKS_KR_DM04', 'ACP_DM05'],
        'B446_01680_M02': ['PKS_KS_DM08', 'PKS_DHt_DM09', 'ACP_DM10'],
        'B446_01680_M03': ['PKS_KS_DM11', 'Trans-AT_docking_DM12', 'PKS_KR_DM13', 'ACP_DM14'],
        'B446_01680_M00': ['PKS_DH2_DM00', 'PKS_KR_DM01', 'ACP_DM02'],
        'B446_01680_M01': ['PKS_KS_DM03', 'PKS_DH2_DM04', 'PKS_KR_DM05', 'cMT_DM06', 'ACP_DM07'],
        'B446_01680_M04': ['PKS_KS_DM15', 'Trans-AT_docking_DM16', 'PKS_KR_DM17', 'ACP_DM18'],
        'B446_01680_M05': ['PKS_KS_DM19'],
        'B446_01660_M00': ['Condensation_LCL_DM00', 'AMP-binding_DM01', 'PCP_DM02']}
    return locustag_module_domain_dict


