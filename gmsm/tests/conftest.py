
import pickle
import pytest
from argparse import Namespace
from Bio import SeqIO
from cobra.io import read_sbml_model, write_sbml_model
from os.path import join, abspath, dirname

data_model_dir = join(dirname(abspath(__file__)), 'data')
data_antismash_dir = join(dirname(abspath(__file__)), 'data_antismash')

@pytest.fixture(scope="function")
def sco_tmp_model():

    # Streptomyces coelicolor A3(2)
    model = read_sbml_model(join(data_model_dir, 'iMK1208Edited4.xml'))

    # Returning model is necessary, otherwise error occurs
    return model


@pytest.fixture(scope="function")
def sci_primary_model():

    # Streptomyces collinus Tu 365
    model = read_sbml_model(join(data_model_dir, 'sci_primary_model.xml'))

    # Returning model is necessary, otherwise error occurs
    return model


@pytest.fixture(scope="function")
def mnxref():

    # MNXref SBML model
    model = pickle.load(open(join(data_model_dir, 'MNXref.p'), 'rb'))

    return model


@pytest.fixture(scope="function")
def options():
    options = Namespace()
    return options


@pytest.fixture(scope="function")
def bbh_dict():
    temp_target_BBH_dict = pickle.load(
            open(join(data_model_dir, 'sco_sci_temp_target_BBH_dict.p'),'rb'))
    return temp_target_BBH_dict


@pytest.fixture(scope="function")
def mnxr_kegg_dict():
    mnxr_kegg_dict = pickle.load(
            open(join(data_model_dir, 'mnxr_kegg_dict.p'),'rb'))
    return mnxr_kegg_dict


@pytest.fixture(scope="function")
def sco_tmp_model_flux():
    template_exrxnid_flux_dict = pickle.load(
            open(join(data_model_dir, 'sco_tempModel_exrxnid_flux_dict.p'),'rb'))
    return template_exrxnid_flux_dict


@pytest.fixture(scope="function")
def seq_record():

    # Streptomyces collinus Tu 365
    seq_record = SeqIO.read(join(data_antismash_dir, 'NC_021985.1.final.gbk'), 'genbank')
    return seq_record


@pytest.fixture(scope="function")
def eficaz_file():
    eficaz_file = join(data_model_dir, 'NSK_all_genomes_ec_test.txt')
    return eficaz_file


@pytest.fixture(scope="function")
def comp_file():
    comp_file = join(data_model_dir, 'Nanno_Compartment_result_dic_v3_test.txt')
    return comp_file


@pytest.fixture(scope="function")
def locustag_domain_dict():
    # Directly copied from raw data - keys not sorted
    locustag_domain_dict = {
        'B446_01685': ['Condensation_LCL_D00', 'AMP-binding_D01', 'ACP_D02', 'PKS_KS_D03'],
        'B446_01680': ['PKS_DH2_D00', 'PKS_KR_D01', 'ACP_D02', 'PKS_KS_D03', 'PKS_DH2_D04', 'PKS_KR_D05', 'cMT_D06', 'ACP_D07', 'PKS_KS_D08', 'PKS_DHt_D09', 'ACP_D10', 'PKS_KS_D11', 'Trans-AT_docking_D12', 'PKS_KR_D13', 'ACP_D14', 'PKS_KS_D15', 'Trans-AT_docking_D16', 'PKS_KR_D17', 'ACP_D18', 'PKS_KS_D19'],
        'B446_01675': ['PKS_DH2_D00', 'ACP_D01', 'PKS_KR_D02', 'PKS_KS_D03', 'PKS_DH_D04', 'PKS_KR_D05', 'ACP_D06', 'PKS_Docking_Cterm_D07'],
        'B446_01485': ['Condensation_LCL_D00', 'AMP-binding_D01', 'PCP_D02', 'Condensation_LCL_D03', 'AMP-binding_D04', 'PCP_D05', 'Condensation_LCL_D06', 'AMP-binding_D07', 'PCP_D08', 'Thioesterase_D09'],
        'B446_01635': ['PKS_AT_D00'],
        'B446_01480': ['Condensation_Starter_D00', 'AMP-binding_D01', 'PCP_D02', 'Epimerization_D03', 'Condensation_DCL_D04', 'AMP-binding_D05', 'PCP_D06'],
        'B446_01525': ['AMP-binding_D00'],
        'B446_01670': ['PKS_Docking_Nterm_D00', 'PKS_KS_D01', 'PKS_AT_D02', 'PKS_DH_D03', 'PKS_KR_D04', 'ACP_D05', 'PKS_KS_D06', 'PKS_AT_D07', 'ACP_D08'],
        'B446_01730': ['NAD_binding_4_D00'],
        'B446_01655': ['PKS_AT_D00', 'PKS_AT_D01'],
        'B446_01565': ['Condensation_LCL_D00', 'AMP-binding_D01', 'PCP_D02', 'Condensation_DCL_D03', 'AMP-binding_D04', 'PCP_D05', 'Thioesterase_D06'],
        'B446_01405': ['Polyketide_cyc_D00'],
        'B446_01400': ['PKS_ER_D00'],
        'B446_01695': ['ACP_D00', 'PKS_KS_D01', 'PKS_DHt_D02', 'PKS_KR_D03', 'ACP_D04', 'PKS_KS_D05'],
        'B446_01690': ['PKS_DHt_D00', 'ACP_D01', 'PKS_KR_D02', 'PKS_KS_D03', 'PKS_DHt_D04', 'PKS_KR_D05', 'cMT_D06', 'ACP_D07', 'PKS_KS_D08', 'PKS_DHt_D09', 'PKS_KR_D10', 'PKS_KS_D11', 'Trans-AT_docking_D12', 'ACP_D13'],
        'B446_01540': ['PCP_D00', 'Epimerization_D01'],
        'B446_01590': ['PKS_ER_D00'],
        'B446_01700': ['ACPS_D00'],
        'B446_01530': ['NRPS-COM_Nterm_D00', 'Cglyc_D01', 'AMP-binding_D02', 'PCP_D03'],
        'B446_01535': ['Condensation_LCL_D00', 'AMP-binding_D01', 'PCP_D02'],
        'B446_01640': ['oMT_D00'],
        'B446_01660': ['Condensation_LCL_D00', 'AMP-binding_D01', 'PCP_D02']}
    return locustag_domain_dict


