
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from gems.primary_model import (
    augPhase_utils,
    prunPhase_utils
)

class TestPrimary_model:
    """Test functions in gems.primary_model"""

    def test_add_nonBBH_rxn(self, model):
        rxnid_to_add_list = []
        mnxr_to_add_list = []
        rxnid_info_dict = []
        # try three 3 cases having a mixture of different types of metabolite IDs
        rxnid_mnxm_coeff_dict = {
            'R00618':{'thmtp': 1.0, 'h': -1.0, 'h2o': 1.0,'pi': -1.0, 'thmpp': -1.0},
            'R00125':{'h': -2.0, 'ap4a': 1.0, 'h2o': 1.0, 'adp': -2.0},
            'R01968':{'dgmp': -1.0, 'h2o': -1.0, 'pi': 1.0, 'dgsn': 1.0}}

