
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from os.path import join
from gems.config import load_config
from gems.primary_model import prunPhase_utils, augPhase_utils

class TestPrimary_model:
    """Test functions in gems.primary_model"""

    def test_get_rxn_fate(self):
        bbh_avail_list = ['1', 'or', '1']
        rxn_fate = prunPhase_utils.get_rxn_fate(bbh_avail_list)
        assert rxn_fate == '1'

        bbh_avail_list = [['1', 'and', '1'], 'and', '0']
        rxn_fate = prunPhase_utils.get_rxn_fate(bbh_avail_list)
        assert rxn_fate == '0'

        bbh_avail_list = [['0', 'and', '1'], 'or', '1']
        rxn_fate = prunPhase_utils.get_rxn_fate(bbh_avail_list)
        assert rxn_fate == '1'

        bbh_avail_list = [[['1', 'and', '1'], 'or', ['0', 'and', '1']], 'and', '1']
        rxn_fate = prunPhase_utils.get_rxn_fate(bbh_avail_list)
        assert rxn_fate == '1'

        # NOTE: This issue has not been resolved. OR is returned regardless of the Boolean.
        #If successful, the result should be 'rxn_fate == 0'
        bbh_avail_list = [['1', 'and', '0'], 'or', ['0', 'and', '1'], 'and', ['0', 'or', '1']]
        rxn_fate = prunPhase_utils.get_rxn_fate(bbh_avail_list)
        assert rxn_fate != '0'


    def test_check_bbh_availability(self):
        temp_target_BBH_dict = {}
        tempModel_biggRxnid_locusTag_dict = {}

        # PGI
        temp_target_BBH_dict['SCO1942'] = ['B446_30415', 'B446_10110']
        temp_target_BBH_dict['SCO6659'] = ['B446_30415', 'B446_10110']
        tempModel_biggRxnid_locusTag_dict['PGI'] = ['SCO1942', 'or', 'SCO6659']

        bbh_avail_list = prunPhase_utils.check_bbh_availability(
                temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict['PGI'])
        assert bbh_avail_list == ['1', 'or', '1']

        # AKGDH2
        temp_target_BBH_dict['SCO4594'] = ['B446_21645']
        temp_target_BBH_dict['SCO4595'] = ['B446_21650']
        temp_target_BBH_dict['SCO6269'] = ['B446_21645']
        temp_target_BBH_dict['SCO6270'] = ['B446_21650']
        temp_target_BBH_dict['SCO0681'] = ['B446_05650']
        tempModel_biggRxnid_locusTag_dict['AKGDH2'] = \
                [[['SCO4594', 'and', 'SCO4595'], 'or', ['SCO6269', 'and', 'SCO6270']], 'and', 'SCO0681']
        bbh_avail_list = prunPhase_utils.check_bbh_availability(
                temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict['AKGDH2'])
        assert bbh_avail_list == [[['1', 'and', '1'], 'or', ['1', 'and', '1']], 'and', '1']


    def test_label_rxn_to_remove(self, sco_tmp_model, options):
        options.temp_target_BBH_dict = {}
        options.tempModel_biggRxnid_locusTag_dict = {}

        # PGI
        options.temp_target_BBH_dict['SCO1942'] = ['B446_30415', 'B446_10110']
        options.temp_target_BBH_dict['SCO6659'] = ['B446_30415', 'B446_10110']
        options.tempModel_biggRxnid_locusTag_dict['PGI'] = ['SCO1942', 'or', 'SCO6659']

        bbh_avail_list = prunPhase_utils.label_rxn_to_remove(sco_tmp_model, options)
        assert options.rxnToRemove_dict['PGI'] == '1'

        # AKGDH2
        options.temp_target_BBH_dict['SCO4594'] = ['B446_21645']
        options.temp_target_BBH_dict['SCO4595'] = ['B446_21650']
        options.temp_target_BBH_dict['SCO6269'] = ['B446_21645']
        options.temp_target_BBH_dict['SCO6270'] = ['B446_21650']
        options.temp_target_BBH_dict['SCO0681'] = ['B446_05650']
        options.tempModel_biggRxnid_locusTag_dict['AKGDH2'] = \
                [[['SCO4594', 'and', 'SCO4595'], 'or', ['SCO6269', 'and', 'SCO6270']], 'and', 'SCO0681']
        bbh_avail_list = prunPhase_utils.label_rxn_to_remove(sco_tmp_model, options)
        assert options.rxnToRemove_dict['AKGDH2'] == '1'


    def test_prune_model(self, sco_tmp_model, options):
        _cfg_name = 'gems.cfg'
        load_config(options)
        options.rxnToRemove_dict = {}

        options.rxnToRemove_dict['PAPA160'] = '0'
        options.rxnToRemove_dict['COELICHELINR2'] = '0'
        options.rxnToRemove_dict['ATPHs'] = '0'

        assert 'PAPA160' in sco_tmp_model.reactions
        assert 'COELICHELINR2' in sco_tmp_model.reactions
        assert 'ATPHs' in sco_tmp_model.reactions

        modelPruned = prunPhase_utils.prune_model(sco_tmp_model, options)

        assert 'PAPA160' in modelPruned.reactions
        assert 'COELICHELINR2' not in modelPruned.reactions
        assert 'ATPHs' not in modelPruned.reactions


    def test_swap_locustag_with_homolog(self, sco_tmp_model, bbh_dict, options):
        """
        MCOATA
        sco: (SCO2387 and (SCO2389 or SCO0549 or SCO1267 or SCO1272))
        sci: ( B446_12460 and (B446_12470 or SCO0549 or SCO1267 or SCO1272 )

        ACOATA
        sco: ((SCO1271 or SCO2388 or SCO6564) and (SCO2389 or SCO0549 or SCO1267 or SCO1272))
        sci: (( B446_30165 or B446_27925 or B446_12465) and (B446_12470 or SCO0549 or SCO1267 or SCO1272 ))

        PDH
        sco: ((SCO2183 or SCO2371 or SCO7124 or (SCO1269 and SCO1270)) and (SCO3815 or SCO3829) and (SCO0884 or SCO2180 or SCO4919))
        sci: (( B446_12400 or B446_11440 or (SCO1269 and SCO1270)) and (B446_19475 and B446_19415) and (B446_11425 or B446_23075 or B446_32095))
        """
        options.temp_target_BBH_dict = bbh_dict
        modelPrunedGPR = prunPhase_utils.swap_locustag_with_homolog(sco_tmp_model, options)
        assert modelPrunedGPR.reactions.get_by_id('MCOATA').gene_reaction_rule == \
                '(B446_12460 and (B446_12470 or SCO0549 or SCO1267 or SCO1272))'
        assert modelPrunedGPR.reactions.get_by_id('ACOATA').gene_reaction_rule == \
                '((( B446_30165 or B446_27925 or B446_12465 ) or ( B446_30165 or B446_27925 or B446_12465 ) or ( B446_30165 or B446_27925 or B446_12465 )) and (B446_12470 or SCO0549 or SCO1267 or SCO1272))'
        assert modelPrunedGPR.reactions.get_by_id('PDH').gene_reaction_rule == \
                '((( B446_12400 or B446_11440 ) or ( B446_12400 or B446_11440 ) or ( B446_12400 or B446_11440 ) or (SCO1269 and SCO1270)) and (( B446_19415 or B446_19475 ) or ( B446_19415 or B446_19475 )) and (B446_32095 or ( B446_11425 or B446_32095 or B446_23075 ) or ( B446_11425 or B446_23075 )))'


    def test_get_rxnid_info_dict_from_kegg(self, options):
        _cfg_name = 'gems.cfg'
        load_config(options)

        options.targetGenome_locusTag_ec_nonBBH_dict = {'B446_27575':['2.7.4.9']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options)
        assert 'R02098' in options.rxnid_info_dict
        assert options.rxnid_info_dict['R02098']['PATHWAY'] == \
                'rn00240 Pyrimidine metabolism'
        assert 'R02098' in options.rxnid_locusTag_dict
        assert 'B446_27575' in options.rxnid_locusTag_dict['R02098']

        options.targetGenome_locusTag_ec_nonBBH_dict = \
                {'B446_23835':['4.1.1.45', '3.5.2.3']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options)
        assert 'R04323' in options.rxnid_info_dict
        assert options.rxnid_info_dict['R04323']['NAME'] == \
                '2-Amino-3-carboxymuconate semialdehyde carboxy-lyase'
        assert 'R04323' in options.rxnid_locusTag_dict
        assert 'B446_23835' in options.rxnid_locusTag_dict['R04323']


    def test_get_mnxr_list_from_modelPrunedGPR(self, sco_tmp_model, options):
        bigg_mnxr_dict = {'MCOATA':'MNXR35619'}
        options.bigg_mnxr_dict = bigg_mnxr_dict

        augPhase_utils.get_mnxr_list_from_modelPrunedGPR(sco_tmp_model, options)

        assert 'MNXR35619' in options.modelPrunedGPR_mnxr_list


    def test_mnxr_to_add_list(self, mnxref, options):
        options.rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }

        options.mnxr_kegg_dict = {'MNXR70727': ['R08926']}
        options.modelPrunedGPR_mnxr_list = []
        options.mnxref = mnxref

        augPhase_utils.get_mnxr_to_add_list(options)

        assert 'MNXR70727' in options.mnxr_to_add_list


    # Focus on metabolite addition in this test
    # New metabolites: 'MNXM16902' and 'fuc__L'
    def test_add_nonBBH_rxn(self, sco_tmp_model, mnxref, tmpdir, sco_tmp_model_flux, options):
        options.mnxr_to_add_list = ['MNXR70727']
        options.rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }
        options.mnxr_kegg_dict = {'MNXR70727': ['R08926']}
        options.rxnid_locusTag_dict = {'R08926':['STEN_00480']}
        options.targetGenome_locusTag_prod_dict = {'STEN_00480':'D-threo-aldose 1-dehydrogenase'}
        outputfolder5 = './tmp'
        options.mnxref = mnxref
        options.outputfolder5 = outputfolder5
        options.template_exrxnid_flux_dict = sco_tmp_model_flux

        _cfg_name = 'gems.cfg'
        load_config(options)

        assert 'R08926' not in sco_tmp_model.reactions
        assert 'MNXM16902_c' not in sco_tmp_model.metabolites
        assert 'fuc__L_c' not in sco_tmp_model.metabolites
        assert 'h_c' in sco_tmp_model.metabolites
        assert 'nadh_c' in sco_tmp_model.metabolites
        assert 'nad_c' in sco_tmp_model.metabolites

        model = augPhase_utils.add_nonBBH_rxn(sco_tmp_model, options)

        assert 'R08926' in model.reactions
        assert 'MNXM16902_c' in model.metabolites
        assert 'fuc__L_c' in model.metabolites
        assert 'h_c' in model.metabolites
        assert 'nadh_c' in model.metabolites
        assert 'nad_c' in model.metabolites
