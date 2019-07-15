
import warnings
from os.path import join
from gmsm.config import load_config
from gmsm.io.input_file_manager import get_eficaz_file, get_locustag_comp_dict
from gmsm.primary_model import prunPhase_utils, augPhase_utils

warnings.filterwarnings("ignore")

class TestPrimary_model:
    """Test functions in gmsm.primary_model"""

    def test_get_rxn_fate(self):
        temp_target_BBH_dict = {}
        tempModel_biggRxnid_locusTag_dict = {}

        # PGI
        temp_target_BBH_dict['SCO1942'] = ['B446_30415', 'B446_10110']
        temp_target_BBH_dict['SCO6659'] = ['B446_30415', 'B446_10110']
        tempModel_biggRxnid_locusTag_dict['PGI'] = ['SCO1942', 'or', 'SCO6659']

        rxn_fate = prunPhase_utils.get_rxn_fate(
                tempModel_biggRxnid_locusTag_dict['PGI'], temp_target_BBH_dict)
        assert rxn_fate == '1'

        # AKGDH2
        temp_target_BBH_dict['SCO4594'] = ['B446_21645']
        temp_target_BBH_dict['SCO4595'] = ['B446_21650']
        temp_target_BBH_dict['SCO6269'] = ['B446_21645']
        temp_target_BBH_dict['SCO6270'] = ['B446_21650']
        temp_target_BBH_dict['SCO0681'] = ['B446_05650']
        tempModel_biggRxnid_locusTag_dict['AKGDH2'] = \
                [[['SCO4594', 'and', 'SCO4595'], 'or', ['SCO6269', 'and', 'SCO6270']], 'and', 'SCO0681']
        rxn_fate = prunPhase_utils.get_rxn_fate(
                tempModel_biggRxnid_locusTag_dict['AKGDH2'], temp_target_BBH_dict)
        assert rxn_fate == '1'

        # LDH_D
        temp_target_BBH_dict['SCO2118'] = ['B446_25870']
        tempModel_biggRxnid_locusTag_dict['LDH_D'] = ['SCO2118', 'or', 'SCO3594']

        rxn_fate = prunPhase_utils.get_rxn_fate(
                tempModel_biggRxnid_locusTag_dict['LDH_D'], temp_target_BBH_dict)
        assert rxn_fate == '1'

        # TGBPA
        temp_target_BBH_dict['SCO5852'] = ['B446_05445']
        tempModel_biggRxnid_locusTag_dict['TGBPA'] = ['SCO5848', 'and', 'SCO5852']

        rxn_fate = prunPhase_utils.get_rxn_fate(
                tempModel_biggRxnid_locusTag_dict['TGBPA'], temp_target_BBH_dict)
        assert rxn_fate == '0'

        # NOTE: This issue has not been resolved. OR is returned regardless of the Boolean.
        #If successful, the result should be 'rxn_fate == 0'.
        # NOTE: locustags should be given.
#        bbh_avail_list = [['1', 'and', '0'], 'or', ['0', 'and', '1'], 'and', ['0', 'or', '1']]
#        rxn_fate = prunPhase_utils.get_rxn_fate(bbh_avail_list)
#        assert rxn_fate != '0'


    def test_label_rxn_to_remove(self, sco_tmp_model, options):
        options.temp_target_BBH_dict = {}
        options.tempModel_biggRxnid_locusTag_dict = {}

        # PGI
        options.temp_target_BBH_dict['SCO1942'] = ['B446_30415', 'B446_10110']
        options.temp_target_BBH_dict['SCO6659'] = ['B446_30415', 'B446_10110']
        options.tempModel_biggRxnid_locusTag_dict['PGI'] = ['SCO1942', 'or', 'SCO6659']

        prunPhase_utils.label_rxn_to_remove(sco_tmp_model, options, options, options)
        assert options.rxnToRemove_dict['PGI'] == '1'

        # AKGDH2
        options.temp_target_BBH_dict['SCO4594'] = ['B446_21645']
        options.temp_target_BBH_dict['SCO4595'] = ['B446_21650']
        options.temp_target_BBH_dict['SCO6269'] = ['B446_21645']
        options.temp_target_BBH_dict['SCO6270'] = ['B446_21650']
        options.temp_target_BBH_dict['SCO0681'] = ['B446_05650']
        options.tempModel_biggRxnid_locusTag_dict['AKGDH2'] = \
                [[['SCO4594', 'and', 'SCO4595'], 'or', ['SCO6269', 'and', 'SCO6270']], 'and', 'SCO0681']
        prunPhase_utils.label_rxn_to_remove(sco_tmp_model, options, options, options)
        assert options.rxnToRemove_dict['AKGDH2'] == '1'

        # TGBPA
        options.temp_target_BBH_dict['SCO5852'] = ['B446_05445']
        options.tempModel_biggRxnid_locusTag_dict['TGBPA'] = ['SCO5848', 'and', 'SCO5852']
        prunPhase_utils.label_rxn_to_remove(sco_tmp_model, options, options, options)
        assert options.rxnToRemove_dict['TGBPA'] == '0'


    def test_prune_model(self, sco_tmp_model, options):
        _cfg_name = 'gmsm.cfg'
        load_config(options)
        options.rxnToRemove_dict = {}

        options.rxnToRemove_dict['PAPA160'] = '0'
        options.rxnToRemove_dict['COELICHELINR2'] = '0'
        options.rxnToRemove_dict['ATPHs'] = '0'

        assert 'PAPA160' in sco_tmp_model.reactions
        assert 'COELICHELINR2' in sco_tmp_model.reactions
        assert 'ATPHs' in sco_tmp_model.reactions

        modelPruned = prunPhase_utils.prune_model(sco_tmp_model, options, options)

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
                '((B446_12465 or B446_27925 or B446_30165) and (B446_12470 or SCO0549 or SCO1267 or SCO1272))'
        assert modelPrunedGPR.reactions.get_by_id('PDH').gene_reaction_rule == \
                '((( B446_11440 or B446_12400 ) or ( B446_11440 or B446_12400 ) or ( B446_11440 or B446_12400 ) or (SCO1269 and SCO1270)) and (B446_19415 or B446_19475) and (B446_11425 or B446_23075 or B446_32095))'


    def test_get_rxnid_from_ECNumber(self, options):
        _cfg_name = 'gmsm.cfg'
        load_config(options)

        rxnid_list = []
        enzymeEC = '2.7.1.22'
        rxnid_list = augPhase_utils.get_rxnid_from_ECNumber(rxnid_list, enzymeEC, options)
        assert 'R02324' in rxnid_list

        rxnid_list = []
        # This EC number is the one deleted from KEGG in 2016.
        # Error caused by this EC number was previously reported
        #for antiSMASH 3.0 (Aug 16, 2016)
        enzymeEC = '2.7.1.69'
        rxnid_list = augPhase_utils.get_rxnid_from_ECNumber(rxnid_list, enzymeEC, options)
        assert not rxnid_list


    def test_get_rxnid_info_dict_from_kegg(self, options):
        _cfg_name = 'gmsm.cfg'
        load_config(options)

        options.targetGenome_locusTag_ec_nonBBH_dict = {'B446_27575':['2.7.4.9']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options, options, options)
        assert 'R02098' in options.rxnid_info_dict
        assert options.rxnid_info_dict['R02098']['PATHWAY'] == \
                'rn00240 Pyrimidine metabolism'
        assert 'R02098' in options.rxnid_locusTag_dict
        assert 'B446_27575' in options.rxnid_locusTag_dict['R02098']

        options.targetGenome_locusTag_ec_nonBBH_dict = \
                {'B446_23835':['4.1.1.45', '3.5.2.3']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options, options, options)
        assert 'R04323' in options.rxnid_info_dict
        assert options.rxnid_info_dict['R04323']['NAME'] == \
                '2-Amino-3-carboxymuconate semialdehyde carboxy-lyase'
        assert 'R04323' in options.rxnid_locusTag_dict
        assert 'B446_23835' in options.rxnid_locusTag_dict['R04323']


    def test_get_mnxr_list_from_modelPrunedGPR(self, sco_tmp_model, options):
        bigg_mnxr_dict = {'MCOATA':'MNXR101421'}
        options.bigg_mnxr_dict = bigg_mnxr_dict

        augPhase_utils.get_mnxr_list_from_modelPrunedGPR(sco_tmp_model, options, options)

        assert 'MNXR101421' in options.modelPrunedGPR_mnxr_list


    def test_mnxr_to_add_list(self, mnxref, options):
        options.rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }

        options.mnxr_kegg_dict = {'MNXR112417': ['R08926']}
        options.modelPrunedGPR_mnxr_list = []
        options.mnxref = mnxref

        augPhase_utils.get_mnxr_to_add_list(options, options)

        assert 'MNXR112417' in options.mnxr_to_add_list


    # Focus on metabolite addition in this test
    # New metabolites: 'MNXM16902' and 'fuc__L'
    def test_add_nonBBH_rxn(self, sco_tmp_model, mnxref, tmpdir, sco_tmp_model_flux, options):
        options.mnxr_to_add_list = ['MNXR112417']
        options.rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }
        options.mnxr_kegg_dict = {'MNXR112417': ['R08926']}
        options.rxnid_locusTag_dict = {'R08926':['STEN_00480']}
        options.targetGenome_locusTag_prod_dict = {'STEN_00480':'D-threo-aldose 1-dehydrogenase'}
        outputfolder5 = './tmp'
        options.mnxref = mnxref
        options.outputfolder5 = outputfolder5
        options.template_exrxnid_flux_dict = sco_tmp_model_flux

        _cfg_name = 'gmsm.cfg'
        load_config(options)

        assert 'R08926' not in sco_tmp_model.reactions
        assert 'MNXM16902_c' not in sco_tmp_model.metabolites
        assert 'MNXM659_c' not in sco_tmp_model.metabolites
        assert 'h_c' in sco_tmp_model.metabolites
        assert 'nadh_c' in sco_tmp_model.metabolites
        assert 'nad_c' in sco_tmp_model.metabolites

        model = augPhase_utils.add_nonBBH_rxn(sco_tmp_model, options, options, options)

        assert 'R08926' in model.reactions
        assert 'MNXM16902_c' in model.metabolites
        assert 'MNXM659_c' in model.metabolites
        assert 'h_c' in model.metabolites
        assert 'nadh_c' in model.metabolites
        assert 'nad_c' in model.metabolites


    def test_get_rxn_newComp_list_from_model(self, sci_primary_model, options):

        options.locustag_comp_dict = {}
        options.locustag_comp_dict['B446_25420'] = ['c']
        rxn_newComp_list = \
                augPhase_utils.get_rxn_newComp_list_from_model(sci_primary_model, options)

        assert len(rxn_newComp_list) == 2
        assert 'ACKr' in rxn_newComp_list


    def test_create_rxn_newComp1(self, sci_primary_model, options):

        options.locustag_comp_dict = {}
        options.locustag_comp_dict['B446_25420'] = ['c']
        rxn_newComp_list = \
                augPhase_utils.get_rxn_newComp_list_from_model(sci_primary_model, options)

        options.outputfolder5 = './tmp'
        model, added_rxn_newComp_list = augPhase_utils.create_rxn_newComp(
                                           rxn_newComp_list, sci_primary_model, options, options)
        assert len(model.reactions) == int(1805)


    def test_create_rxn_newComp2(self, sci_primary_model, options):

        options.locustag_comp_dict = {}
        options.locustag_comp_dict['B446_25420'] = ['p', 'm']
        rxn_newComp_list = \
                augPhase_utils.get_rxn_newComp_list_from_model(sci_primary_model, options)

        assert 'ACKrp' not in sci_primary_model.reactions
        assert 'PPAKrp' not in sci_primary_model.reactions
        assert 'ACKrm' not in sci_primary_model.reactions
        assert 'PPAKrm' not in sci_primary_model.reactions

        assert 'ac_p' not in sci_primary_model.metabolites
        assert 'ac_m' not in sci_primary_model.metabolites
        assert 'ppap_p' not in sci_primary_model.metabolites
        assert 'ppap_m' not in sci_primary_model.metabolites

        assert len(sci_primary_model.reactions) == int(1805)
        assert len(sci_primary_model.metabolites) == int(1582)

        options.outputfolder5 = './tmp'
        model, added_rxn_newComp_list = augPhase_utils.create_rxn_newComp(
                                           rxn_newComp_list, sci_primary_model, options, options)

        assert 'ACKrp' in model.reactions
        assert 'PPAKrp' in model.reactions
        assert 'ACKrm' in model.reactions
        assert 'PPAKrm' in model.reactions
        assert 'ACKrp' in added_rxn_newComp_list
        assert 'PPAKrp' in added_rxn_newComp_list
        assert 'ACKrm' in added_rxn_newComp_list
        assert 'PPAKrm' in added_rxn_newComp_list

        assert 'ac_p' in model.metabolites
        assert 'ac_m' in model.metabolites
        assert 'ppap_p' in model.metabolites
        assert 'ppap_m' in model.metabolites

        assert len(model.reactions) == int(1809)
        assert len(model.metabolites) == int(1594)
        assert len(added_rxn_newComp_list) == 4


    def test_remove_inactive_rxn_newComp(self, sci_primary_model, tmpdir, options):

        # Reaction 'CSND' is an inactive reaction in 'sci_primary_model'
        added_rxn_newComp_list = ['CSND']
        options.outputfolder5 = './tmp'

        assert len(sci_primary_model.reactions) == int(1805)
        assert 'CSND' in sci_primary_model.reactions

        model = augPhase_utils.remove_inactive_rxn_newComp(
                                                            added_rxn_newComp_list,
                                                            sci_primary_model,
                                                            options, options)

        assert len(model.reactions) == 1804
        assert 'CSND' not in model.reactions
        assert 'CSND' in options.inactive_rxn_newComp_list
