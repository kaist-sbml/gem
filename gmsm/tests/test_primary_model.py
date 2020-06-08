
import warnings
from os.path import join, abspath, dirname
from gmsm.config import load_config
from gmsm.io.input_file_manager import get_eficaz_file, get_locustag_comp_dict
from gmsm.primary_model import prunPhase_utils, augPhase_utils, run_primary_modeling

warnings.filterwarnings("ignore")

class TestPrimary_model:
    """Test functions in gmsm.primary_model"""

    #-------------------------------------------
    # test functions of prunPhase_utils.py
    #-------------------------------------------

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

        # 'AND' and 'OR' within a parenthesis
        bbh_avail_list = [['1', 'and', '0'], 'or', ['0', 'and', '1'], 'and', ['0', 'or', '1']]
        rxn_fate = prunPhase_utils.get_rxn_fate(bbh_avail_list, temp_target_BBH_dict)
        assert rxn_fate == '0'


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
                'B446_12460 and (B446_12470 or SCO0549 or SCO1267 or SCO1272)'
        assert modelPrunedGPR.reactions.get_by_id('ACOATA').gene_reaction_rule == \
                '(( B446_30165 or B446_27925 or B446_12465 ) or ( B446_30165 or B446_27925 or B446_12465 ) or ( B446_30165 or B446_27925 or B446_12465 )) and (B446_12470 or SCO0549 or SCO1267 or SCO1272)'
        assert modelPrunedGPR.reactions.get_by_id('PDH').gene_reaction_rule == \
                '(( B446_12400 or B446_11440 ) or ( B446_12400 or B446_11440 ) or ( B446_12400 or B446_11440 ) or (SCO1269 and SCO1270)) and (( B446_19415 or B446_19475 ) or ( B446_19415 or B446_19475 )) and (B446_32095 or ( B446_11425 or B446_32095 or B446_23075 ) or ( B446_11425 or B446_23075 ))'

    #-------------------------------------------
    # test functions of augPhase_utils.py
    #-------------------------------------------

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


    def test_get_rxnInfo_from_rxnid(self, options):
        _cfg_name = 'gmsm.cfg'
        load_config(options)

        rxnid1 = 'R00001' # polyphosphate polyphosphohydrolase
        rxnid2 = 'R00771' # D-glucose-6-phosphate aldose-ketose-isomerase
        rxnid_info1 = augPhase_utils.get_rxnInfo_from_rxnid(rxnid1, options)
        rxnid_info2 = augPhase_utils.get_rxnInfo_from_rxnid(rxnid2, options)
        assert rxnid_info1 == None
        assert rxnid_info2 == {'NAME': 'D-glucose-6-phosphate aldose-ketose-isomerase', 'ENZYME': '5.3.1.9', 'PATHWAY': 'rn00500 Starch and sucrose metabolism', 'DEFINITION': 'D-Glucose 6-phosphate <=> D-Fructose 6-phosphate', 'EQUATION': 'C00092 <=> C00085'}


    def test_load_cache(self):
        cache_ec_rxn_dict = {}
        cache_rxnid_info_dict = {}
        cache_dumped_ec_list = []
        cache_dumped_rxnid_list = []
        nonexistence_dict = {}
        nonexistence_list = []
        nonpickle = ""

        data_model_dir = join(dirname(abspath(__file__)), 'data')
        cache_ec_rxn_dict_dir = join(data_model_dir, 'cache_ec_rxn_dict.p')
        cache_rxnid_info_dict_dir = join(data_model_dir, 'cache_rxnid_info_dict.p')
        cache_dumped_ec_list_dir = join(data_model_dir, 'cache_dumped_ec_list.p')
        cache_dumped_rxnid_list_dir = join(data_model_dir, 'cache_dumped_rxnid_list.p')
        nonexistence_dict_dir = join(data_model_dir, 'nonexistence_dict.p')
        nonexistence_list_dir = join(data_model_dir, 'nonexistence_list.p')
        nonpickle_dir = join(data_model_dir, 'nonpickle.png')

        cache_ec_rxn_dict = augPhase_utils.load_cache(
                cache_ec_rxn_dict_dir, cache_ec_rxn_dict)
        cache_rxnid_info_dict = augPhase_utils.load_cache(
                cache_rxnid_info_dict_dir, cache_rxnid_info_dict)
        cache_dumped_ec_list = augPhase_utils.load_cache(
                cache_dumped_ec_list_dir, cache_dumped_ec_list)
        cache_dumped_rxnid_list = augPhase_utils.load_cache(
                cache_dumped_rxnid_list_dir, cache_dumped_rxnid_list)
        nonexistence_dict = augPhase_utils.load_cache(
                nonexistence_dict_dir, nonexistence_dict)
        nonexistence_list = augPhase_utils.load_cache(
                nonexistence_list_dir, nonexistence_list)
        nonpickle = augPhase_utils.load_cache(
                nonpickle_dir, nonpickle)

        assert cache_ec_rxn_dict != {}
        assert cache_rxnid_info_dict != {}
        assert cache_dumped_ec_list != []
        assert cache_dumped_rxnid_list != []
        assert nonexistence_dict == {}
        assert nonexistence_list == []
        assert nonpickle == ""


    def test_get_targetGenome_locusTag_ec_nonBBH_dict(self, options):
        options.nonBBH_list = ['B446_27575', 'B446_23835']
        options.targetGenome_locusTag_ec_dict = \
                              {'B446_23835':['4.1.1.45', '3.5.2.3']}

        augPhase_utils.get_targetGenome_locusTag_ec_nonBBH_dict(options, options, options)

        assert options.targetGenome_locusTag_ec_nonBBH_dict == {'B446_23835':['4.1.1.45', '3.5.2.3']}


    def test_edit_mnxr_kegg_dict(self, mnxr_kegg_dict, options):
        options.mnxr_kegg_dict = mnxr_kegg_dict
        mnxr_kegg_len = len(options.mnxr_kegg_dict)
        keggid = 'R08385'
        keggid2 = 'R04783'

        augPhase_utils.edit_mnxr_kegg_dict(keggid, options)
        augPhase_utils.edit_mnxr_kegg_dict(keggid2, options)

        assert 'R08385' not in options.mnxr_kegg_dict.values()
        assert 'R04783' not in options.mnxr_kegg_dict.values()
        assert len(options.mnxr_kegg_dict) == mnxr_kegg_len - 1


    def test_get_rxnid_locusTag_dict(self):
        rxnid_locusTag_dict = {}
        rxnid = 'R04558'
        locusTag = 'SCO2048'
        locusTag2 = 'SCO2051'

        rxnid_locustag_dict = \
        augPhase_utils.get_rxnid_locusTag_dict(rxnid_locusTag_dict, rxnid, locusTag)

        rxnid_locustag_dict = \
        augPhase_utils.get_rxnid_locusTag_dict(rxnid_locusTag_dict, rxnid, locusTag2)

        assert rxnid_locustag_dict == {'R04558' : ['SCO2048', 'SCO2051']}


    def test_get_rxnid_info_dict_from_kegg(self, mnxr_kegg_dict, options):
        _cfg_name = 'gmsm.cfg'
        load_config(options)

        options.mnxr_kegg_dict = mnxr_kegg_dict

        # This dictionary is for a case : EC_number info from a cache file
        options.targetGenome_locusTag_ec_nonBBH_dict = {'B446_27575':['2.7.4.9']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options, options, options)
        assert 'R02098' in options.rxnid_info_dict
        assert options.rxnid_info_dict['R02098']['PATHWAY'] == \
                'rn00240 Pyrimidine metabolism'
        assert 'R02098' in options.rxnid_locusTag_dict
        assert 'B446_27575' in options.rxnid_locusTag_dict['R02098']

        # This dictionary is for a case : EC_number info from a cache file
        options.targetGenome_locusTag_ec_nonBBH_dict = \
                {'B446_23835':['4.1.1.45', '3.5.2.3']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options, options, options)
        assert 'R04323' in options.rxnid_info_dict
        assert options.rxnid_info_dict['R04323']['NAME'] == \
                '2-Amino-3-carboxymuconate semialdehyde carboxy-lyase'
        assert 'R04323' in options.rxnid_locusTag_dict
        assert 'B446_23835' in options.rxnid_locusTag_dict['R04323']

        # This dictionary is for a case : EC number not available at KEGG
        options.targetGenome_locusTag_ec_nonBBH_dict = \
                {'B446_07840':['3.1.22.4']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options, options, options)

        # This dictionary is for a case : EC number info fetched from KEGG
        # Depending on the accumulation of the cache file,
        # it may be the case of EC_number info from a cache file
        options.targetGenome_locusTag_ec_nonBBH_dict = \
                {'B446_00385':['3.5.1.90']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options, options, options)
        assert 'R05226' in options.rxnid_info_dict
        assert options.rxnid_info_dict['R05226']['NAME'] == \
                'adenosylcobinamide amidohydrolase'
        assert 'R05226' in options.rxnid_locusTag_dict
        assert 'B446_00385' in options.rxnid_locusTag_dict['R05226']

        # This dictionary is for a case :
        # rxnid not in cache_rxnid_info_dict but in cache_dumped_rxnid_list
        options.targetGenome_locusTag_ec_nonBBH_dict = \
                {'B446_06350':['2.7.1.11']}
        augPhase_utils.get_rxnid_info_dict_from_kegg(options, options, options)
        assert 'R01843' not in options.rxnid_info_dict
        assert 'R01843' in options.rxnid_locusTag_dict
        assert 'B446_06350' in options.rxnid_locusTag_dict['R01843']
        assert 'MNXR102510' not in options.mnxr_kegg_dict.keys()
        assert 'R01843' not in options.mnxr_kegg_dict.values()


    def test_get_mnxr_list_from_modelPrunedGPR(self, sco_tmp_model, options):
        bigg_mnxr_dict = {'MCOATA':'MNXR101421'}
        options.bigg_mnxr_dict = bigg_mnxr_dict

        augPhase_utils.get_mnxr_list_from_modelPrunedGPR(sco_tmp_model, options, options)

        assert 'MNXR101421' in options.modelPrunedGPR_mnxr_list


    def test_get_mnxr_to_add_list(self, mnxref, options):
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


    def test_add_nonBBH_rxn(self, sco_tmp_model, mnxref, tmpdir, sco_tmp_model_flux, options):
        # Focus on metabolite addition in this test
        # New metabolites: 'MNXM16902' and 'fuc__L'
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
        options.rxnid_locusTag_dict = {'R08926': ['STEN_00480']}
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
        assert 'fuc__L_c' in model.metabolites
        assert 'h_c' in model.metabolites
        assert 'nadh_c' in model.metabolites
        assert 'nad_c' in model.metabolites

        # Focus on passing GPR association with subunit
        options.mnxr_to_add_list = ['MNXR102634']
        options.rxnid_info_dict = {
            'R03660': {
                'ENZYME': '6.1.1.20',
                'DEFINITION': 'ATP + L-Phenylalanine + tRNA(Phe) <=> AMP + Diphosphate + L-Phenylalanyl-tRNA(Phe)',
                'EQUATION': 'C00002 + C00079 + C01648 <=> C00020 + C00013 + C03511',
                'NAME': 'L-Phenylalanine:tRNA(Ala) ligase (AMP-forming)',
                'PATHWAY': 'rn00970 Aminoacyl-tRNA biosynthesis'}}
        options.mnxr_kegg_dict = {'MNXR102634': ['R03660']}
        options.rxnid_locusTag_dict = {'R03660': ['SCO1594','SCO1595']}
        options.targetGenome_locusTag_prod_dict = {'SCO1594': 'phenylalanyl-tRNA synthetase subunit beta', 'SCO1595': 'phenylalanyl-tRNA synthetase subunit alpha'}

        assert 'R03660' not in sco_tmp_model.reactions

        model = augPhase_utils.add_nonBBH_rxn(sco_tmp_model, options, options, options)
        assert 'R03660' in model.reactions


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
        
        #According to cobrapy updates, the number of reactions are changed 1805 to 2009
        assert len(model.reactions) == int(2009)


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
        
        #According to cobrapy updates, the number of reactions are changed 1805 to 2009
        #According to cobrapy updates, the number of metabolites are changed 1582 to 1786
        assert len(sci_primary_model.reactions) == int(2009)
        assert len(sci_primary_model.metabolites) == int(1786)

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

        assert len(model.reactions) == int(2013)
        assert len(model.metabolites) == int(1798)
        assert len(added_rxn_newComp_list) == 4


    def test_remove_inactive_rxn_newComp(self, sci_primary_model, tmpdir, options):

        # Reaction 'CSND' is an inactive reaction in 'sci_primary_model'
        added_rxn_newComp_list = ['CSND']
        options.outputfolder5 = './tmp'

        #According to cobrapy updates, the number of reactions are changed 1805 to 2009
        assert len(sci_primary_model.reactions) == int(2009)
        assert 'CSND' in sci_primary_model.reactions

        model = augPhase_utils.remove_inactive_rxn_newComp(
                                                            added_rxn_newComp_list,
                                                            sci_primary_model,
                                                            options, options)

        assert len(model.reactions) == int(2008)
        assert 'CSND' not in model.reactions
        assert 'CSND' in options.inactive_rxn_newComp_list

    #-------------------------------------------
    # test functions of run_primary_modeling.py
    #-------------------------------------------

    # This function is just for checking the lines in the python file.
    # Above functions are already testing the validity of each of functions.
    def test_run_primary_modeling(self, sco_tmp_model, options):

        options.tempModel_biggRxnid_locusTag_dict = {}
        options.temp_target_BBH_dict = {}
        options.nonBBH_list = []
        options.bigg_mnxr_dict = {}
        options.comp = ['c']
        options.locustag_comp_dict = {}
        options.outputfolder5 = './tmp'
        model = run_primary_modeling.run_prunPhase(sco_tmp_model, options, options, options, options)
        run_primary_modeling.run_augPhase(model, options, options, options, options, options)
