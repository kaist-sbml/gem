
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from os.path import join
from gems.config import load_config
from gems.primary_model import prunPhase_utils, augPhase_utils

class TestPrimary_model:
    """Test functions in gems.primary_model"""


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


    def test_mnxr_to_add_list(self, mnxref, options):
        rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }

        modelPrunedGPR_mnxr_list = []

        options.rxnid_info_dict = rxnid_info_dict
        options.mnxref = mnxref
        options.modelPrunedGPR_mnxr_list = modelPrunedGPR_mnxr_list

        augPhase_utils.get_mnxr_to_add_list(options)

        assert 'MNXR70727' in options.mnxr_to_add_list


    # Focus on metabolite addition in this test
    # New metabolites: 'MNXM16902' and 'fuc_DASH_L'
    def test_add_nonBBH_rxn(self, sco_tmp_model, mnxref, tmpdir, sco_tmp_model_flux, options):
        mnxr_to_add_list = ['MNXR70727']
        rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }
        kegg_mnxr_dict = {'R08926':'MNXR70727'}
        rxnid_locusTag_dict = {'R08926':['STEN_00480']}
        targetGenome_locusTag_prod_dict = {'STEN_00480':'D-threo-aldose 1-dehydrogenase'}
        outputfolder5 = './tmp'

        options.mnxref = mnxref
        options.mnxr_to_add_list = mnxr_to_add_list
        options.rxnid_info_dict = rxnid_info_dict
        options.kegg_mnxr_dict = kegg_mnxr_dict
        options.rxnid_locusTag_dict = rxnid_locusTag_dict
        options.targetGenome_locusTag_prod_dict = targetGenome_locusTag_prod_dict
        options.outputfolder5 = outputfolder5
        options.template_exrxnid_flux_dict = sco_tmp_model_flux

        _cfg_name = 'gems.cfg'
        load_config(options)

        assert 'R08926' not in sco_tmp_model.reactions
        assert 'MNXM16902_c' not in sco_tmp_model.metabolites
        assert 'fuc_DASH_L_c' not in sco_tmp_model.metabolites
        assert 'h_c' in sco_tmp_model.metabolites
        assert 'nadh_c' in sco_tmp_model.metabolites
        assert 'nad_c' in sco_tmp_model.metabolites

        model = augPhase_utils.add_nonBBH_rxn(sco_tmp_model, options)

        assert 'R08926' in model.reactions
        assert 'MNXM16902_c' in model.metabolites
        assert 'fuc_DASH_L_c' in model.metabolites
        assert 'h_c' in model.metabolites
        assert 'nadh_c' in model.metabolites
        assert 'nad_c' in model.metabolites
