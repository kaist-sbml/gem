
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from os.path import join
from gems.config import load_config
from gems.primary_model import prunPhase_utils, augPhase_utils

class TestPrimary_model:
    """Test functions in gems.primary_model"""

    """
    MCOATA
    sco: (SCO2387 and (SCO2389 or SCO0549 or SCO1267 or SCO1272))
    sci_current: ((B446_12460 and B446_12470) or SCO0549 or SCO1267 or SCO1272)
    sci_correct: ( B446_12460 and (B446_12470 or SCO0549 or SCO1267 or SCO1272 )

    ###
    SCO2387 B446_12460  0.0 1315.000000 316 83.860000
    SCO2389 B446_12470  2e-49   382.000000  82  93.900000
    ###

    ACOATA
    sco: ((SCO1271 or SCO2388 or SCO6564) and (SCO2389 or SCO0549 or SCO1267 or SCO1272))
    sci_current: (B446_30165 or B446_27925 or B446_12465 or (B446_12470 and B446_30165 and B446_27925 and B446_12465) or SCO0549 or SCO1267 or SCO1272)
    sci_correct: (( B446_30165 or B446_27925 or B446_12465) and (B446_12470 or SCO0549 or SCO1267 or SCO1272 ))

    ###
    SCO1271 B446_12465  5e-85   668.000000  338 42.010000
    SCO1271 B446_30165  2e-54   459.000000  317 39.750000
    SCO1271 B446_27925  1e-41   372.000000  326 34.050000

    SCO2388 B446_12465  0.0 1619.000000 351 89.170000
    SCO2388 B446_30165  6e-72   579.000000  334 41.320000
    SCO2388 B446_27925  3e-63   522.000000  334 38.620000

    SCO6564 B446_30165  0.0 1357.000000 314 86.310000
    SCO6564 B446_12465  5e-76   607.000000  342 39.770000
    SCO6564 B446_27925  1e-56   475.000000  324 37.650000
    ###

    ###
    SCO2389 B446_12470  2e-49   382.000000  82  93.900000
    ###

    PDH
    sco: ((SCO2183 or SCO2371 or SCO7124 or (SCO1269 and SCO1270)) and (SCO3815 or SCO3829) and (SCO0884 or SCO2180 or SCO4919))
    sci_current: (B446_12400 or B446_11440 or B446_12400 or B446_11440 or (B446_32095 and B446_19415 and B446_19475) or B446_11425 or B446_32095 or B446_23075)
    sci_correct: (( B446_12400 or B446_11440 or (SCO1269 and SCO1270)) and (B446_19475 and B446_19415) and (B446_11425 or B446_23075 or B446_32095))


    ###
    SCO2183 B446_11440  0.0 4392.000000 899 91.100000
    SCO2183 B446_12400  0.0 2716.000000 872 60.670000

    SCO2371 B446_12400  0.0 4424.000000 878 92.940000
    SCO2371 B446_11440  0.0 2813.000000 898 60.910000

    SCO7124 B446_11440  0.0 3403.000000 891 72.500000
    SCO7124 B446_12400  0.0 2604.000000 865 57.920000
    ###

    ###
    SCO3815 B446_19475  0.0 1795.000000 488 77.870000
    SCO3815 B446_19415  8e-81   660.000000  467 40.260000

    SCO3829 B446_19415  0.0 1557.000000 484 74.170000
    SCO3829 B446_19475  5e-78   646.000000  491 38.490000
    ###

    ###
    SCO0884 B446_32095  0.0 1923.000000 478 80.750000
    SCO0884 B446_11425  8e-33   320.000000  480 28.960000

    SCO2180 B446_11425  0.0 2208.000000 462 92.420000
    SCO2180 B446_23075  5e-47   428.000000  389 30.330000
    SCO2180 B446_32095  3e-31   308.000000  473 25.580000

    SCO4919 B446_23075  0.0 2225.000000 465 93.980000
    SCO4919 B446_11425  2e-41   387.000000  460 26.960000
    ###
    """

    def test_get_gpr_fromString_toList(self):
        MCOATA_gpr = '(SCO2387 and (SCO2389 or SCO0549 or SCO1267 or SCO1272))'
        ACOATA_gpr = '((SCO1271 or SCO2388 or SCO6564) and (SCO2389 or SCO0549 or SCO1267 or SCO1272))'
        PDH_gpr = '((SCO2183 or SCO2371 or SCO7124 or (SCO1269 and SCO1270)) and (SCO3815 or SCO3829) and (SCO0884 or SCO2180 or SCO4919))'
        PPCOAC_gpr = '((SCO2776 and SCO2777) or (SCO4380 and SCO4381) or ((SCO4921 or SCO6271) and (SCO4925 and SCO4926)))'

        gpr_list1 = prunPhase_utils.get_gpr_fromString_toList(MCOATA_gpr)
        gpr_list2 = prunPhase_utils.get_gpr_fromString_toList(ACOATA_gpr)
        gpr_list3 = prunPhase_utils.get_gpr_fromString_toList(PDH_gpr)
        gpr_list4 = prunPhase_utils.get_gpr_fromString_toList(PPCOAC_gpr)

        assert len(gpr_list1) == 2
        assert len(gpr_list2) == 2
        assert len(gpr_list3) == 4
        assert len(gpr_list4) == 4


#    def test_swap_locustag_with_homolog(self, sco_tmp_model, options):
#        gpr = prunPhase_utils.swap_locustag_with_homolog(sco_tmp_model, options)
#        assert isinstance(gpr, basestring)

    # Focus on metabolite addition in this test
    # New metabolites: 'MNXM38659' and 'fuc_DASH_L'
    def test_add_nonBBH_rxn(self, sco_tmp_model, tmpdir, sco_tmp_model_flux, options):
        rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }
        rxnid_mnxm_coeff_dict = {
            'R08926':{
                'h': -1.0, 'MNXM38659': -1.0, 'nadh': -1.0,
                'fuc_DASH_L': 1.0, 'nad': 1.0}}

        bigg_mnxm_compound_dict = {'fuc_DASH_L':'MNXM659'}
        mnxm_compoundInfo_dict = {
                'MNXM659':['L-fucose', 'C6H12O5'],
                'MNXM38659':['6-deoxy-D-glucono-1,5-lactone', 'C6H10O5']}
        rxnid_locusTag_dict = {'R08926':['STEN_00480']}
        targetGenome_locusTag_prod_dict = {'STEN_00480':'D-threo-aldose 1-dehydrogenase'}
        outputfolder5 = './tmp'

        options.rxnid_info_dict = rxnid_info_dict
        options.rxnid_mnxm_coeff_dict = rxnid_mnxm_coeff_dict
        options.bigg_mnxm_compound_dict = bigg_mnxm_compound_dict
        options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict
        options.rxnid_locusTag_dict = rxnid_locusTag_dict
        options.targetGenome_locusTag_prod_dict = targetGenome_locusTag_prod_dict
        options.outputfolder5 = outputfolder5
        options.template_exrxnid_flux_dict = sco_tmp_model_flux

        _cfg_name = 'gems.cfg'
        load_config(options)

        assert 'R08926' not in sco_tmp_model.reactions # To be added to the model
        assert 'MNXM38659_c' not in sco_tmp_model.metabolites # To be added to the model
        assert 'fuc_DASH_L_c' not in sco_tmp_model.metabolites # To be added to the model
        assert 'h_c' in sco_tmp_model.metabolites
        assert 'nadh_c' in sco_tmp_model.metabolites
        assert 'nad_c' in sco_tmp_model.metabolites

        augPhase_utils.add_nonBBH_rxn(sco_tmp_model, options)

        assert 'R08926' in sco_tmp_model.reactions # Should be available in the model
        assert 'MNXM38659_c' in sco_tmp_model.metabolites # Should be available in the model
        assert 'fuc_DASH_L_c' in sco_tmp_model.metabolites # Should be available in the model
        assert 'h_c' in sco_tmp_model.metabolites
        assert 'nadh_c' in sco_tmp_model.metabolites
        assert 'nad_c' in sco_tmp_model.metabolites
