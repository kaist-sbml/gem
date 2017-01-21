
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from gems.primary_model import augPhase_utils

class TestPrimary_model:
    """Test functions in gems.primary_model"""

    def test_add_nonBBH_rxn(self, model, options):
        # All available in Namespace 'options'
        rxnid_info_dict = {
            'R08926':{
                'ENZYME': '1.1.1.122',
                'DEFINITION': '6-Deoxy-L-galactose + NAD+ <=> L-Fucono-1,5-lactone + NADH + H+',
                'EQUATION': 'C01019 + C00003 <=> C18028 + C00004 + C00080',
                'NAME': 'L-fucose:NAD+ 1-oxidoreductase',
                'PATHWAY': 'rn00051 Fructose and mannose metabolism'}
                }
        # New metabolites: 'MNXM38659' and 'fuc_DASH_L'
        rxnid_mnxm_coeff_dict = {
            'R08926':{
                'h': -1.0, 'MNXM38659': -1.0, 'nadh': -1.0,
                'fuc_DASH_L': 1.0, 'nad': 1.0}}

        bigg_mnxm_compound_dict = {}
        mnxm_compoundInfo_dict = {}
        rxnid_locusTag_dict = {}
        targetGenome_locusTag_prod_dict = {}
        outputfolder5 = ''

        options.rxnid_info_dict = rxnid_info_dict
        options.rxnid_mnxm_coeff_dict = rxnid_mnxm_coeff_dict
        options.bigg_mnxm_compound_dict = bigg_mnxm_compound_dict
        options.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict
        options.rxnid_locusTag_dict = rxnid_locusTag_dict
        options.targetGenome_locusTag_prod_dict = targetGenome_locusTag_prod_dict
        options.outputfolder5 = outputfolder5

        assert 'R08926' not in model.reactions # To be added to the model
        assert 'MNXM38659_c' not in model.metabolites # To be added to the model
        assert 'fuc_DASH_L_c' not in model.metabolites # To be added to the model
        assert 'h_c' in model.metabolites
        assert 'nadh_c' in model.metabolites
        assert 'nad_c' in model.metabolites

        #augPhase_utils.add_nonBBH_rxn(model, options)

#        assert 'R08926' in model.reactions # Should be available in the model
#        assert 'MNXM38659_c' in model.metabolites # Should be available in the model
#        assert 'fuc_DASH_L_c' in model.metabolites # Should be available in the model
#        assert 'h_c' in model.metabolites
#        assert 'nadh_c' in model.metabolites
#        assert 'nad_c' in model.metabolites
