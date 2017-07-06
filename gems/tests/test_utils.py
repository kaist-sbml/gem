
import warnings
from cobra import Reaction, Metabolite
from gems.utils import check_duplicate_rxn, compare_rxns

warnings.filterwarnings("ignore")

class TestUtils:
    """Test functions in gems.utils"""

    def test_check_duplicate_rxn1(self, sco_tmp_model, options):

        ala__L = Metabolite('ala_DASH_L_c')
        h2o = Metabolite('h2o_c')
        nad = Metabolite('nad_c')
        h = Metabolite('h_c')
        nadh = Metabolite('nadh_c')
        nh4 = Metabolite('nh4_c')
        pyr = Metabolite('pyr_c')

        rxn = Reaction('test_rxn')
        rxn.add_metabolites({
            ala__L: -1.0,
            h2o: -1.0,
            nad: -1.0,
            h: 1.0,
            nadh: 1.0,
            nh4: 1.0,
            pyr: 1.0
            })

        res = check_duplicate_rxn(sco_tmp_model, rxn)
        assert res == 'duplicate'


    def test_check_duplicate_rxn2(self, sco_tmp_model, options):

        ala__L = Metabolite('ala_DASH_L_c')
        h2o = Metabolite('h2o_c')
        nad = Metabolite('nad_c')
        nadh = Metabolite('nadh_c')
        nh4 = Metabolite('nh4_c')
        pyr = Metabolite('pyr_c')

        rxn = Reaction('test_rxn')
        rxn.add_metabolites({
            ala__L: -1.0,
            h2o: -1.0,
            nad: -1.0,
            nadh: 1.0,
            nh4: 1.0,
            pyr: 1.0
            })

        res = check_duplicate_rxn(sco_tmp_model, rxn)
        assert res == 'unique'


    def test_compare_rxns1(self, sco_tmp_model, options):

        #'ALAD_L': 'ala_DASH_L_c + h2o_c + nad_c --> h_c + nadh_c + nh4_c + pyr_c'
        rxn1 = sco_tmp_model.reactions[0]

        #'ALAR' = 'ala_DASH_L_c <=> ala_DASH_D_c'
        rxn2 = sco_tmp_model.reactions[1]
        res = compare_rxns(rxn1, rxn2)

        assert res == 'different'


    def test_compare_rxns2(self, sco_tmp_model, options):

        #'ALAD_L': 'ala_DASH_L_c + h2o_c + nad_c --> h_c + nadh_c + nh4_c + pyr_c'
        rxn1 = sco_tmp_model.reactions[0]

        ala__L = Metabolite('ala_DASH_L_c')
        h2o = Metabolite('h2o_c')
        nad = Metabolite('nad_c')
        h = Metabolite('h_c')
        nadh = Metabolite('nadh_c')
        nh4 = Metabolite('nh4_c')
        pyr = Metabolite('pyr_c')

        rxn = Reaction('test_rxn')
        rxn.add_metabolites({
            ala__L: -1.0,
            h2o: -1.0,
            nad: -1.0,
            h: 1.0,
            nadh: 1.0,
            nh4: 1.0,
            pyr: 1.0
            })

        res = compare_rxns(rxn1, rxn)

        assert res == 'same'
