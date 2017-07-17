
import sys
import warnings
from gems.config import load_config
from os.path import abspath, dirname, join

warnings.filterwarnings("ignore")

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
from scripts import input1_manager

class TestInput1_manager:
    """Test functions in scripts/input1_manager.py"""

    def test_get_nonstd_model(self, sco_tmp_model, options):

        # Compare model with metabolite IDs not corrected with
        #the one having metabolite IDs corrected via 'cobra.io.read_legacy_sbml'
        assert 'dad_DASH_2_c' in sco_tmp_model.metabolites
        assert 'dad__2_c' not in sco_tmp_model.metabolites
        assert 'dad_2_c' not in sco_tmp_model.metabolites

        input1_tmp_dir = join(dirname(abspath(__file__)), 'data')
        model = input1_manager.get_nonstd_model(input1_tmp_dir, options)

        assert 'dad_DASH_2_c' in model.metabolites
        assert 'dad__2_c' not in model.metabolites
        assert 'dad_2_c' not in model.metabolites


    def test_fix_nonstd_model(self, options):

        input1_tmp_dir = join(dirname(abspath(__file__)), 'data')
        model = input1_manager.get_nonstd_model(input1_tmp_dir, options)

        assert 'EX_glc_LPAREN_e_RPAREN_' in model.reactions
        assert 'EX_glc__D_e' not in model.reactions
        assert 'dad_DASH_2_c' in model.metabolites
        assert 'dad__2_c' not in model.metabolites
        assert 'dad_2_c' not in model.metabolites

        input1_tmp_dir = './tmp'
        options.acc_number = None

        _cfg_name = 'gems.cfg'
        load_config(options)

        bigg_old_new_dict = {}
        bigg_old_new_dict['EX_glc_LPAREN_e_RPAREN_'] = 'EX_glc__D_e'
        bigg_old_new_dict['dad_DASH_2_c'] = 'dad_2'

        model, model_info_dict = input1_manager.fix_nonstd_model(
                                    bigg_old_new_dict,
                                    input1_tmp_dir,
                                    model,
                                    options)

        assert 'EX_glc_LPAREN_e_RPAREN_' not in model.reactions
        assert 'EX_glc__D_e' in model.reactions
        assert 'dad_DASH_2_c' not in model.metabolites
        assert 'dad__2_c' not in model.metabolites
        assert 'dad_2_c' in model.metabolites


    def test_get_gpr_fromString_toList(self):

        gpr = '(A and B)'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr)
        assert gpr_list == ['A', 'and', 'B']

        # From: 3HAD100 in iJN746 (Pseudomonas putida KT2440)
        # Underscore should be considered in pyparsing.Word
        gpr = 'PP_1602 or PP_4174'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr)
        assert gpr_list == ['PP_1602', 'or', 'PP_4174']

        gpr = '((A and B) or (C and D))'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr)
        assert gpr_list == [['A', 'and', 'B'], 'or', ['C', 'and', 'D']]

        gpr = '((A and (B1 and B2)) or (C and D))'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr)
        assert gpr_list == [['A', 'and', ['B1', 'and', 'B2']], 'or', ['C', 'and', 'D']]

        # GPR with three nested lists fails
        gpr = '((A and (B1 or B2)) or (C and D))'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr)
        assert gpr_list == [['A', 'and', ['B1', 'or', 'B2']], 'or', ['C', 'and', 'D']]

