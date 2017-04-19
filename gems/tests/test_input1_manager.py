
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import sys
from os.path import abspath, dirname, join

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
from scripts import input1_manager

class TestInput1_manager:
    """Test functions in scripts/input1_manager.py"""

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

