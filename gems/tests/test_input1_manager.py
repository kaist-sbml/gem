
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import sys
from os.path import abspath, dirname, join

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
from scripts import input1_manager

class TestInput1_manager:
    """Test functions in scripts/input1_manager.py"""

    def test_get_gpr_fromString_toList(self):

        gpr1 = '(A and B)'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr1)
        assert gpr_list == [['A','B']]

        gpr1 = '((A and B) or (C and D))'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr1)
        assert gpr_list == [['A','B'], ['C','D']]

        gpr1 = '((A and (B1 and B2)) or (C and D))'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr1)
        assert gpr_list == [['A','B1','B2'], ['C','D']]

        # GPR with three nested lists fails
        gpr1 = '((A and (B1 or B2)) or (C and D))'
        gpr_list = input1_manager.get_gpr_fromString_toList(gpr1)
        assert gpr_list != [['A',['B1','B2']], ['C','D']]

