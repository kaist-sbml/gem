
from os.path import join
from gems.config import load_config
from gems.io.input_file_manager import get_eficaz_file, get_locustag_comp_dict
from os.path import join, abspath, dirname

class TestInput_file_manager:
    """Test functions in gems.io.input_file_manager"""

#    def test_get_eficaz_file(self):
        #TODO: Define


    def test_get_locustag_comp_dict(self, options):

        options.comp = \
                join(dirname(abspath(__file__)),
                     'data',
                     'Nanno_Compartment_result_dic_v3_test.txt')

        get_locustag_comp_dict(options)

        assert len(options.locustag_comp_dict) == 8
        assert options.locustag_comp_dict['NSK_00004-RA'] == 'h'

