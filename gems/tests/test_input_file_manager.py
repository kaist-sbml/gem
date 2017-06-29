
from os.path import join
from gems.config import load_config
from gems.io.input_file_manager import get_eficaz_file, get_locustag_comp_dict
from os.path import join, abspath, dirname

data_model_dir = join(dirname(abspath(__file__)), 'data')

class TestInput_file_manager:
    """Test functions in gems.io.input_file_manager"""

    def test_get_eficaz_file(self, options):

        options.eficaz_file = join(data_model_dir, 'NSK_all_genomes_ec.txt')
        options.targetGenome_locusTag_ec_dict = {}
        get_eficaz_file(options)

        assert len(options.targetGenome_locusTag_ec_dict) == 2
        assert options.targetGenome_locusTag_ec_dict['NSK_04159-RA'] == ['2.6.1.83']


    def test_get_locustag_comp_dict(self, options):

        options.comp = join(data_model_dir, 'Nanno_Compartment_result_dic_v3_test.txt')
        get_locustag_comp_dict(options)

        assert len(options.locustag_comp_dict) == 8
        assert options.locustag_comp_dict['NSK_00004-RA'] == 'h'

