
from os.path import join
from gems.config import load_config
from gems.io.input_file_manager import get_eficaz_file, get_locustag_comp_dict

class TestInput_file_manager:
    """Test functions in gems.io.input_file_manager"""

    def test_get_eficaz_file(self, eficaz_file, options):

        options.eficaz_file = eficaz_file
        options.targetGenome_locusTag_ec_dict = {}
        get_eficaz_file(options)

        assert len(options.targetGenome_locusTag_ec_dict) == 2
        assert options.targetGenome_locusTag_ec_dict['NSK_00005-RA'] == ['2.7.1.83']


    def test_get_locustag_comp_dict(self, comp_file, options):

        options.comp = comp_file
        get_locustag_comp_dict(options)

        assert len(options.locustag_comp_dict) == 8
        assert options.locustag_comp_dict['NSK_00004-RA'] == 'h'

