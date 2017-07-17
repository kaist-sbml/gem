
from os.path import join
from gmsm.config import load_config

class TestConfig:
    """Test functions in gmsm.config"""

    def test_load_config(self, options):

        _cfg_name = 'gmsm.cfg'
        load_config(options)

        assert options.blastp
        assert options.blastp.evalue
        assert options.cobrapy
        assert options.cobrapy.gapfill_iter
