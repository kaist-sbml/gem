
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from os.path import join
from gems.config import load_config

class TestConfig:
    """Test functions in gems.config"""

    def test_load_config(self, options):

        _cfg_name = 'gems.cfg'
        load_config(options)

        assert options.blastp
        assert options.blastp.evalue
        assert options.cobrapy
        assert options.cobrapy.gapfill_iter
