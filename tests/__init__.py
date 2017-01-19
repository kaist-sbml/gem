
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import os
try:
    import pytest
except ImportError:
    pytest = None

def test_gems():
    if pytest:
        pytest.main(['--pyargs', 'gems'] )
    else:
        raise ImportError("pytest is not installed")
