
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

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
