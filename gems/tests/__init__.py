
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import warnings
from cobra.io import read_sbml_model, write_sbml_model
from os.path import join, abspath, dirname

try:
    import pytest
except ImportError:
    pytest = None

def test_gems():
    # Suppress warning messages from cobrapy when tests fail
    # The messages do not help the testing
    warnings.filterwarnings("ignore")

    if pytest:
        # Arguement 'gems' is needed to test only the GEMS files
        pytest.main(['--pyargs', 'gems', '--basetemp=tmp', '--cov=gems', '-v'])
    else:
        raise ImportError("pytest is not installed")
