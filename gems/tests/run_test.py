
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

import pytest

@pytest.mark.skip(reason = "This is a command for testing")
def run_test_gems():

    from __init__ import test_gems
    test_gems()

if __name__ == '__main__':
    run_test_gems()
