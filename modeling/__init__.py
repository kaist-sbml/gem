
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
from io.eficaz import utils

def check_prereqs():
    "Check if all required files and applications are around"

    # Tuple is ( binary_name, optional)
    _required_binaries = [
        ('blastp', False),
        ('makeblastdb', False),
        ('eficaz2.5', False)
    ]

    failure_messages = []

    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    try:
        import cobra
        logging.debug("Found cobra version %s", cobra.__version__)
    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))

    try:
        import libsbml
        logging.debug("Found libsmbl version %s", libsbml.getLibSBMLVersion())
    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))

    for msg in failure_messages:
        logging.error(msg)

