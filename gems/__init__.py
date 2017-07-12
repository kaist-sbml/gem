
import logging
import os
import utils


__version__ = '0.3.3'


def check_prereqs(options):
    "Check if all required files and applications are around"

    # Tuple is ( binary_name, optional)
    _required_binaries = [
        ('blastp', False),
        ('makeblastdb', False),
        ('eficaz2.5', False)
    ]

    failure_messages = []

    for binary_name, optional in _required_binaries:
        binary_path = utils.locate_executable(binary_name)

        if binary_path:
            if binary_name == 'makeblastdb':
                options.makeblastdb_path = binary_path
            elif binary_name == 'blastp':
                options.blastp_path = binary_path
            elif binary_name == 'eficaz2.5':
                options.eficaz_path = binary_path
        elif binary_path is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)
            if binary_name == 'eficaz2.5':
                options.eficaz_path = binary_path

    try:
        import cobra
        version = cobra.__version__
        logging.debug("Found cobra version %s", version)
    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))

    try:
        import libsbml
        logging.debug("Found libsmbl version %s", libsbml.getLibSBMLVersion())
    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))

    for msg in failure_messages:
        logging.error(msg)

