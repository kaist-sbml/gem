
#import filecmp
import logging
#import os
#import shutil
#import sys
import utils
#from os.path import join, abspath, dirname

__version__ = '0.7.3'


def check_prereqs(run_ns):
    "Check if all required files and applications are around"

    # Tuple is ( binary_name, optional)
    _required_binaries = [
        ('diamond', False),
        ('eficaz2.5', False)
    ]

    failure_messages = []

    for binary_name, optional in _required_binaries:
        binary_path = utils.locate_executable(binary_name)

        if binary_path:
            if binary_name == 'eficaz2.5':
                run_ns.eficaz_path = binary_path
        elif binary_path is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)
            if binary_name == 'eficaz2.5':
                run_ns.eficaz_path = binary_path

    try:
        import cobra
        version = cobra.__version__
        logging.debug("Found cobra version %s", version)

        # NOTE: Issue #333 was broken again in cobra >= 0.6.2:
        #https://github.com/opencobra/cobrapy/issues/333
        #Following fixation is automatically incorporated before running this code:
        #https://github.com/opencobra/cobrapy/commit/ac2f2e8bd31c982e87d5f385f7a8eea35e7ba811

        # NOTE: This issue has been fixed in cobra 0.8.2.
        #gmsm_dir = join(abspath(os.pardir), 'gmsm')
        #good_sbml_dir = join(gmsm_dir, 'scripts', 'bigg', 'sbml.py')
        #bad_sbml_dir = join(dirname(abspath(cobra.__file__)), 'io', 'sbml.py')

        #if filecmp.cmp(good_sbml_dir, bad_sbml_dir) == True:
        #    logging.debug("'sbml.py' in cobra is OK")
        #elif filecmp.cmp(good_sbml_dir, bad_sbml_dir) == False:
        #    shutil.copyfile(good_sbml_dir, bad_sbml_dir)
        #    logging.debug("'sbml.py' in cobra contains an error in SBML file writing.")
        #    logging.debug("Good version of 'sbml.py' has been copied in %s." %bad_sbml_dir)
        #    logging.error("Please re-run the program to have this fix to take effect.")
        #    sys.exit(0)

    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))

    try:
        import libsbml
        logging.debug("Found libsmbl version %s", libsbml.getLibSBMLVersion())
    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))

    for msg in failure_messages:
        logging.error(msg)

