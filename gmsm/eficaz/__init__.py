# vim: set fileencoding=utf-8 :
#
#
# Copyright (C) 2014 Tilmann Weber, Hyun Uk Kim
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


"""
Predict EC numbers by EFICAz 2.5a
"""

import glob
import logging
import multiprocessing
import os
import re
import shutil
import sys
import tempfile
import time
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from gmsm import utils


name = "eficaz"
short_description = name.capitalize()

EFICAzBinary = "eficaz2.5"

class EFICAzECPrediction:

    def __init__(self, run_ns, io_ns, inputfile):

        # Assign variables
        self.inputfile = inputfile
        self.run_ns = run_ns
        self.io_ns = io_ns

        # Variables to store EC prediction
        self.EC4Dict = {}
        self.EC3Dict = {}
        self.EC4InfoDict = {}
        self.EC3InfoDict = {}

        # Dictionary to store the fasta file per Chunkdir
        # self.ChunkFilenames['/path/to/Chunk1'] = '/path/to/fastafile'
        self.ChunkFilenames = {}

        self.basedirName = os.path.abspath(os.path.join(io_ns.outputfolder1, "EFICAz"))
        try:
            os.makedirs(self.basedirName)
        except OSError:
            if os.path.exists(self.basedirName):
                # We are nearly safe
                logging.debug("Outputfolder %s already exists." % self.basedirName)
            else:
                logging.exception("Cannot create EFICAz output directory %s" % self.basedirName)
                sys.exit(1)

        tempdir = tempfile.mkdtemp(prefix='antiSMASH_ECpred')
        self.tempdirname = tempdir

    def _getMultiFastaList(self, io_ns):
        allFastaList = []

        for gene_id, fasta_seq in io_ns.targetGenome_locusTag_aaSeq_dict.iteritems():
            if "-" in str(fasta_seq):
                fasta_seq = Seq(str(fasta_seq).replace("-",""), generic_protein)

            # Never write empty fasta entries
            elif len(fasta_seq) == 0:
                logging.debug("No translation for %s, skipping" % gene_id)
                continue

            allFastaList.append(">%s\n%s\n" % (gene_id, fasta_seq))

        return allFastaList


    def _prepareInput(self, run_ns, io_ns):
        """Generate $run_ns.cpus chunks of multi-Fasta-files; each in it's own subdirectory named "Chunk0000x";
        returns: list of directorynames"""

        logging.debug("Preparing input files for EFICAz")
        InputDirList = []
        allFastaList = self._getMultiFastaList(io_ns)

        logging.debug("Number of sequences: %s", len(allFastaList))

        maxChunks = self.run_ns.cpus


        if len(allFastaList) < maxChunks:
            maxChunks = len(allFastaList)

        if maxChunks == 0:
            logging.warn('No input files for %s', self.inputfile)
            return []
        equalpartsizes = int(len(allFastaList) / maxChunks)

        # Generate directory structure and write chunks;
        # for debug purposes use outputfolder; later on we should move to a temporary directory

        for i in range(maxChunks):
            if i == 0:
                fastaChunk = allFastaList[:equalpartsizes]
            elif i == (self.run_ns.cpus-1):
                fastaChunk = allFastaList[(i*equalpartsizes):]
            else:
                fastaChunk = allFastaList[(i*equalpartsizes):((i+1)*equalpartsizes)]


            # setup separate directories for EFICAz
            chunkDirName = "{basedir}{sep}Chunk{chunk_no:05d}".format(basedir=self.tempdirname, sep=os.sep,chunk_no=i+1)
            # logging.debug("Trying to create folder: %s" % chunkDirName)
            try:
                os.mkdir(chunkDirName)
            except OSError:
                if os.path.exists(chunkDirName):
                    # We are nearly safe
                    logging.debug("Outputfolder %s already exists." % chunkDirName)
                else:
                    logging.exception("Cannot create directory %s." % chunkDirName)
                    sys.exit(1)
            InputDirList.append(chunkDirName)

            chunkFileName = "{dirname}{sep}{seqid}_{chunk_no:05d}.fasta".format(
                    dirname=chunkDirName, sep=os.sep, seqid=self.inputfile, chunk_no=i+1)
            try:
                f  = open(chunkFileName, "w")
            except OSError:
                logging.exception("Cannot create fasta file %s" % chunkFileName)
                sys.exit(1)

            self.ChunkFilenames[chunkDirName] = os.path.abspath(chunkFileName)
            for seq in fastaChunk:
                f.write(seq)
            f.close()
        self.InputDirList = InputDirList
        return InputDirList

    def _runEFICAz(self, chunkDir):
        cwd = os.getcwd()
        try:
            os.chdir(chunkDir)
        except OSError:
            logging.exception("Can't chdir to %s" % chunkDir)
            sys.exit(1)

        fastafile = os.path.basename(self.ChunkFilenames[chunkDir])
        ecpredfile = fastafile+".ecpred"
        # Only perform calculations if result file does not already exist (from previous run)
        if not os.path.isfile(os.path.join(self.basedirName, ecpredfile)):
            EFICAzExecutable = utils.locate_executable(EFICAzBinary)
            cmdline = [EFICAzExecutable, fastafile]
            logging.debug("Executing %s in directory %s" % (" ".join(cmdline), chunkDir))

            try:
                utils.execute(cmdline)
            except:
                logging.exception('Cannot execute EFICAz!')
                sys.exit(1)
        else:
            # As this method is executed in an own thread, it does not have the ability to change
            # the variables within th eobject;
            # As a workaround we just copy the "old" file to the tempdir...
            try:
                shutil.copy(os.path.abspath(os.path.join(self.basedirName, ecpredfile)), self.ChunkFilenames[chunkDir]+".ecpred")
            except:
                logging.exception(
                        "Could not copy existing eficaz result file %s to tempfile %s",
                        os.path.isfile(os.path.abspath(self.basedirName, ecpredfile)),
                        self.ChunkFilenames[chunkDir]+".ecpred" )
                sys.exit(1)

        os.chdir(cwd)


    def _execute_EFICAz_processes(self, directorynames):

        processList = []

        for directoryname in directorynames:
            processList.append(multiprocessing.Process(target=self._runEFICAz, args = (directoryname, )))

        for process in processList:
            process.start()
        time.sleep(10)
        while True:
            processrunning = "n"
            for process in processList:
                if process.is_alive():
                    processrunning = "y"
            if processrunning == "y":
                time.sleep(5)
            else:
                break
        for process in processList:
            process.join()


        logging.debug("After joining the processes EC4Dict has %s entries" % len(self.EC4Dict.keys()))

    def _parseEFICAzResults(self, chunkDirs):

        for chunkDir in chunkDirs:

            # logging.debug("ChunkFilenames[%s]=%s", chunkDir, self.ChunkFilenames[chunkDir])
            ecpredfile = self.ChunkFilenames[chunkDir]+".ecpred"
            try:
                f = open(ecpredfile,"r")
            except OSError as e:
                logging.error("No EFICAz outputfile %s found. Skipping Chunk...\nOSError: %s", ecpredfile, e)
                continue
            except IOError as e:
                logging.error("No EFICAz outputfile %s found. Skipping Chunk...\nIOError: %s", ecpredfile, e)
                continue
            EC4Pred = {}
            EC4Info = {}
            EC3Pred = {}
            EC3Info = {}

            for line in f.read().splitlines():
                # First get antiSMASH-ID
                (antiSMASH_ID, eficazResultString) = line.split(',', 1)
                eficazResultString = eficazResultString.strip()
                if eficazResultString == 'No EFICAz EC assignment':
                    #logging.debug("No EC assignment found for %s" % antiSMASH_ID)
                    continue

                if eficazResultString.strip().startswith("3EC"):
                    #logging.debug("3EC: %s" % eficazResultString)
                    r = re.match('3EC: (\d+\.\d+\.\d+), (.*)', eficazResultString)
                    if r:
                        EC = r.group(1) + ".-"
                        ECDesc = r.group(2)
                        if not EC3Pred.has_key(antiSMASH_ID):
                            EC3Pred[antiSMASH_ID] = []
                            EC3Info[antiSMASH_ID] = []
                        EC3Pred[antiSMASH_ID].append(EC)
                        EC3Info[antiSMASH_ID].append(ECDesc)
                        continue

                if eficazResultString.strip().startswith("4EC"):
                    r = re.match('4EC: (\d+\.\d+\.\d+\.\d+), (.*)', eficazResultString)
                    if r:
                        EC = r.group(1)
                        ECDesc = r.group(2)
                        if not EC4Pred.has_key(antiSMASH_ID):
                            EC4Pred[antiSMASH_ID] = []
                            EC4Info[antiSMASH_ID] = []
                        EC4Pred[antiSMASH_ID].append(EC)
                        EC4Info[antiSMASH_ID].append(ECDesc)
                        continue

                logging.warn("Could not parse line %s:" % line)
            f.close()

            self.EC4Dict.update(EC4Pred)
            self.EC4InfoDict.update(EC4Info)
            self.EC3Dict.update(EC3Pred)
            self.EC3InfoDict.update(EC3Info)

            logging.debug("EC4Pred has %s entries for chunk" % len(self.EC4Dict.keys()))

    def _copyFiles(self, chunkDirs):
        "Copy the input and output files into outputfolder"

        logging.debug(
            "Copying the eficaz input/result files from tempdir %s to outputfolder %s",
            self.tempdirname, self.basedirName)

        for chunkDir in chunkDirs:
            try:
                # logging.debug("Copying input fasta file from %s to outputfolder", chunkDir)
                shutil.copy(self.ChunkFilenames[chunkDir], self.basedirName)
            except:
                logging.error("Could not copy eficaz input file %s to destination %s", \
                              self.ChunkFilenames[chunkDir], self.basedirName)

            try:
                # logging.debug("Copying results file from %s to outputfolder", chunkDir)
                shutil.copy(self.ChunkFilenames[chunkDir]+".ecpred", self.basedirName)
            except:
                logging.error("Could not copy eficaz result file %s to destination %s", \
                              self.ChunkFilenames[chunkDir]+".ecpred", self.basedirName)
        # And finally remove temporary directory
        logging.debug("Removing temp dir %s", self.tempdirname)
        shutil.rmtree(self.tempdirname)


    def runECpred(self, run_ns, io_ns):
        "Runs the EFICAz EC number predictions"
        chunkDirs = self._prepareInput(run_ns, io_ns)
        if len(chunkDirs) > 0:
            logging.debug("Split inputs to %s directories; first one is %s" % (len(chunkDirs), chunkDirs[0]))
            self._execute_EFICAz_processes(chunkDirs)
            self._parseEFICAzResults(chunkDirs)
            self._copyFiles(chunkDirs)
        else:
            logging.warn("ECpredictor: No protein coding sequences found for in record: %s" % self.inputfile)

    def getEC3(self, antiSMASH_ID):
        """Return list of EC3 numbers for antiSMASH_ID"""

        if self.EC3Dict.has_key(antiSMASH_ID):
            return self.EC3Dict[antiSMASH_ID]
        else:
            return None

    def getEC3Info(self, antiSMASH_ID):
        """Return list of infos for EC3 number prediction for antiSMASH_ID"""

        if self.EC3InfoDict.has_key(antiSMASH_ID):
            return self.EC3InfoDict[antiSMASH_ID]
        else:
            return None

    def getEC4(self, antiSMASH_ID):
        """Return list of EC4 numbers for antiSMASH_ID"""

        if self.EC4Dict.has_key(antiSMASH_ID):
            return self.EC4Dict[antiSMASH_ID]
        else:
            return None

    def getEC4Info(self, antiSMASH_ID):
        """Return list of infos for EC4 number prediction for antiSMASH_ID"""

        if self.EC4InfoDict.has_key(antiSMASH_ID):
            return self.EC4InfoDict[antiSMASH_ID]
        else:
            return None

    def getEC4Dict(self):
        """Return dictionary of list for 4-digit EC numbers

        Example:
        a = EFICAzObject.getEC4Dict
        will result in:
        a[antiSMASH_ID] = ['1.2.3.4', '5.6.7.8']"""

        return self.EC4Dict

    def getEC3Dict(self):
        """Return dictionary of 3-digit EC numbers

        Example:
        a = EFICAzObject.getEC3Dict
        will result in:
        a[antiSMASH_ID] = ['1.2.3.x', '5.6.7.x']"""

        return self.EC3Dict


    def getEC4InfoDict(self):
        """Return dictionary of description for 4-digit EC assignment

        Example:
        a = EFICAzObject.getEC4ToolDict
        will result in:
        a[antiSMASH_ID] = 'EFICAz_components: CHIEFc_SVM; PFAM_SVM, MTTSI_bin: 6, Precision (mean; SD): 0.991; 0.094' """

        return self.EC4InfoDict

    def getEC3InfoDict(self):
        """Return dictionary of description for 3-digit EC assignment

        Example:
        a = EFICAzObject.getEC3ToolDict
        will result in:
        a[antiSMASH_ID] = 'EFICAz_components: CHIEFc_SVM; PFAM_SVM, MTTSI_bin: 6, Precision (mean; SD): 0.991; 0.094' """

        return self.EC3InfoDict


def getECs1(run_ns, io_ns, seq_record):
    logging.debug("Predicting EC numbers with EFICAz (getECs1)")

    if not 'cpus' in run_ns:
        run_ns.cpus = 1

    inputfile = os.path.basename(run_ns.input).split('.')[0]
    EFICAzECs = EFICAzECPrediction(run_ns, io_ns, inputfile)
    EFICAzECs.runECpred(run_ns, io_ns)

    logging.debug("Found %s predictions for EC4" % len(EFICAzECs.getEC4Dict().keys()))
    logging.debug("Found %s predictions for EC3" % len(EFICAzECs.getEC3Dict().keys()))

    for feature in utils.get_cds_features(seq_record):
        featureID = utils.get_gene_id(feature)

        notes = []

        if feature.qualifiers.has_key("note"):
            notes = feature.qualifiers['note']

        if EFICAzECs.getEC4(featureID):
            logging.debug("Annotating %s" % featureID)
            if feature.qualifiers.has_key('EC_number'):
                logging.warn(
                    'ECpredictor[eficaz]: Overwriting existing EC annotation: %s with %s'
                    %(", ".join(feature.qualifiers['EC_number']),
                    ", ".join(EFICAzECs.getEC4(featureID)))
                    )

            feature.qualifiers['EC_number'] = EFICAzECs.getEC4(featureID)
            notes.append(
                    "EFICAz EC number prediction: EC4: {0}; {1}".format(
                    ", ".join(EFICAzECs.getEC4(featureID)),
                    "; ".join(EFICAzECs.getEC4Info(featureID))))

        # Only annotate 3 digit EC if no 4 digit EC is available
        if (EFICAzECs.getEC3(featureID) and not EFICAzECs.getEC4(featureID)):
            if feature.qualifiers.has_key('EC_number'):
                if not re.search(
                        "\d+\.\d+\.\d+\.\d+", " ".join(feature.qualifiers['EC_number'])):
                    logging.warn(
                        'ECpredictor[eficaz]: Overwriting existing EC annotation: %s  with %s'
                        %(", ".join(feature.qualifiers['EC_number']),
                        ", ".join(EFICAzECs.getEC3(featureID)))
                        )
                    feature.qualifiers['EC_number'] = EFICAzECs.getEC3(featureID)

        if EFICAzECs.getEC3Info(featureID):
            notes.append(
                "EFICAz EC number prediction: EC3: {0}; {1}".format(", ".join(
                EFICAzECs.getEC3(featureID)),
                "; ".join(EFICAzECs.getEC3Info(featureID)))
                )
            if not feature.qualifiers.has_key('EC_number'):
                feature.qualifiers['EC_number'] = EFICAzECs.getEC3(featureID)

        feature.qualifiers['note'] = notes
    logging.debug("Finished EC number prediction with EFICAz")

    #Write output gbk file
    output_gbk = inputfile + '_ec.gbk'
    SeqIO.write(seq_record, os.path.join(io_ns.outputfolder1, output_gbk), 'genbank')


def getECs2(run_ns, io_ns):
    logging.debug("Predicting EC numbers with EFICAz (getECs2)")

    if not 'cpus' in run_ns:
        run_ns.cpus = 1

    inputfile = os.path.basename(run_ns.input).split('.')[0]
    EFICAzECs = EFICAzECPrediction(run_ns, io_ns, inputfile)
    EFICAzECs.runECpred(run_ns, io_ns)

    logging.debug("Found %s predictions for EC4" % len(EFICAzECs.getEC4Dict().keys()))
    logging.debug("Found %s predictions for EC3" % len(EFICAzECs.getEC3Dict().keys()))

    for locustag in io_ns.targetGenome_locusTag_aaSeq_dict:
        if EFICAzECs.getEC4(locustag):
            ecnum = EFICAzECs.getEC4(locustag)

            if locustag not in io_ns.targetGenome_locusTag_ec_dict:
                io_ns.targetGenome_locusTag_ec_dict[locustag] = []

            for ec in ecnum:
                if ec not in io_ns.targetGenome_locusTag_ec_dict[locustag]:
                    io_ns.targetGenome_locusTag_ec_dict[locustag].append(ec)

    logging.debug("len(io_ns.targetGenome_locusTag_ec_dict.keys): %s",
                  len(io_ns.targetGenome_locusTag_ec_dict.keys()))

    #Write output file
    ecpredfile = inputfile + '_ec.txt'
    outputfiles = glob.glob(os.path.abspath(os.path.join(io_ns.outputfolder1, "EFICAz", '*.ecpred')))

    with open(os.path.join(io_ns.outputfolder1, ecpredfile), 'w') as outfile:
        for outputfile in outputfiles:
            with open(outputfile) as infile:
                outfile.write(infile.read())

