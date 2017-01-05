
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import pickle
from blastp_utils import (
    make_blastDB,
    run_blastp,
    parseBlaspResults,
    makeBestHits_dict,
    getBBH,
    get_nonBBH
)

def get_homologs(options):
    logging.info("Searching bidirectional best hits (BBHs)..")
    logging.info("Generating a DB for the target genome..")
    make_blastDB(options)

    logging.info("Running BLASTP #1: target genome --> template model genes..")
    run_blastp(
            target_fasta='./%s/targetGenome_locusTag_aaSeq.fa' %options.outputfolder2,
            blastp_result='./%s/blastp_targetGenome_against_tempGenome.txt'
            %options.outputfolder2,
            db_dir = '%s/tempBlastDB' %(options.input1),
            evalue=1e-30)

    logging.info("Running BLASTP #2: template model genes --> target genome..")
    run_blastp(
            target_fasta=options.temp_fasta,
            blastp_result='./%s/blastp_tempGenome_against_targetGenome.txt'
            %options.outputfolder2,
            db_dir = './%s/targetBlastDB' %options.outputfolder2,
            evalue=1e-30)

    logging.debug("Parsing the results of BLASTP #1..")
    blastpResults_dict1 = parseBlaspResults(
            './%s/blastp_targetGenome_against_tempGenome.txt'
            %options.outputfolder2,
            './%s/blastp_targetGenome_against_tempGenome_parsed.txt'
            %options.outputfolder2)

    logging.debug("Parsing the results of BLASTP #2..")
    blastpResults_dict2 = parseBlaspResults(
            './%s/blastp_tempGenome_against_targetGenome.txt'
            %options.outputfolder2,
            './%s/blastp_tempGenome_against_targetGenome_parsed.txt'
            %options.outputfolder2)

    logging.debug("Selecting the best hits for BLASTP #1..")
    bestHits_dict1 = makeBestHits_dict(
            './%s/blastp_targetGenome_against_tempGenome_parsed.txt'
            %options.outputfolder2)

    logging.debug("Selecting the best hits for BLASTP #2..")
    bestHits_dict2 = makeBestHits_dict(
            './%s/blastp_tempGenome_against_targetGenome_parsed.txt'
            %options.outputfolder2)

    logging.debug("Selecting the bidirectional best hits..")
    getBBH(bestHits_dict1, bestHits_dict2, options)

    logging.debug("Selecting genes that are not bidirectional best hits..")
    get_nonBBH(options)

