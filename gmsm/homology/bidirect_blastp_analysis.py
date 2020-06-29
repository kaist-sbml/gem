
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import pickle
from gmsm.homology.blastp_utils import (
    make_blastDB,
    run_blastp,
    parseBlaspResults,
    makeBestHits_dict,
    getBBH,
    get_nonBBH
)

def get_homologs(io_ns, homology_ns):
    logging.info("Searching bidirectional best hits (BBHs)..")
    logging.info("Generating a DB for the target genome..")
    make_blastDB(io_ns)

    logging.info("Running BLASTP #1: target genome --> template model genes..")
    run_blastp(
            target_fasta='./%s/targetGenome_locusTag_aaSeq.fa' %io_ns.outputfolder2,
            blastp_result='./%s/blastp_targetGenome_against_tempGenome.txt'
            %io_ns.outputfolder2,
            db_dir = './%s/tempBlastDB.dmnd' %io_ns.outputfolder2)

    logging.info("Running BLASTP #2: template model genes --> target genome..")
    run_blastp(
            target_fasta=io_ns.temp_fasta,
            blastp_result='./%s/blastp_tempGenome_against_targetGenome.txt'
            %io_ns.outputfolder2,
            db_dir = './%s/targetBlastDB.dmnd' %io_ns.outputfolder2)
    
    logging.debug("Parsing the results of BLASTP #1..")
    blastpResults_dict1 = parseBlaspResults(
            './%s/blastp_targetGenome_against_tempGenome.txt'
            %io_ns.outputfolder2,
            './%s/blastp_targetGenome_against_tempGenome_parsed.txt'
            %io_ns.outputfolder2)

    logging.debug("Parsing the results of BLASTP #2..")
    blastpResults_dict2 = parseBlaspResults(
            './%s/blastp_tempGenome_against_targetGenome.txt'
            %io_ns.outputfolder2,
            './%s/blastp_tempGenome_against_targetGenome_parsed.txt'
            %io_ns.outputfolder2)

    logging.debug("Selecting the best hits for BLASTP #1..")
    bestHits_dict1 = makeBestHits_dict(
            './%s/blastp_targetGenome_against_tempGenome_parsed.txt'
            %io_ns.outputfolder2)

    logging.debug("Selecting the best hits for BLASTP #2..")
    bestHits_dict2 = makeBestHits_dict(
            './%s/blastp_tempGenome_against_targetGenome_parsed.txt'
            %io_ns.outputfolder2)

    logging.debug("Selecting the bidirectional best hits..")
    getBBH(bestHits_dict1, bestHits_dict2, homology_ns)

    logging.debug("Selecting genes that are not bidirectional best hits..")
    get_nonBBH(io_ns, homology_ns)

