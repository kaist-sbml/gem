
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
    logging.info("Searching bidirectional homolgs..")
    logging.info("Generating a DB for the genes from the target genome..")
    make_blastDB(options)

    logging.info("Running BLASTP #1: genes in the target genome against genes in the template model..")
    run_blastp(target_fasta='./%s/1_blastp_results/targetGenome_locusTag_aaSeq.fa' %options.output, blastp_result='./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, db_dir = '%s/tempBlastDB' %(options.input1), evalue=1e-30)

    logging.info("Running BLASTP #2: genes in the template model against genes in the target genome..")
    run_blastp(target_fasta=options.temp_fasta, blastp_result='./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, db_dir = './%s/1_blastp_results/targetBlastDB' %options.output, evalue=1e-30)

    logging.info("Parsing the results of BLASTP #1..")
    blastpResults_dict1 = parseBlaspResults('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome.txt' %options.output, './%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)

    logging.info("Parsing the results of BLASTP #2..")
    blastpResults_dict2 = parseBlaspResults('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome.txt' %options.output, './%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

    logging.info("Selecting the best hits for BLASTP #1..")
    bestHits_dict1 = makeBestHits_dict('./%s/1_blastp_results/blastp_targetGenome_against_tempGenome_parsed.txt' %options.output)

    logging.info("Selecting the best hits for BLASTP #2..")
    bestHits_dict2 = makeBestHits_dict('./%s/1_blastp_results/blastp_tempGenome_against_targetGenome_parsed.txt' %options.output)

    logging.info("Selecting the bidirectional best hits..")
    getBBH(bestHits_dict1, bestHits_dict2, options)

    logging.info("Selecting genes that are not bidirectional best hits..")
    get_nonBBH(options)

