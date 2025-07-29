
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import os
import subprocess
import multiprocessing
from collections import defaultdict

#Make database files using fasta files
def make_blastDB(io_ns):
    db_dir = './%s/targetBlastDB' %io_ns.outputfolder2
    subprocess.call("./bin/diamond makedb --in %s -d %s"%(io_ns.target_fasta, db_dir), shell=True, stderr=subprocess.STDOUT)
    
    #Checks if DB is properly created; otherwise shutdown
    if os.path.isfile('./%s/targetBlastDB.dmnd' %io_ns.outputfolder2) == False:
        logging.debug("Error in make_blastDB: blast DB not created")
    else:
        logging.debug("targetBlastDB.dmnd created")
        
    db_dir = './%s/tempBlastDB' %io_ns.outputfolder2
    subprocess.call("./bin/diamond makedb --in %s -d %s"%(io_ns.temp_fasta, db_dir), shell=True, stderr=subprocess.STDOUT)
    
    #Checks if DB is properly created; otherwise shutdown
    if os.path.isfile('./%s/tempBlastDB.dmnd' %io_ns.outputfolder2) == False:
        logging.debug("Error in make_blastDB: blast DB not created")
    else:
        logging.debug("tempBlastDB.dmnd created")


#Output: b0002,ASPK|b0002,0.0,100.00,820
#"1e-30" is set as a threshold for bidirectional best hits
def run_blastp(target_fasta, blastp_result, db_dir):
    subprocess.call("./bin/diamond blastp -d %s -q %s -o %s --evalue 1e-30 --id 30 --outfmt 6 qseqid sseqid evalue score length pident"%(db_dir, target_fasta, blastp_result), shell=True, stderr=subprocess.STDOUT)

#Input: Results file from "run_blastp"
#Output: '\t' inserted between each element of the input
def parseBlaspResults(inputFile, outputFile):
    blastpResults_dict = {}
    fp = open(inputFile,'r')
    fp2 = open(outputFile,'w')
    itemnum=0
    for line in fp:
        key = itemnum
        sptList = line.strip().split('\t')
        qseqid = sptList[0].strip()
        sseqid = sptList[1].strip()
        evalue = sptList[2].strip()
        score = float(sptList[3].strip())
        length = int(sptList[4].strip())
        pident = float(sptList[5].strip())

        blastpResults_dict[key] = {"query_locusTag": qseqid, "db_locusTag": sseqid, "evalue": evalue, "score": score, "length": length, "identity": pident}
        print('%s\t%s\t%s\t%f\t%d\t%f' % (qseqid, sseqid, evalue, score, length, pident), file=fp2)
        itemnum += 1
    fp.close()
    fp2.close()
    return blastpResults_dict


#searching the best hit of a particular gene to a target genome
#Input: tab-separated file from parseBlaspResults() (qseqid  sseqid  evalue  score  length  pident)
#Output: { qseqid: [sseqid1, sseqid2, …], … } (lists sorted by ascending e-value)
def makeBestHits_dict(inputFile):
    # 1) Keep only the smallest e-value for each (qseqid, sseqid) pair
    temp = defaultdict(dict)
    with open(inputFile) as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            qseqid, sseqid, evalue_str, *_ = line.split('\t')
            evalue = float(evalue_str)
            prev = temp[qseqid].get(sseqid)
            # if this pair is new or has a smaller e-value, update it
            if prev is None or evalue < prev:
                temp[qseqid][sseqid] = evalue

    # 2) Sort subjects by ascending e-value and extract sseqid lists
    bestHits_dict = {}
    for qseqid, subj_dict in temp.items():
        # subj_dict.items() -> [(sseqid, evalue), ...]
        sorted_subjs = sorted(subj_dict.items(), key=lambda x: x[1])
        # extract sorted sseqid list
        bestHits_dict[qseqid] = [s for s, _ in sorted_subjs]

    return bestHits_dict


#Finding bidirectional best hits
#Input: two dict data from "selectBestHits" (e.g.,bestHits_dict)
def getBBH(dic1, dic2, homology_ns):
    """
    dic1: query → [subject1, subject2, …]  
          (hits list, where the first element is the best hit)
    dic2: subject → [query1, query2, …]  
          (reverse hits list)
    homology_ns: namespace object in which to store results
    """
    target_bbh_list = []
    temp_target_bbh_dict = {}

    for query, subjects in dic1.items():
        if not subjects:
            continue
        # extract the best hit from the forward direction
        best_subject = subjects[0]

        # look up the reverse hits for that subject
        reverse_hits = dic2.get(best_subject, [])
        if not reverse_hits:
            continue

        # only accept as BBH if the top reverse hit matches the original query
        if reverse_hits[0] == query:
            # add to the list of BBH queries
            target_bbh_list.append(query)
            # record the reciprocal mapping for the subject
            temp_target_bbh_dict.setdefault(best_subject, []).append(query)

    homology_ns.targetBBH_list = target_bbh_list
    homology_ns.temp_target_BBH_dict = temp_target_bbh_dict


#A set of locusTag not included in BBH_list were considered nonBBH_list.
#Their respective reactions, if available, are added to the model in augPhase.
#def get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list):
def get_nonBBH(io_ns, homology_ns):
    nonBBH_list = []

    for locusTag in io_ns.targetGenome_locusTag_ec_dict.keys():
        if locusTag not in homology_ns.targetBBH_list:
            nonBBH_list.append(locusTag)

    nonBBH_list = sorted(set(nonBBH_list))
    homology_ns.nonBBH_list = nonBBH_list

