
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import os
import subprocess
import multiprocessing

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
        print >>fp2, '%s\t%s\t%s\t%f\t%d\t%f' % (qseqid, sseqid, evalue, score, length, pident)
        itemnum += 1
    fp.close()
    fp2.close()
    return blastpResults_dict


#searching the best hit of a particular gene to a target genome
#Input: Results file from "parseBlaspResults"
#Output: query_locusTag '\t' db_locusTag
def makeBestHits_dict(inputFile):
    bestHits_dict = {}
    fp1 = open(inputFile,'r')

    for line in fp1:
        sptList = line.strip().split('\t')
        query_locusTag = sptList[0].strip()
        db_locusTag = sptList[1].strip()

        if query_locusTag not in bestHits_dict.keys():
            #Value is in List to enable "append" below
            bestHits_dict[query_locusTag] = ([db_locusTag])

        #This additional condition is necessary because blastp strangely produces
        #a redundant set of gene pairs:
        #e.g., SCO5892	['SAV_7184', 'SAV_419', 'SAV_419', 'SAV_419', 'SAV_419', ...]
        elif '%s' %(db_locusTag) not in bestHits_dict[query_locusTag]:
            bestHits_dict[query_locusTag].append((db_locusTag))

    fp1.close()
    return bestHits_dict


#Finding bidirectional best hits
#Input: two dict data from "selectBestHits" (e.g.,bestHits_dict)
def getBBH(dic1, dic2, homology_ns):
    targetBBH_list = []
    temp_target_BBH_dict = {}

    for target_locusTag in dic1.keys():
        for temp_locusTag in dic1[target_locusTag]:
            if temp_locusTag in dic2.keys():
                for target_locusTag2 in dic2[temp_locusTag]:

                    #The BBH case
                    if target_locusTag == target_locusTag2:
                        if target_locusTag not in targetBBH_list:
                            targetBBH_list.append(target_locusTag)

                        #Some genes in template model have more than one homologous
                        #in a target genome
                        if temp_locusTag not in temp_target_BBH_dict.keys():
                            temp_target_BBH_dict[temp_locusTag] = \
                            ([target_locusTag])
                        else:
                            temp_target_BBH_dict[temp_locusTag].append((target_locusTag))

    homology_ns.targetBBH_list = targetBBH_list
    homology_ns.temp_target_BBH_dict = temp_target_BBH_dict


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

