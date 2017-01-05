
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import os
import subprocess


#Make database files using fasta files
def make_blastDB(options):
    db_dir = './%s/targetBlastDB' %options.outputfolder2
    DBprogramName = './gems/homology/blastpfiles/makeblastdb.exe'
    subprocess.call([DBprogramName,'-in',options.target_fasta,'-out',db_dir,'-dbtype','prot'])

    #Checks if DB is properly created; otherwise shutdown
    if os.path.isfile('./%s/targetBlastDB.psq' %options.outputfolder2) == False:
	logging.debug("Error in make_blastDB: blast DB not created")


#Output: b0002,ASPK|b0002,0.0,100.00,820
#"1e-30" is set as a threshold for bidirectional best hits
def run_blastp(target_fasta = '', blastp_result = '', db_dir = '', evalue = 1e-30):
    BLASTPprogramName = './gems/homology/blastpfiles/blastp.exe'
    subprocess.call([BLASTPprogramName,'-query',target_fasta,'-out',blastp_result,'-db',db_dir,'-evalue', str(evalue),'-outfmt',"10 qseqid sseqid evalue score length pident"])


#Input: Results file from "run_blastp"
#Output: '\t' inserted between each element of the input
def parseBlaspResults(inputFile, outputFile):
    blastpResults_dict = {}
    fp = open(inputFile,'r')
    fp2 = open(outputFile,'w')
    itemnum=0
    for line in fp:
        key = itemnum
        sptList = line.strip().split(',')
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
def getBBH(dic1, dic2, options):
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

    options.targetBBH_list = targetBBH_list
    options.temp_target_BBH_dict = temp_target_BBH_dict


#A set of locusTag not included in BBH_list were considered nonBBH_list.
#Their respective reactions, if available, are added to the model in augPhase.
#def get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list):
def get_nonBBH(options):
    nonBBH_list = []

    for locusTag in options.targetGenome_locusTag_ec_dict.keys():
	if locusTag not in options.targetBBH_list:
            nonBBH_list.append(locusTag)

    nonBBH_list = sorted(set(nonBBH_list))
    options.nonBBH_list = nonBBH_list

