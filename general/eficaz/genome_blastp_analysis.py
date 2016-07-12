'''
2014-2016
Hyun Uk Kim
'''

from Bio import SeqIO
import glob
import os
import pickle
import subprocess
import sys

def pickle_template_genes_list(dirname):
    fp1 = open('./%s/gap-filling_genes.txt' %dirname,"r")
    locusTag_list = []

    while 1:
        locusTag = fp1.readline()

        if not locusTag: break

	locusTag = locusTag.strip()
	print locusTag

	locusTag_list.append(locusTag)

    pickle.dump(locusTag_list, open('./%s/locusTag_list.p' %dirname,'wb'))
    fp1.close()


def get_template_genomeInfo(dirname, FileType, locusTag_list):
    fp = open('./%s/genome_locusTag_aaSeq.fa' %dirname,'w')
    genome_locusTag_aaSeq_dict = {}
    genome_locusTag_ec_dict = {}
    genome_locusTag_prod_dict = {}

    gbk_files = glob.glob("./%s/*.gb" %dirname)
    print "gbk files identified:"
    print gbk_files

    #Reads GenBank file
    for gbk_file in gbk_files:
        record = SeqIO.read(gbk_file, FileType)
        for feature in record.features:
            if feature.type == 'CDS':

                #Retrieving "locus_tag (i.e., ORF name)" for each CDS
                locusTag = feature.qualifiers['locus_tag'][0]

                if locusTag in locusTag_list:

                    #Some locus_tag's have multiple same qualifiers (e.g., EC_number)
	            for item in feature.qualifiers:

                        #Note that the numbers of CDS and "translation" do not match.
                        #There are occasions that CDS does not have "translation".
                        if item == 'translation':

                            #Retrieving "translation (i.e., amino acid sequences)" for each CDS
                            translation = feature.qualifiers.get('translation')
                            genome_locusTag_aaSeq_dict[locusTag] = translation[0]
                            print >>fp, '>%s\n%s' % (str(locusTag), str(translation[0]))

                        #Used to find "and" relationship in the GPR association
	                if item == 'product':
		            product = feature.qualifiers.get('product')[0]
		            genome_locusTag_prod_dict[locusTag] = product

                        #Watch multiple EC_number's
	     	        if item == 'EC_number':
	                    ecnum = feature.qualifiers.get('EC_number')
		            genome_locusTag_ec_dict[locusTag] = ecnum

    #Check if the gbk file has EC_number
    #Additional conditions should be given upon setup of in-house EC_number assigner
    print "\n", "len(genome_locusTag_ec_dict.keys):"
    print len(genome_locusTag_ec_dict)
    print "len(genome_locusTag_prod_dict.keys):"
    print len(genome_locusTag_prod_dict)

    if len(genome_locusTag_ec_dict) == 0:
	print "Error: no EC_number in the submitted gbk file"
        #FIXME: Don't use sys.exit
	sys.exit(1)
    fp.close()
    return genome_locusTag_ec_dict

 
def get_genome_info_from_gbk(dirname, FileType):
    fp = open('./%s/genome_locusTag_aaSeq.fa' %dirname,'w')
    genome_locusTag_aaSeq_dict = {}
    genome_locusTag_ec_dict = {}
    genome_locusTag_prod_dict = {}

    gbk_files = glob.glob("./%s/*.gb" %dirname)
    print "gbk files identified:"
    print gbk_files

    #Reads GenBank file
    for gbk_file in gbk_files:
        record = SeqIO.read(gbk_file, FileType)
        for feature in record.features:
            if feature.type == 'CDS':

                #Retrieving "locus_tag (i.e., ORF name)" for each CDS
                locusTag = feature.qualifiers['locus_tag'][0]

                #Some locus_tag's have multiple same qualifiers (e.g., EC_number)
	        for item in feature.qualifiers:

                    #Note that the numbers of CDS and "translation" do not match.
                    #There are occasions that CDS does not have "translation".
                    if item == 'translation':

                        #Retrieving "translation (i.e., amino acid sequences)" for each CDS
                        translation = feature.qualifiers.get('translation')
                        genome_locusTag_aaSeq_dict[locusTag] = translation[0]
                        print >>fp, '>%s\n%s' % (str(locusTag), str(translation[0]))

                    #Used to find "and" relationship in the GPR association
	            if item == 'product':
		        product = feature.qualifiers.get('product')[0]
		        genome_locusTag_prod_dict[locusTag] = product

                    #Watch multiple EC_number's
	     	    if item == 'EC_number':
	                ecnum = feature.qualifiers.get('EC_number')
		        genome_locusTag_ec_dict[locusTag] = ecnum

    #Check if the gbk file has EC_number
    #Additional conditions should be given upon setup of in-house EC_number assigner
    print "\n", "len(genome_locusTag_ec_dict.keys):"
    print len(genome_locusTag_ec_dict)
    print "len(genome_locusTag_prod_dict.keys):"
    print len(genome_locusTag_prod_dict)

    if len(genome_locusTag_ec_dict) == 0:
	print "Error: no EC_number in the submitted gbk file"
        #FIXME: Don't use sys.exit
	sys.exit(1)
    fp.close()
    return genome_locusTag_ec_dict 


#Looks for .fa and .gbk  files in the pre-defined folder
def get_fasta(dirname):
    for root, _, files in os.walk('./%s' %dirname):
        for f in files:
            if f.endswith('.fa'):
                fasta = os.path.join(root, f)
	        return fasta
	else:
            #FIXME: Don't use sys.exit
	    sys.exit(1)

#making database files using fasta files
def make_blastDB(dirname, query_fasta):
    db_dir = './%s/blastDB' %dirname
    DBprogramName = './blastpfiles/makeblastdb.exe'
    subprocess.call([DBprogramName,'-in',query_fasta,'-out',db_dir,'-dbtype','prot'])  

    #Checks if DB is properly created; otherwise shutdown
    if os.path.isfile('./%s/blastDB.psq' %dirname) == False:
	print "error in make_blastDB: blast DB not created"
        #FIXME: Don't use sys.exit
	sys.exit(1)


#Output: b0002,ASPK|b0002,0.0,100.00,820
#"1e-30" is set as a threshold for bidirectional best hits
def run_blastp(target_fasta = '', blastp_result = '', db_dir = '', evalue = 1e-30):
    BLASTPprogramName = './blastpfiles/blastp.exe'  
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
def getBBH(evalue_input, dic1,dic2, dirname1, dirname2):
    fp1 = open("./blastp_results_evalue_%s/bidirectional_best_hits1.txt" %evalue_input, 'w')
    fp2 = open("./blastp_results_evalue_%s/bidirectional_best_hits2.txt" %evalue_input, 'w')

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

                        #Some genes in template model have more than one homologous gene 
                        #in a target genome
		        if temp_locusTag not in temp_target_BBH_dict.keys():
			    temp_target_BBH_dict[temp_locusTag] = ([target_locusTag])
			else:
			    temp_target_BBH_dict[temp_locusTag].append((target_locusTag))

    fp1.write("locus_tag in %s" %dirname1+"\t"+"locus_tag(s) in %s" %dirname2+"\n")
    fp2.write("locus_tag in %s" %dirname1+"\t"+"locus_tag(s) in %s" %dirname2+"\n")

    for temp_locustag, target_locustag in temp_target_BBH_dict.iteritems():
        print temp_locustag, target_locustag
        print >>fp1, '%s\t%s' %(temp_locustag, target_locustag)

    for temp_locustag, target_locustag in temp_target_BBH_dict.iteritems():
        for target_gene in target_locustag:
            print >>fp2, '%s\t%s' %(temp_locustag, target_gene)

    fp1.close()
    fp2.close()
    return targetBBH_list


#A set of locusTag not included in BBH_list were considered nonBBH_list.
#Their respective reactions, if available, are added to the model in augPhase.
def get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list):
    nonBBH_list = []

    for locusTag in targetGenome_locusTag_ec_dict.keys():
	if locusTag not in targetBBH_list:
	    nonBBH_list.append(locusTag)

    nonBBH_list = sorted(set(nonBBH_list))
    return nonBBH_list

