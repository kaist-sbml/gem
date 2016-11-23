'''
2014-2016
Hyun Uk Kim
'''

from genome_blastp_analysis import (
    pickle_template_genes_list,
    get_new_old_locustag_from_gbk,
    get_template_genomeInfo,
    get_genome_info_from_gbk,
    get_fasta,
    make_blastDB,
    run_blastp,
    parseBlaspResults,
    makeBestHits_dict,
    getBBH,
    get_nonBBH
)
import os
import pickle
import sys
import time

start = time.time()

evalue_input = sys.argv[1]

dirname1 = 'genome_query' #template genome
dirname2 = 'genome_db' #target genome

#Create output folders
folder = 'blastp_results_evalue_%s' %evalue_input

if not os.path.isdir(folder):
   os.makedirs(folder)


###################################################################
#print "Creating and loading a pickle file having a list of locus_tag's to study.."
#pickle_template_genes_list(dirname1)
#locusTag_list = pickle.load(open('./%s/locusTag_list.p' %dirname1,'rb'))
###################################################################

###################################################################
#print "Extract individual information from gbk files"
#get_new_old_locustag_from_gbk(dirname2, "genbank")
###################################################################

                           
###################################################################
#print "\n", "Reading genbank file of the template genome.."
#tempGenome_locusTag_ec_dict = get_template_genomeInfo(dirname1, "genbank", locusTag_list)
#targetGenome_locusTag_ec_dict = get_genome_info_from_gbk(dirname1, "genbank")

print "\n", "Looking for a fasta file of a template genome.."
template_fasta = get_fasta(dirname1)

print "\n", "Generating a DB for the genes from the template genome.."
make_blastDB(dirname1, query_fasta=template_fasta)


print "\n", "Reading genbank file of the target genome.."
targetGenome_locusTag_ec_dict = get_genome_info_from_gbk(dirname2, "genbank")

print "\n", "Looking for a fasta file of a target genome..", "\n"
target_fasta = get_fasta(dirname2)

print "\n", "Generating a DB for the genes from the target genome.."
make_blastDB(dirname2, query_fasta=target_fasta)
###################################################################

###################################################################
print "\n", "Running BLASTP #1: genes in target genome against genes in template genome.."
run_blastp(target_fasta='./%s/genome_locusTag_aaSeq.fa' %dirname2, blastp_result='./blastp_results_evalue_%s/blastp_%s_against_%s.txt' %(evalue_input, dirname2, dirname1), db_dir = './%s/blastDB' %dirname1, evalue=evalue_input)

print "\n", "Running BLASTP #2: genes in template genome against genes in target genome.."
run_blastp(target_fasta='./%s/genome_locusTag_aaSeq.fa' %dirname1, blastp_result='./blastp_results_evalue_%s/blastp_%s_against_%s.txt' %(evalue_input, dirname1, dirname2), db_dir = './%s/blastDB' %dirname2, evalue=evalue_input)
###################################################################

###################################################################
print "\n", "Parsing the results of BLASTP #1.."
blastpResults_dict1 = parseBlaspResults('./blastp_results_evalue_%s/blastp_%s_against_%s.txt' %(evalue_input, dirname2, dirname1), './blastp_results_evalue_%s/blastp_%s_against_%s_parsed.txt' %(evalue_input, dirname2, dirname1))

print "\n", "Parsing the results of BLASTP #2.."
blastpResults_dict2 = parseBlaspResults('./blastp_results_evalue_%s/blastp_%s_against_%s.txt' %(evalue_input, dirname1, dirname2), './blastp_results_evalue_%s/blastp_%s_against_%s_parsed.txt' %(evalue_input, dirname1, dirname2))

print "\n", "Selecting the best hits for BLASTP #1.."
bestHits_dict1 = makeBestHits_dict('./blastp_results_evalue_%s/blastp_%s_against_%s_parsed.txt' %(evalue_input, dirname2, dirname1))

print "\n", "Selecting the best hits for BLASTP #2.."
bestHits_dict2 = makeBestHits_dict('./blastp_results_evalue_%s/blastp_%s_against_%s_parsed.txt' %(evalue_input, dirname1, dirname2))
###################################################################

###################################################################
print "\n", "Selecting the bidirectional best hits.."
targetBBH_list = getBBH(evalue_input, bestHits_dict1, bestHits_dict2, dirname1, dirname2)

#print "Selecting genes that are not bidirectional best hits in %s.." %dirname1
#nonBBH_list = get_nonBBH(tempGenome_locusTag_ec_dict, targetBBH_list)

#print "Selecting genes that are not bidirectional best hits in %s.." %dirname2
#nonBBH_list = get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list)
###################################################################

print "\n", "Elapsed time:", time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
