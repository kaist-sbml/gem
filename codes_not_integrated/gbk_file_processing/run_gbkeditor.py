'''
2015 Hyun Uk Kim
'''

from gbk_editor import read_gbk, get_aa_sequence, count_cds_qualifiers, get_targetGenomeInfo
from Bio import SeqIO
import time

start = time.time()

fileName, outfile, seq_records = read_gbk()
#get_aa_sequence(outfile, seq_records)
num_cds, num_ec, num_protid, num_trans = count_cds_qualifiers(seq_records)
get_targetGenomeInfo(seq_records, fileName)

#Retrieving aa sequence often skips some CDSs
#while num_cds - num_trans > 5:
#    get_aa_sequence(outfile, seq_records)
#    num_cds, num_ec, num_protid, num_trans = count_cds_qualifiers(seq_records)

#Writing an updaged gbk file
#print "Writing output file with translation"
#SeqIO.write(seq_records, outfile, "genbank")

print "elpased time:", time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
