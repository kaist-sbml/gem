'''
2015 Hyun Uk Kim
'''

from gbk_editor import read_gbk, get_aa_sequence, count_cds_ec_translation
from Bio import SeqIO
import time

start = time.time()

outfile, seq_records = read_gbk()
get_aa_sequence(outfile, seq_records)
num_cds, num_ec, num_trans = count_cds_ec_translation(seq_records)

#Retrieving aa sequence often skips some CDSs
while num_cds - num_trans > 0:
    get_aa_sequence(outfile, seq_records)
    num_cds, num_ec, num_trans = count_cds_ec_translation(seq_records)

#Writing an updaged gbk file
print "Writing output file with translation"
SeqIO.write(seq_records, outfile, "genbank")

print "elpased time:", time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
