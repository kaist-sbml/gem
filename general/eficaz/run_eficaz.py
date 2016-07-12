#!/usr/bin/env python
'''
2016 Hyun Uk Kim

Objective: To calculate EC_numbers for the genes in FASTA format
'''

#import helperlibs.bio.seqio
#import antismash.generic_modules.ecpredictor
#import antismash.config
from eficaz import getECs 
import logging
import sys
from Bio import SeqIO 
from argparse import Namespace
#from antismash import utils
import os

fasta_file = sys.argv[1]

#Remove everything except for locustag ID in the header line
def read_fasta(fasta_file):
    fp1 = open('targetGenome_locusTag_aaSeq.fa', 'w')

    fasta_dict = {}

    for locustag in fasta_dict.keys():
        print >>fp1, '>%s\n%s' % (locustag, fasta_dict[locustag].seq)

    return fasta_dict

def get_targetGenome_locusTag_ec_dict(fasta_file):
    seq_record = fasta_dict
    options = Namespace()
#    antismash.config.load_config(options)
#    options.statusfile = "test.status"
    options.outputfoldername = "."
    options.cpus = 8

#    outputfoldername = "."
#    cpus = 8
    
    options.ecpred = "eficaz"

    #seq_record contains all the CDS features within, hence for loop is not necessary
    seq_record = list(SeqIO.parse(fasta_file, "fasta"))
    print len(seq_record)
    print type(seq_record)
    print seq_record.id
#    for seq_record in SeqIO.parse(fasta_file, "fasta"):
#        getECs(seq_record, options)


if __name__ == "__main__":
   fasta_dict = read_fasta(fasta_file)
   get_targetGenome_locusTag_ec_dict(fasta_file)

