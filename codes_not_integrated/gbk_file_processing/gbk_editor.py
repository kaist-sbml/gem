'''
2015 Hyun Uk Kim
'''

import sys
from Bio import SeqIO, Entrez

#Email is mandatory from June 1, 2010
#In case of excessive usage, NCBI will attempt to contact a user at the e-mail address provided prior to blocking access to the E-utilities
Entrez.email = "ehukim@kaist.ac.kr"

def read_gbk():

    GenbankFile = ""
    
    if len(sys.argv)>1:
        GenbankFile = sys.argv[1]
    else:
        ".gb(k) file is needed as an input"
        sys.exit(1)

    #Input and Output files can be located in the same directory: antiSMASH/EFICAz/
    if sys.argv[1].endswith('.gbk'):
        fileName = sys.argv[1][:-4]
    elif sys.argv[1].endswith('.gb'):
        fileName = sys.argv[1][:-3]
    outFile = "/data1/user_home/edhyunukkim/aSmodule_codes_not_integrated/gbk_file_processing/%s_trans.gbk" %fileName
    print outFile
    seq_records = SeqIO.read(GenbankFile, "genbank")

    return outFile, seq_records

def get_aa_sequence(outfile, seq_records):
    #Checks the presence of the qualifier "translation" and inserts one if absent
    out_recs = []

    for feature in seq_records.features:
        if feature.type == 'CDS':
            #Retrieving protein sequence from Entrez (protein) using protein_id
            print feature.qualifiers['locus_tag'][0]
	    if "translation" not in feature.qualifiers and "protein_id" in feature.qualifiers:
                print feature.qualifiers['protein_id'][0]
                try: 
	            handle = Entrez.efetch(db="protein", id=feature.qualifiers['protein_id'][0], rettype="fasta", retmode="text")
	            for protein_seq in SeqIO.parse(handle, 'fasta'):
	                feature.qualifiers['translation'] = protein_seq.seq
	            handle.close()
                except:
                    print "Cannot access Entrez Protein"
                    print "Writing output file with translations fetched up to current point"
                    SeqIO.write(seq_records, outfile, "genbank")

def count_cds_ec_translation(seq_records):
    num_cds = 0
    num_ec = 0
    num_trans = 0

    #Counts the number of CDS assigned with EC_number and translation
    for feature in seq_records.features:
        if feature.type == 'CDS':
            num_cds += 1
            if 'EC_number' in feature.qualifiers:
                num_ec += 1
            if 'translation' in feature.qualifiers:
                num_trans += 1
    
    print "\n", "number of CDS:", num_cds
    print "Number of CDS assigned with EC_number:", num_ec
    print "Number of CDS assigned with translation:", num_trans, "\n"

    return num_cds, num_ec, num_trans

