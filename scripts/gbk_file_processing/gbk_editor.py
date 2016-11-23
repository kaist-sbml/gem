'''
2015 Hyun Uk Kim

Objective: To handle non-standard .gbk files by inserting "translation" from NCBI
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
    outFile = "/data1/user_home/edhyunukkim/antismash_modeling/codes_not_integrated/gbk_file_processing/%s_trans.gbk" %fileName
    print outFile
    seq_records = SeqIO.read(GenbankFile, "genbank")

    return fileName, outFile, seq_records

def get_aa_sequence(outfile, seq_records):
    #Checks the presence of the qualifier "translation" and inserts one if absent
    out_recs = []

    for feature in seq_records.features:
        if feature.type == 'CDS':
            #Retrieving protein sequence from Entrez (protein) using protein_id
	    if "translation" not in feature.qualifiers:
                if "protein_id" in feature.qualifiers:
                    print feature.qualifiers['locus_tag'][0], feature.qualifiers['protein_id'][0]
                    try: 
	                handle = Entrez.efetch(db="protein", id=feature.qualifiers['protein_id'][0], rettype="fasta", retmode="text")
	                for protein_seq in SeqIO.parse(handle, 'fasta'):
	                    feature.qualifiers['translation'] = protein_seq.seq
	                handle.close()
                    except:
                        print "Cannot access Entrez Protein"
                        print "Writing output file with translations fetched up to current point"
                        SeqIO.write(seq_records, outfile, "genbank")

                else:
                    print feature.qualifiers['locus_tag'][0], ": This CDS does not have protein_id"

def count_cds_qualifiers(seq_records):
    num_cds = 0
    num_ec = 0
    num_protid = 0
    num_trans = 0

    #Counts the number of CDS assigned with EC_number and translation
    for feature in seq_records.features:
        if feature.type == 'CDS':
            num_cds += 1
            if 'EC_number' in feature.qualifiers:
                num_ec += 1
            if 'protein_id' in feature.qualifiers:
                num_protid += 1
            if 'translation' in feature.qualifiers:
                num_trans += 1
    
    print "\n", "No. CDS:", num_cds
    print "No. CDS assigned with EC_number:", num_ec
    print "No. CDS assigned with protein_id:", num_protid
    print "No. CDS assigned with translation:", num_trans, "\n"

    return num_cds, num_ec, num_protid, num_trans

#This function was adopted from metabolicmodel / prunPhase.py
def get_targetGenomeInfo(seq_records, fileName):
    fp1 = open('./%s_ec.txt' %(fileName),'w')
    #fp2 = open('./%s_prod.txt' %(fileName),'w')
    targetGenome_locusTag_ec_dict = {}
    targetGenome_locusTag_prod_dict = {}

    ec_all = []
    ec3 = []
    ec4 = []

    for feature in seq_records.features:
        if feature.type == 'CDS':

            #Retrieving "locus_tag (i.e., ORF name)" for each CDS
            locusTag = feature.qualifiers['locus_tag'][0]

            #Some locus_tag's have multiple same qualifiers (e.g., EC_number)
	    for item in feature.qualifiers:

                #Used to find "and" relationship in the GPR association
	        if item == 'product':
		    product = feature.qualifiers.get('product')[0]
		    targetGenome_locusTag_prod_dict[locusTag] = product

                #Watch multiple EC_number's
		if item == 'EC_number':
	            ecnums = feature.qualifiers.get('EC_number')
		    targetGenome_locusTag_ec_dict[locusTag] = ecnums
                    for ecnum in ecnums:
                        ec_all.append(ecnum)
                        if ecnum[-1] != '-':
                            ec4.append(ecnum)
                        else:
                            ec3.append(ecnum)

    #Check if the gbk file has EC_number
    #Additional conditions should be given upon setup of in-house EC_number assigner
    print "len(targetGenome_locusTag_ec_dict.keys):"
    print len(targetGenome_locusTag_ec_dict)
    print "len(targetGenome_locusTag_prod_dict.keys):"
    print len(targetGenome_locusTag_prod_dict)

    print "\n", "No. all EC_number:", len(ec_all)
    print "No. EC_number with 3 numbers:", len(ec3)
    print "No. EC_number with 4 numbers:", len(ec4), "\n"

    for locus_tag in targetGenome_locusTag_ec_dict.keys():
        print >>fp1, '%s\t%s\t%s' %(locus_tag, targetGenome_locusTag_prod_dict[locus_tag], targetGenome_locusTag_ec_dict[locus_tag])

    #for locus_tag in targetGenome_locusTag_prod_dict.keys():
    #    print >>fp2, '%s\t%s' %(locus_tag, targetGenome_locusTag_prod_dict[locus_tag])

    fp1.close()
    #fp2.close()
