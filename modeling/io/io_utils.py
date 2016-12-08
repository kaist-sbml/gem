
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import os
from Bio import SeqIO
from eficaz.__init__ import getECs


#Look for .xml and .gb(k) files in the pre-defined folder
def get_temp_fasta(options):
    for root, _, files in os.walk('./modeling/io/data/input1/%s/' %(options.orgName)):
        for f in files:
            if f.endswith('.fa'):
                tempFasta = os.path.join(root, f)
                options.input1 = root
                options.temp_fasta = tempFasta


def get_targetGenomeInfo(options, file_type):

    fp = open('./%s/1_blastp_results/targetGenome_locusTag_aaSeq.fa'
            %options.outputfolder,'w')

    targetGenome_locusTag_aaSeq_dict = {}
    targetGenome_locusTag_ec_dict = {}
    targetGenome_locusTag_prod_dict = {}

    #Read GenBank file
    try:
        seq_record = SeqIO.read(options.outputfolder+'/'+options.input, file_type)
    except ValueError:
        logging.debug("Warning: ValueError occurred in SeqIo.read")
        seq_record = SeqIO.parse(options.outputfolder+'/'+options.input,
                     file_type).next()

    if options.eficaz:
        getECs(seq_record, options)

    total_cluster = 0

    for feature in seq_record.features:
        if feature.type == 'CDS':

            #Retrieving "locus_tag (i.e., ORF name)" for each CDS
            locusTag = feature.qualifiers['locus_tag'][0]

            #Some locus_tag's have multiple same qualifiers (e.g., EC_number)
            for item in feature.qualifiers:

                #Note that the numbers of CDS and "translation" do not match.
                #There are occasions that CDS does not have "translation".

#These if statements may be removed:e.g., sec_met_rxn_generation.py
                if item == 'translation':

                    #Retrieving "translation (i.e., amino acid sequences)" for each CDS
                    translation = feature.qualifiers.get('translation')
                    targetGenome_locusTag_aaSeq_dict[locusTag] = translation[0]
                    print >>fp, '>%s\n%s' % (str(locusTag), str(translation[0]))

                #Used to find "and" relationship in the GPR association
                if item == 'product':
                    product = feature.qualifiers.get('product')[0]
                    targetGenome_locusTag_prod_dict[locusTag] = product

                #Watch multiple EC_number's
		if item == 'EC_number':
                    ecnum = feature.qualifiers.get('EC_number')
                    targetGenome_locusTag_ec_dict[locusTag] = ecnum

        if feature.type == 'cluster':
            total_cluster += 1

    #Check if the gbk file has EC_number
    #Additional conditions should be given upon setup of in-house EC_number assigner
    logging.debug("len(targetGenome_locusTag_ec_dict.keys): %s" %len(targetGenome_locusTag_ec_dict))
    logging.debug("len(targetGenome_locusTag_prod_dict.keys):%s" %len(targetGenome_locusTag_prod_dict))

    options.targetGenome_locusTag_ec_dict = targetGenome_locusTag_ec_dict
    options.targetGenome_locusTag_prod_dict = targetGenome_locusTag_prod_dict
    options.total_cluster = total_cluster
    options.seq_record = seq_record

    fp.close()


#Look for .fa and .gbk  files in the pre-defined folder
def get_target_fasta(options):
    for root, _, files in os.walk('./%s/1_blastp_results' %options.outputfolder):
        for f in files:
            if f.endswith('.fa'):
                target_fasta = os.path.join(root, f)
                options.target_fasta = target_fasta

    if not options.target_fasta:
        logging.warning("FASTA file %s for bidirectional blastp hits not found.")

