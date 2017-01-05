
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import os
from Bio import SeqIO
from eficaz.__init__ import getECs


#Look for pre-stored fasta file of the template model
def get_temp_fasta(options):
    for root, _, files in os.walk('./gems/io/data/input1/%s/' %(options.orgName)):
        for f in files:
            if f.endswith('.fa'):
                tempFasta = os.path.join(root, f)
                options.input1 = root
                options.temp_fasta = tempFasta


def get_targetGenomeInfo(options, file_type):

    fp = open('./%s/targetGenome_locusTag_aaSeq.fa' %options.outputfolder2,'w')

    targetGenome_locusTag_aaSeq_dict = {}
    targetGenome_locusTag_ec_dict = {}
    targetGenome_locusTag_prod_dict = {}

    #Read GenBank file
    try:
        seq_record = SeqIO.read(options.input, file_type)
    except ValueError:
        logging.debug("Warning: ValueError occurred in SeqIo.read")
        seq_record = SeqIO.parse(options.outputfolder+'/'+options.input,
                     file_type).next()

    if options.eficaz:
        getECs(seq_record, options)

    total_cluster = 0
    locus_tag_list = []
    number_product_list = []
    number_ec_list = []

    for feature in seq_record.features:
        if feature.type == 'CDS':

            #Retrieving "locus_tag (i.e., ORF name)" for each CDS
            locusTag = feature.qualifiers['locus_tag'][0]
            locus_tag_list.append(locusTag)

            #Note that the numbers of CDS and "translation" do not match.
            #Some CDSs do not have "translation".
            if feature.qualifiers.get('translation'):
                translation = feature.qualifiers.get('translation')[0]
                targetGenome_locusTag_aaSeq_dict[locusTag] = translation
                print >>fp, '>%s\n%s' % (str(locusTag), str(translation))

            #Used to find "and" relationship in the GPR association
            if feature.qualifiers.get('product'):
                #It is confirmed that each locus_tag has a single '/product' annotation.
                #Thus, it's OK to use '[0]'.
                product = feature.qualifiers.get('product')[0]
                targetGenome_locusTag_prod_dict[locusTag] = product

            #Multiple 'EC_number's may exit for a single CDS.
            #Nver use '[0]' for the 'qualifiers.get' list.
            if feature.qualifiers.get('EC_number'):
                ecnum = feature.qualifiers.get('EC_number')
                targetGenome_locusTag_ec_dict[locusTag] = ecnum

        if feature.type == 'cluster':
            total_cluster += 1

    #Number of 'locus_tag's obtained above may be different from
    #the number directly obtained from genbank file
    #because some 'locus_tag's exit in the features 'tRNA' and 'rRNA',
    #which are not considered herein.
    logging.debug("Number of 'locus_tag's: %s" %len(locus_tag_list))

    #Same above comment for 'product's
    for locus_tag in targetGenome_locusTag_prod_dict.keys():
        number_product_list.append(targetGenome_locusTag_prod_dict[locus_tag])
    logging.debug("Number of 'product's: %s" %len(number_product_list))

    for locus_tag in targetGenome_locusTag_ec_dict.keys():
        for ec in targetGenome_locusTag_ec_dict[locus_tag]:
            number_ec_list.append(ec)
    logging.debug("Number of 'EC_number's: %s" %len(number_ec_list))

    logging.debug(
            "len(targetGenome_locusTag_prod_dict.keys):%s"
            %len(targetGenome_locusTag_prod_dict.keys()))

    logging.debug(
            "len(targetGenome_locusTag_ec_dict.keys): %s"
            %len(targetGenome_locusTag_ec_dict.keys()))

    options.targetGenome_locusTag_ec_dict = targetGenome_locusTag_ec_dict
    options.targetGenome_locusTag_prod_dict = targetGenome_locusTag_prod_dict
    options.total_cluster = total_cluster
    options.seq_record = seq_record

    fp.close()


def get_target_fasta(options):
    for root, _, files in os.walk('./%s' %options.outputfolder2):
        for f in files:
            if f.endswith('.fa'):
                target_fasta = os.path.join(root, f)
                options.target_fasta = target_fasta

    if not options.target_fasta:
        logging.warning("FASTA file %s for bidirectional blastp hits not found.")

