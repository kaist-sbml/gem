
import logging
import os

def get_features_from_gbk(seq_record, options):

    for feature in seq_record.features:
        if feature.type == 'CDS':

            #Retrieving "locus_tag (i.e., ORF name)" for each CDS
            locusTag = feature.qualifiers['locus_tag'][0]

            #Note that the numbers of CDS and "translation" do not match.
            #Some CDSs do not have "translation".
            if feature.qualifiers.get('translation'):
                translation = feature.qualifiers.get('translation')[0]
                options.targetGenome_locusTag_aaSeq_dict[locusTag] = translation

            #Used to find "and" relationship in the GPR association
            if feature.qualifiers.get('product'):
                #It is confirmed that each locus_tag has a single '/product' annotation.
                #Thus, it's OK to use '[0]'.
                product = feature.qualifiers.get('product')[0]
                options.targetGenome_locusTag_prod_dict[locusTag] = product

            #Multiple 'EC_number's may exit for a single CDS.
            #Nver use '[0]' for the 'qualifiers.get' list.
            if feature.qualifiers.get('EC_number'):
                ecnum = feature.qualifiers.get('EC_number')
                options.targetGenome_locusTag_ec_dict[locusTag] = ecnum

        if feature.type == 'cluster':
            options.total_cluster += 1


def get_features_from_fasta(seq_record, options):
    locusTag = seq_record.id
    options.targetGenome_locusTag_aaSeq_dict[locusTag] = seq_record.seq
    options.targetGenome_locusTag_prod_dict[locusTag] = seq_record.description


def get_target_fasta(options):

    if options.targetGenome_locusTag_aaSeq_dict:

        target_fasta_dir = os.path.join(
                options.outputfolder2, 'targetGenome_locusTag_aaSeq.fa')
        with open(target_fasta_dir,'w') as f:
            for locusTag in options.targetGenome_locusTag_aaSeq_dict.keys():
                print >>f, '>%s\n%s' \
                %(str(locusTag), str(options.targetGenome_locusTag_aaSeq_dict[locusTag]))
        options.target_fasta = target_fasta_dir
    else:
        logging.warning("FASTA file 'targetGenome_locusTag_aaSeq.fa' not found")


#Look for pre-stored fasta file of the template model
def get_temp_fasta(options):
    for root, _, files in os.walk('./gmsm/io/data/input1/%s/' %(options.orgName)):
        for f in files:
            if f.endswith('.fa'):
                tempFasta = os.path.join(root, f)
                options.input1 = root
                options.temp_fasta = tempFasta

    if options.temp_fasta:
        logging.debug("FASTA file for '%s' found", options.orgName)
    else:
        logging.warning("FASTA file for '%s' not found", options.orgName)
