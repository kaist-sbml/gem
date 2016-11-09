# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
import logging
import os
import sys
from os import path
import subprocess
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings
# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO
import re

class Storage(dict):
    """Simple storage class"""
    def __init__(self, indict=None):
        if indict is None:
            indict = {}
        dict.__init__(self, indict)
        self.__initialized = True

    def __getattr__(self, attr):
        try:
            return self.__getitem__(attr)
        except KeyError:
            raise AttributeError(attr)

    def __setattr__(self, attr, value):
        if '_Storage__initialized' not in self.__dict__:
            return dict.__setattr__(self, attr, value)
        elif attr in self:
            dict.__setattr__(self, attr, value)
        else:
            self.__setitem__(attr, value)


def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features

def get_all_features_of_type_with_query(seq_record, feature_type, query_tag, query_value):
    """Return all features of type 'type' which contain a 'query_tag' with value 'query_value'
    Note: query has to be exact!"""

    features_with_type = get_all_features_of_type(seq_record, feature_type)
    features = []
    for feature_to_test in features_with_type:

        if (query_tag in feature_to_test.qualifiers) and (query_value in feature_to_test.qualifiers[query_tag]):
            features.append(feature_to_test)
    return features

def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")

def get_withincluster_cds_features(seq_record):
    features = get_cds_features(seq_record)
    clusters = get_cluster_features(seq_record)
    withinclusterfeatures = []
    for feature in features:
        for cluster in clusters:
            if not (cluster.location.start <= feature.location.start <= cluster.location.end or \
               cluster.location.start <= feature.location.end <= cluster.location.end):
                continue
            if feature not in withinclusterfeatures:
                withinclusterfeatures.append(feature)
    return withinclusterfeatures

def get_cluster_cds_features(cluster, seq_record):
    clustercdsfeatures = []
    for feature in seq_record.features:
        if feature.type != 'CDS':
            continue
        if cluster.location.start <= feature.location.start <= cluster.location.end or \
           cluster.location.start <= feature.location.end <= cluster.location.end:
            clustercdsfeatures.append(feature)
    return clustercdsfeatures

def get_cluster_aSDomain_features(cluster, seq_record):
    aSDomainfeatures = []
    for feature in seq_record.features:
        if feature.type != 'aSDomain':
            continue
        if cluster.location.start <= feature.location.start <= cluster.location.end or \
           cluster.location.start <= feature.location.end <= cluster.location.end:
            aSDomainfeatures.append(feature)
    return aSDomainfeatures

def features_overlap(a, b):
    "Check if two features have overlapping locations"
    astart = min([a.location.start, a.location.end])
    aend = max([a.location.start, a.location.end])
    bstart = min([b.location.start, b.location.end])
    bend = max([b.location.start, b.location.end])
    return (astart >= bstart and astart <= bend) or \
           (aend >= bstart and aend <= bend) or \
           (bstart >= astart and bstart <= aend) or \
           (bend >= astart and bend <= aend)

def get_secmet_cds_features(seq_record):
    features = get_cds_features(seq_record)
    secmet_features = []
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            for annotation in feature.qualifiers['sec_met']:
                if annotation.startswith('Type: '):
                    secmet_features.append(feature)
                    break
    return secmet_features

def get_pksnrps_cds_features(seq_record):
    features = get_cds_features(seq_record)
    pksnrpscoregenes = []
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            for annotation in feature.qualifiers['sec_met']:
                if annotation.startswith('NRPS/PKS Domain:'):
                    pksnrpscoregenes.append(feature)
                    break
    return pksnrpscoregenes

def get_nrpspks_domain_dict(seq_record):
    domaindict = {}
    features = get_cds_features(seq_record)
    for feature in features:
        domainlist = []
        if 'sec_met' in feature.qualifiers:
            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
            for domain in domains:
                hit_id =  domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[0]
                domstart = domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[2].partition("-")[0]
                domend = domain.partition("NRPS/PKS Domain: ")[2].partition("). ")[0].rpartition("-")[2]
                evalue = domain.partition("E-value: ")[2].partition(". Score:")[0]
                bitscore = domain.partition("Score: ")[2].partition(";")[0]
                domainlist.append([hit_id, int(domstart), int(domend), evalue, float(bitscore)])
            if len(domainlist) > 0:
                domaindict[get_gene_id(feature)] = domainlist
    return domaindict

def get_nrpspks_substr_spec_preds(seq_record):
    substr_spec_preds = Storage()
    substr_spec_preds.consensuspreds = {}
    substr_spec_preds.nrps_svm_preds = {}
    substr_spec_preds.nrps_code_preds = {}
    substr_spec_preds.minowa_nrps_preds = {}
    substr_spec_preds.pks_code_preds = {}
    substr_spec_preds.minowa_pks_preds = {}
    substr_spec_preds.minowa_cal_preds = {}
    substr_spec_preds.kr_activity_preds = {}
    substr_spec_preds.kr_stereo_preds = {}
    features = get_cds_features(seq_record)
    for feature in features:
        nrat, nra, nrcal, nrkr = 0, 0, 0, 0
        if 'sec_met' in feature.qualifiers:
            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
            for domain in domains:
                if "AMP-binding" in domain or "A-OX" in domain:
                    nra += 1
                    domainname = get_gene_id(feature) + "_A" + str(nra)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    nrps_svm_pred = predictionstext.partition(" (NRPSPredictor2 SVM)")[0]
                    nrps_code_pred = predictionstext.partition(" (NRPSPredictor2 SVM), ")[2].partition(" (Stachelhaus code)")[0]
                    minowa_nrps_pred = predictionstext.partition("(Stachelhaus code), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.nrps_svm_preds[domainname] = nrps_svm_pred
                    substr_spec_preds.nrps_code_preds[domainname] = nrps_code_pred
                    substr_spec_preds.minowa_nrps_preds[domainname] = minowa_nrps_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_AT" in domain:
                    nrat += 1
                    domainname = get_gene_id(feature) + "_AT" + str(nrat)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    pks_code_pred = predictionstext.partition(" (PKS signature)")[0]
                    minowa_pks_pred = predictionstext.partition("(PKS signature), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.pks_code_preds[domainname] = pks_code_pred
                    substr_spec_preds.minowa_pks_preds[domainname] = minowa_pks_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "CAL_domain" in domain:
                    nrcal += 1
                    domainname = get_gene_id(feature) + "_CAL" + str(nrcal)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    minowa_cal_pred = predictionstext.partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.minowa_cal_preds[domainname] = minowa_cal_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_KR" in domain:
                    nrkr += 1
                    domainname = get_gene_id(feature) + "_KR" + str(nrkr)
                    activityprediction = domain.partition("Predicted KR activity: ")[2].partition(";")[0]
                    stereoprediction = domain.partition("Predicted KR stereochemistry: ")[2].partition(";")[0]
                    substr_spec_preds.kr_activity_preds[domainname] = activityprediction
                    substr_spec_preds.kr_stereo_preds[domainname] = stereoprediction
    return substr_spec_preds

def get_smcog_annotations(seq_record):
    smcogdict = {}
    smcogdescriptions = {}
    features = get_cds_features(seq_record)
    for feature in features:
        if 'note' in feature.qualifiers:
            notes = feature.qualifiers['note']
            for note in notes:
                if "smCOG: " in note:
                    smcogid = note.partition("smCOG: ")[2].partition(":")[0]
                    smcog_descr = note.partition("smCOG: ")[2].partition(":")[2].partition("(Score:")[0]
                    smcogdict[get_gene_id(feature)] = smcogid
                    smcogdescriptions[smcogid] = smcog_descr
    return smcogdict, smcogdescriptions


def get_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    return get_all_features_of_type(seq_record, "cluster")

def get_sorted_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    clusters = get_all_features_of_type(seq_record, "cluster")
    if len(clusters) == 0:
        return []
    numberdict = {}
    for cluster in clusters:
        numberdict[get_cluster_number(cluster)] = cluster
    return [numberdict[clusternr] for clusternr in numberdict.keys()]


def get_cluster_number(cluster):
    "Get the integer representation of the Cluster number qualifier"
    if 'note' in cluster.qualifiers and "Cluster number: " in cluster.qualifiers['note'][0]:
        clusternr = int(cluster.qualifiers['note'][0].partition("Cluster number: ")[2])
        return clusternr
    else:
        return 0

def get_cluster_type(cluster):
    "Get product type of a gene cluster"
    return cluster.qualifiers['product'][0]

def get_cluster_by_nr(seq_record, queryclusternr):
    "Return all cluster features for a seq_record of a certain type"
    clusters = get_all_features_of_type(seq_record, "cluster")
    for cluster in clusters:
        if "Cluster number: " in cluster.qualifiers['note'][0]:
            clusternr = int(cluster.qualifiers['note'][0].partition("Cluster number: ")[2])
            if clusternr == queryclusternr:
                return cluster

def get_cluster_features_of_type(seq_record, clustertype):
    "Return all cluster features for a seq_record of a certain type"
    clusters = get_all_features_of_type(seq_record, "cluster")
    return [cluster for cluster in clusters if clustertype in cluster.qualifiers['product'][0]]

def locate_executable(name):
    "Find an executable in the path and return the full path"
    # In windows, executables tend to end on .exe
    if sys.platform == 'win32':
        name += ".exe"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.X_OK):
            logging.debug("Found executable %r", name)
            return name
    for p in os.environ["PATH"].split(os.pathsep):
        full_name = path.join(p, name)
        if path.isfile(full_name) and os.access(full_name, os.X_OK):
            logging.debug("Found executable %r", full_name)
            return full_name

    return None

def locate_file(name):
    "Find a file and return the full path"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.R_OK):
            logging.debug("Found file %r", name)
            return name
    return None

# Ignore the pylint warning about input being redifined, as we're just
# following the subprocess names here.
# pylint: disable=redefined-builtin
def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError, e:
        logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, e)
        raise
# pylint: enable=redefined-builtin

def cmp_feature_location(a, b):
    "Compare two features by their start/end locations"
    ret = cmp(a.location.start, b.location.start)
    if ret != 0:
        return ret
    return cmp(a.location.end, b.location.end)

def sort_features(seq_record):
    "Sort features in a seq_record by their position"
    #Check if all features have a proper location assigned
    for feature in seq_record.features:
        if feature.location is None:
            if feature.id != "<unknown id>":
                logging.error("Feature '%s' has no proper location assigned", feature.id)
            elif "locus_tag" in feature.qualifiers:
                logging.error("Feature '%s' has no proper location assigned", feature.qualifiers["locus_tag"][0])
            else:
                logging.error("File contains feature without proper location assignment")
            sys.exit(0) #FIXME: is sys.exit(0) really what we want to do here?
    #Sort features by location
    seq_record.features.sort(cmp=cmp_feature_location)

def fix_locus_tags(seq_record, config):
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    if 'next_locus_tag' not in config:
        config.next_locus_tag = 1

    cds_list = get_cds_features(seq_record)
    for feature in cds_list:
        if get_gene_id(feature) == "no_tag_found":
            feature.qualifiers['locus_tag'] = ['AUTOORF_%05d' % config.next_locus_tag]
            config.next_locus_tag += 1
        #Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars  = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        if 'locus_tag' in feature.qualifiers:
            for char in feature.qualifiers['locus_tag'][0]:
                if char in illegal_chars:
                    feature.qualifiers['locus_tag'][0] = feature.qualifiers['locus_tag'][0].replace(char, "_")
        if 'gene' in feature.qualifiers:
            for char in feature.qualifiers['gene'][0]:
                if char in illegal_chars:
                    feature.qualifiers['gene'][0] = feature.qualifiers['gene'][0].replace(char, "_")
        if 'protein_id' in feature.qualifiers:
            for char in feature.qualifiers['protein_id'][0]:
                if char in illegal_chars:
                    feature.qualifiers['protein_id'][0] = feature.qualifiers['protein_id'][0].replace(char, "_")

def _shorten_ids(idstring, options):
    contigstrmatch = re.search(r"onti?g?(\d+)\b", idstring)
    if not contigstrmatch:
        # if there is a substring "[Ss]caf(fold)XXX" use this number
        contigstrmatch = re.search(r"caff?o?l?d?(\d+)\b", idstring)
    if not contigstrmatch:
        # if there is a substring " cXXX" use this number
        contigstrmatch = re.search(r"\bc(\d+)\b", idstring)
    if contigstrmatch:
        contig_no = int(contigstrmatch.group(1))
    else:
        # if the contig number cannot be parsed out, just count the contigs from 1 to n
        contig_no = options.orig_record_idx

    return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])



def get_gene_acc(feature):
    "Get the gene accesion number; if not defined, use content of locus_tag or gene qualifier instead"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]

    # I needed to include this for clusterblast to work as well with the data from db (which all has protein_id) as with user data (which hasn't)
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return "no_accession"

def get_gene_accs(feature):
    "Get array of the protein_id tags"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id']
    elif 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag']
    elif 'gene' in feature.qualifiers:
        return feature.qualifiers['gene']
    return ["no_accession"]

def get_ncbi_gi(feature):
    """Get the NCBI gi from feature
    returns gi if found, None if not present"""

    gi = None
    if 'db_xref' in feature.qualifiers:
        for db_xref in feature.qualifiers['db_xref']:
            if 'GI:' in db_xref:
                gi = db_xref[3:]
    return gi

def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"

def get_gene_annotation(feature):
    "Get the gene annotation from the product qualifier"
    if 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    return "unannotated orf"

def get_gene_accession(feature):
    "Get the gene accession number from protein_id, locus_tag or gene name, in that order"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return "no_tag_found"

def get_locus_tag(feature):
    "Get the locus tag of a feature, no fallback"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    else:
        logging.exception("Feature found without locus tag! Feature start: %s end: %s", str(feature.location.start), str(feature.location.end))
        raise ValueError

def get_full_path(current_file, file_to_add):
    "Get the full path of file_to_add in the same directory as current_file"
    return path.join(path.dirname(path.abspath(current_file)), file_to_add)

def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_id = get_gene_id(feature)
        feature_by_id[gene_id] = feature
    return feature_by_id


# Should be refactored!!
# TODO: Add tests

def get_feature_dict_locus_tag(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        locus_tag = get_locus_tag(feature)
        feature_by_id[locus_tag] = feature
    return feature_by_id

def get_feature_dict_protein_id(seq_record):
    """Get a dictionary mapping features to their protein_ids"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_ids = get_gene_accs(feature)

        # As the HMMer result file does only contain result lines for 1 protein_id;
        # if an entry has 2 (or more) protein_id entries, lookups failed in the original version
        # Therefore I changed this method to use the protein_id list instead of just
        # extracting the [0] element - althouth we loose some memory here due to redundancy.
        for gene_id in gene_ids:
            feature_by_id[gene_id] = feature
    return feature_by_id



def get_aa_sequence(feature, to_stop=False):
    """Extract sequence from specific CDS feature in sequence record"""
    fasta_seq = feature.qualifiers['translation'][0]
    if "*" in fasta_seq:
        if to_stop:
            fasta_seq = fasta_seq.split('*')[0]
        else:
            fasta_seq = fasta_seq.replace("*","X")
    if "-" in fasta_seq:
        fasta_seq = fasta_seq.replace("-","")
    return fasta_seq

def writefasta(names, seqs, filename):
    "Write sequence to a file"
    e = 0
    f = len(names) - 1
    out_file = open(filename,"w")
    while e <= f:
        out_file.write(">")
        out_file.write(names[e])
        out_file.write("\n")
        out_file.write(seqs[e])
        out_file.write("\n")
        e += 1
    out_file.close()

def sortdictkeysbyvaluesrev(indict):
    items = [(value, key) for key, value in indict.items()]
    items.sort()
    items.reverse()
    return [key for value, key in items]

def sortdictkeysbyvaluesrevv(indict):
    values = indict.values()
    values.sort()
    values.reverse()
    return values

def get_genefinding_basedir(options):
    "Identify the basedir for glimmer files"
    basedir = ''
    if 'glimmer' in options:
        if 'basedir' in options.glimmer:
            basedir = options.glimmer.basedir
    return basedir

def zip_dir(dir_path, archive, prefix_to_remove=""):
    """Recursively add a directory's contents to a zip archive"""
    entries = os.listdir(dir_path)
    for entry in entries:
        entry = path.join(dir_path, entry)
        if path.isdir(entry):
            zip_dir(entry, archive, prefix_to_remove)
        else:
            arcname = entry.replace(prefix_to_remove + os.sep, "", 1)
            archive.write(entry, arcname)



def get_git_version():
    """Get the sha1 of the current git version"""
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        out, _, retcode = execute(args)
        return out.strip()
    except OSError:
        pass

    return ""

