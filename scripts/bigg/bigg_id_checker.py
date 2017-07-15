#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import cobra
import glob
import logging
import os
import zipfile
from os.path import join, abspath, dirname

bigg_dir = join(dirname(abspath(__file__)))

def unzip_files():
    zip = zipfile.ZipFile(join(bigg_dir, 'bigg.zip'))
    zip.extractall(bigg_dir)
    zip.close()


def get_txt_files():
   file_list = glob.glob(join(bigg_dir, '*.txt'))
   logging.debug("BiGG txt files found: %s", file_list)
   return file_list


def parse_bigg_metab_file(file_list):
    bigg_metab_list = []

    for filename in file_list:
        if 'metabolites' in filename:
            f = open(filename,'r')
            f.readline()

    for line in f:
        try:
            metab_id_list = line.split('\t')
            bigg_metab_id = metab_id_list[0].strip()
            bigg_metab_list.append(bigg_metab_id)
            logging.debug("%s added to 'bigg_metab_list'", bigg_metab_id)
        except:
            logging.debug('Cannot parse: %s' %line)

    f.close()
    return bigg_metab_list


def parse_bigg_rxn_file(file_list):
    bigg_rxn_list = []

    for filename in file_list:
        if 'reactions' in filename:
            f = open(filename,'r')
            f.readline()

    for line in f:
        try:
            rxn_id_list = line.split('\t')
            bigg_rxn_id = rxn_id_list[0].strip()
            bigg_rxn_list.append(bigg_rxn_id)
            logging.debug("%s added to 'bigg_rxn_list'", bigg_rxn_id)
        except:
            logging.debug('Cannot parse: %s' %line)

    f.close()
    return bigg_rxn_list


def check_bigg_id_consistency(bigg_metab_list, bigg_rxn_list):
    model = cobra.io.read_sbml_model('model.xml')

    metab_cnt = 0
    for i in range(len(model.metabolites)):
        if model.metabolites[i].id in bigg_metab_list:
            metab_cnt += 1

    rxn_cnt = 0
    for i in range(len(model.reactions)):
        if model.reactions[i].id in bigg_rxn_list:
            rxn_cnt += 1

    consist_metab_ratio = metab_cnt/len(model.metabolites)
    consist_rxn_ratio = rxn_cnt/len(model.reactions)

    logging.debug("Ratio of metabolites with consistent BiGG IDs: %s (%s / %s)",
                    consist_metab_ratio, metab_cnt, len(model.metabolites))

    logging.debug("Ratio of reactions with consistent BiGG IDs: %s (%s / %s)",
                    consist_rxn_ratio, rxn_cnt, len(model.reactions))


def remove_files(file_list):

    for filename in file_list:
        os.remove(filename)


if __name__ == '__main__':
    import time

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

    # Logfile setup
    logger = logging.getLogger('')
    fomatter = logging.Formatter('[%(levelname)s|%(filename)s:%(lineno)s] > %(message)s')
    fh = logging.FileHandler(join(bigg_dir, 'bigg_id_checker.log'), mode = 'w')
    fh.setFormatter(fomatter)
    logger.addHandler(fh)

    unzip_files()
    file_list = get_txt_files()
    bigg_metab_list = parse_bigg_metab_file(file_list)
    bigg_rxn_list = parse_bigg_rxn_file(file_list)
    check_bigg_id_consistency(bigg_metab_list, bigg_rxn_list)
    remove_files(file_list)
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
