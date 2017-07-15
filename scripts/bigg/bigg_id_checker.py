#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cobra
import copy
import glob
import logging
import os
import pickle
import shutil
import sys
import zipfile
from cobra import Model, Reaction, Metabolite
from os.path import join, abspath, dirname

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import gems


input2_dir = join(os.pardir, 'gems', 'io', 'data', 'input2')
input2_tmp_dir = join(dirname(abspath(__file__)), 'input2_data')


# mnxm_compoundInfo_dict =
#{'MNXM128019': ['Methyl trans-p-methoxycinnamate', 'C11H12O3']}
def read_chem_prop(self, filename):
    mnxm_compoundInfo_dict = {}

    f = open(filename,'r')
    f.readline()

    for line in f:
        try:
            metab_prop_list = line.split('\t')
            mnxm_id = metab_prop_list[0].strip()
            mnxm_name = metab_prop_list[1].strip()
            mnxm_formula = metab_prop_list[2].strip()
            mnxm_compoundInfo_dict[mnxm_id] = [mnxm_name]
            mnxm_compoundInfo_dict[mnxm_id].append(mnxm_formula)
        except:
            logging.debug('Cannot parse MNXM: %s' %line)

    f.close()
    self.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

    return mnxm_compoundInfo_dict


def unzip_files():
    tsv_files = glob.glob(join(input2_tmp_dir, '*.tsv'))
    if len(tsv_files) ==  0:
        zip = zipfile.ZipFile(join(input2_tmp_dir, 'mnxref.zip'))
        zip.extractall(input2_tmp_dir)
        zip.close()


def create_zip_file():
    input2_tmp_dir_list = glob.glob(join(input2_tmp_dir, '*.*'))
    input2_tmp_dir_list2 = []

    for output in input2_tmp_dir_list:
        if '.tsv' not in output and \
                '.zip' not in output and \
                'bigg_old_new_dict.p' not in output:
            input2_tmp_dir_list2.append(output)

    zip = zipfile.ZipFile(join(input2_tmp_dir, 'mnxref_input2_data.zip'),
                            'w',
                            zipfile.ZIP_DEFLATED)

    for output in input2_tmp_dir_list2:
        zip.write(output, os.path.basename(output))

    zip.close()
    return input2_tmp_dir_list2


def remove_tsv_files(input2_tmp_dir_list2):
    tsv_files = glob.glob(join(input2_tmp_dir, '*.tsv'))
    for tsv_file in tsv_files:
        os.remove(tsv_file)

    for output in input2_tmp_dir_list2:
        os.remove(join(input2_tmp_dir, output))


if __name__ == '__main__':
    import time

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

    unzip_tsv_files()
    run_ParseMNXref()
    input2_tmp_dir_list2 = create_zip_file()
    remove_tsv_files(input2_tmp_dir_list2)

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
