#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ast
import cobra
#import copy
#import glob
import logging
import os
#import pickle
#import shutil
import sys
import urllib2
#import zipfile
from cobra import Model, Reaction, Metabolite
from os.path import join, abspath, dirname

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import gmsm


input2_dir = join(os.pardir, 'gmsm', 'io', 'data', 'input2')
input2_tmp_dir = join(dirname(abspath(__file__)), 'input2_data')


class CreateBiGGModel(object):

    def get_bigg_rxns(self):

        url = 'http://bigg.ucsd.edu/api/v2/universal/reactions'
        logging.debug('URL for downloading a full list of reactions from BiGG Models:')
        logging.debug(url)

        bigg_rxns = urllib2.urlopen(url).read()

        bigg_rxns_fixed = bigg_rxns.replace('null', '""')

        with open(join(input2_tmp_dir, 'bigg_rxns.txt'), 'wb') as f:
            f.write(bigg_rxns_fixed)

        bigg_rxns_dict = ast.literal_eval(bigg_rxns_fixed)

        self.bigg_rxns_dict = bigg_rxns_dict


    def get_rxn_info(self):

        for rxn in self.bigg_rxns_dict['results']:
            rxn_info = ''
            rxn_info_dict = {}

            bigg_id = rxn['bigg_id']
            url = ''.join(['http://bigg.ucsd.edu/api/v2/universal/reactions/', bigg_id])

            rxn_info = urllib2.urlopen(url).read()

            rxn_info_fixed = rxn_info.replace('null', '""')
            rxn_info_fixed = rxn_info.replace('false', 'False')

            rxn_info_dict = ast.literal_eval(rxn_info_fixed)
            logging.debug('%s: %s', bigg_id, rxn_info_dict['metabolites'])


    #TODO: Should be modified, may be combined with the function above
    def get_bigg_metabolites():
        url = 'http://bigg.ucsd.edu/api/v2/universal/metabolites'
        logging.debug('URL for downloading a model from the BiGG Models:')
        logging.debug(url)

        bigg_metabs = urllib2.urlopen(url).read()

        with open(join(input2_tmp_dir, 'bigg_metabs.txt'), 'wb') as f:
            f.write(bigg_metabs)

    #TODO: Should be fixed
    def read_reac_prop(self, filename):
        reaction_info = {}

        sptlist = line.split('\t')
        reaction_id = sptlist[0].strip()
        reversibility = True

        stoich_dict = {}
        equation = sptlist[1].strip()

        reactant_info, product_info = self.parse_equation(equation)

        for metabolite in reactant_info:
            stoich_dict[metabolite] = reactant_info[metabolite]
        for metabolite in product_info:
            stoich_dict[metabolite] = product_info[metabolite]

        flag = True
        for each_reactant in reactant_info.keys():
            if each_reactant in product_info.keys():
                flag = False
                break

        if flag == False:
            continue

        ec_number_list = ec_number.split(';')

        if len(ec_number_list) > 0:
            reaction_info[reaction_id] = {
                            'stoichiometry': stoich_dict,
                            'reversibility': reversibility}

        self.reaction_info = reaction_info


    def get_cobra_reactions(self):
        logging.debug('Creating MNXR reactions: time-consuming')

        cobra_reactions = []
        cobra_metabolites = []

        cnt = 0
        for each_reaction in self.reaction_info:
            cnt += 1
            logging.debug('Total reaction number %s; Reaction number covered %s; %s' \
                    %(len(self.reaction_info), cnt, each_reaction))
            reaction_name = each_reaction
            metabolites = self.reaction_info[each_reaction]['stoichiometry']
            reaction_reversibility = self.reaction_info[each_reaction]['reversibility']

            compartment_list = ['c']
            compartment_list = list(set(compartment_list))
            flag = True

            for each_compartment in compartment_list:
                if each_compartment == 'e':
                    continue

                new_reaction_metabolite_obj = {}
                for each_metabolite in metabolites:
                    # Convert MNXM to BiGG IDs
                    if each_metabolite in self.mnxm_bigg_compound_dict.keys() \
                            and each_metabolite in self.mnxm_compoundInfo_dict.keys():
                        metabolite_id = '%s_%s' \
                            %(self.mnxm_bigg_compound_dict[each_metabolite.strip()],
                                each_compartment)
                        metabolite_obj = Metabolite(
                            str(metabolite_id),
                            name = self.mnxm_compoundInfo_dict[each_metabolite][0],
                            formula = self.mnxm_compoundInfo_dict[each_metabolite][1],
                            compartment = str(each_compartment))
                    elif each_metabolite not in self.mnxm_bigg_compound_dict.keys() \
                            and each_metabolite in self.mnxm_compoundInfo_dict.keys():
                        metabolite_id = '%s_%s' \
                            %(each_metabolite.strip(), each_compartment)
                        metabolite_obj = Metabolite(
                            str(metabolite_id),
                            name = self.mnxm_compoundInfo_dict[each_metabolite][0],
                            formula = self.mnxm_compoundInfo_dict[each_metabolite][1],
                            compartment = str(each_compartment))
                    elif each_metabolite not in self.mnxm_bigg_compound_dict.keys() \
                            and each_metabolite not in self.mnxm_compoundInfo_dict.keys():
                        metabolite_id = '%s_%s' \
                            %(each_metabolite.strip(), each_compartment)
                        metabolite_obj = Metabolite(
                            str(metabolite_id),
                            name = '',
                            formula = '',
                            compartment = str(each_compartment))

                    # Insert other db IDs
                    metabolite_obj.notes = {}
                    metabolite_obj.notes['MNXM'] = each_metabolite
                    if each_metabolite in self.mnxm_bigg_compound_dict:
                        metabolite_obj.notes['BiGG'] = \
                                self.mnxm_bigg_compound_dict[each_metabolite]
                    if each_metabolite in self.mnxm_kegg_compound_dict:
                        metabolite_obj.notes['KEGG'] = \
                                self.mnxm_kegg_compound_dict[each_metabolite]

                    coeff = metabolites[each_metabolite]
                    new_reaction_metabolite_obj[metabolite_obj] = float(coeff)

                if reaction_reversibility == True:
                    lb = -1000.0
                    ub = 1000.0
                    rev = 1
                else:
                    lb = 0.0
                    ub = 1000.0
                    rev = 0

                if flag == True:
                    new_reaction_id = str(each_reaction)
                    reaction_obj = Reaction(new_reaction_id)

                    # Insert bigg or kegg reaction IDs
                    try:
                        reaction_obj.name = ';'.join(self.mnxr_xref_dict[each_reaction])
                    except:
                        reaction_obj.name = ''

                    reaction_obj.subsystem = ''
                    reaction_obj.lower_bound = lb
                    reaction_obj.upper_bound = ub
                    reaction_obj.objective_coefficient = 0
                    reaction_obj.reversibility = rev
                    reaction_obj.gene_reaction_rule = ''
                    reaction_obj.add_metabolites(new_reaction_metabolite_obj)

                    # Currently writing to sbml not supported
                    reaction_obj.notes = {}
                    reaction_obj.notes['EC_number'] = ';'.join(ec_number_list)
                    reaction_obj.notes['Balance'] = mass_balance
                    cobra_reactions.append(copy.deepcopy(reaction_obj))

        self.cobra_reactions = cobra_reactions


    def make_cobra_model(self, mnx_file):
        cobra_model = Model('mnxref_model')
        self.read_reac_prop(mnx_file)
        self.get_cobra_reactions()

        reaction_list = []
        for reaction_id in cobra_model.reactions:
            reaction_list.append(reaction_id.id)

        for each_cobra_reaction in self.cobra_reactions:
            if each_cobra_reaction.id not in reaction_list:
                cobra_model.add_reaction(each_cobra_reaction)

        for i in range(len(cobra_model.metabolites)):
            metab = cobra_model.metabolites[i]
            metab.id = cobra.io.sbml.fix_legacy_id(metab.id)

        cobra_model = gmsm.utils.stabilize_model(
                cobra_model, input2_tmp_dir, 'MNXref', diff_name=True)

        logging.debug('%i reactions in model' % len(cobra_model.reactions))
        logging.debug('%i metabolites in model' % len(cobra_model.metabolites))
        logging.debug('%i genes in model' % len(cobra_model.genes))

        return cobra_model


def unzip_bigg():
    zip = zipfile.ZipFile(join(input2_tmp_dir, 'bigg.zip'))
    zip.extractall(input2_tmp_dir)


def run_CreateBiGGModel():
    bigg_modeler = CreateBiGGModel()

    bigg_modeler.get_bigg_rxns()
    bigg_modeler.get_rxn_info()
#    mnxm_compoundInfo_dict = bigg_modeler.read_bigg_metabs(join(input2_tmp_dir, 'bigg_models_metabolites.txt'))
#    mnxr_kegg_dict, bigg_mnxr_dict = bigg_modeler.read_bigg_rxns(join(input2_tmp_dir, 'bigg_models_reactions.txt'))
#    cobra_model = bigg_modeler.make_cobra_model(join(input2_tmp_dir, 'reac_prop.tsv'))

    # Write SBML file
#    cobra.io.write_sbml_model(cobra_model,
#            join(input2_tmp_dir, 'MNXref.xml'), use_fbc_package=False)

    # Copy pickles to the destination
    # NOTE: Currently disabled as the bigg model is not used for gapfilling in GMSM
    #shutil.copyfile(join(input2_dir, 'MNXref.p'),
    #        join(os.pardir, 'gmsm', 'tests', 'data', 'MNXref.p'))


def remove_txt_files():
    tsv_files = glob.glob(join(input2_tmp_dir, '*.tsv'))
    for tsv_file in tsv_files:
        os.remove(tsv_file)


if __name__ == '__main__':
    import time

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

#    bigg_rxns_dict = get_bigg_reactions()
#    get_rxn_info(bigg_rxns_dict)
#    get_bigg_metabolites()
#    unzip_bigg()
    run_CreateBiGGModel()
#    remove_txt_files()

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
