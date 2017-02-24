# -*- coding: utf-8 -*-

import copy
import glob
import os
import pickle
import zipfile
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model
from os.path import join, abspath, dirname

mnxref_dir = join(dirname(abspath(__file__)), 'input2_data')

class ParseMNXM(object):

    #TODO: Try to import cobrapy function rather than explicitly writing here
    #Based on fix_legacy_id(id, use_hyphens=False, fix_compartments=False) of cobrapy:
    def replace_special_characters_compoundid(biggid):
        biggid = biggid.replace('-', '_DASH_')
        biggid = biggid.replace('/', '_FSLASH_')
        biggid = biggid.replace("\\",'_BSLASH_')
        biggid = biggid.replace('(', '_LPAREN_')
        biggid = biggid.replace('[', '_LSQBKT_')
        biggid = biggid.replace(']', '_RSQBKT_')
        biggid = biggid.replace(')', '_RPAREN_')
        biggid = biggid.replace(',', '_COMMA_')
        biggid = biggid.replace('.', '_PERIOD_')
        biggid = biggid.replace("'", '_APOS_')
        biggid = biggid.replace('&', '&amp;')
        biggid = biggid.replace('<', '&lt;')
        biggid = biggid.replace('>', '&gt;')
        biggid = biggid.replace('"', '&quot;')

        return biggid

    # This function cannot parse the initial version of NMXref data
    def read_metabolite_info(self, filename):
        bigg_mnxm_compound_dict = {}
        mnxm_bigg_compound_dict = {}
        kegg_mnxm_compound_dict = {}
        mnxm_kegg_compound_dict = {}

        f = open(filename,'r')
        f.readline()

        for line in f:
            try:
                metab_info_list = line.split('\t')
                xref = metab_info_list[0].strip()
                xref_list = xref.split(':')
                xref_db = xref_list[0].strip()
                xref_id = xref_list[1].strip()
                mnxm = metab_info_list[1].strip()

                if xref_db == 'bigg':
                    biggid = replace_special_characters_compoundid(xref_id)
                    bigg_mnxm_compound_dict[biggid] = mnxm
                    mnxm_bigg_compound_dict[mnxm] = biggid
                elif xref_db == 'kegg':
                    kegg_mnxm_compound_dict[xref_id] = mnxm

                    # Following conditions give priority to compoundID starting with 'C'
                    if 'D' not in xref_id \
                            and 'E' not in xref_id \
                            and 'G' not in xref_id \
                            and mnxm not in mnxm_kegg_compound_dict.keys():
                        mnxm_kegg_compound_dict[mnxm] = xref_id

                logging.debug('%s; %s; %s' %(xref_db, xref_id, mnxm))
            except:
                logging.debug('Cannot parse MNXM: %s' %line)

        f.close()
        return bigg_mnxm_compound_dict, mnxm_bigg_compound_dict,\
                kegg_mnxm_compound_dict, mnxm_kegg_compound_dict


class ParseMNXR(object):

    def parse_equation(self, equation):
        equation = equation.replace('>', '')
        equation = equation.replace('<', '')
        spt_equation = equation.split('=')
        reactant_str = spt_equation[0].strip()
        product_str = spt_equation[1].strip()

        reactant_info = {}
        product_info = {}
        metabolite_info = {}
        reactant_list = reactant_str.split('+')
        for each_metabolite in reactant_list:
            spt_metabolite = each_metabolite.strip().split(' ')

            if len(spt_metabolite) == 1:
                metabolite_name = spt_metabolite[0].strip()
                coeff = 1.0
                reactant_info[metabolite_name] = coeff * -1
            elif len(spt_metabolite) == 2:
                coeff = spt_metabolite[0].strip()
                metabolite_name = spt_metabolite[1].strip()
                reactant_info[metabolite_name] = float(coeff) * -1

        product_list = product_str.split('+')
        for each_metabolite in product_list:
            spt_metabolite = each_metabolite.strip().split(' ')

            if len(spt_metabolite) == 1:
                metabolite_name = spt_metabolite[0].strip()
                coeff = 1.0
                product_info[metabolite_name] = float(coeff)
            elif len(spt_metabolite) == 2:
                coeff = spt_metabolite[0].strip()
                metabolite_name = spt_metabolite[1].strip()
                product_info[metabolite_name] = float(coeff)

        return reactant_info, product_info

    def read_reaction_info(self, filename):
        reaction_info = {}
        mass_balance = ''
        ec_number = ''

        fp = open(filename, 'r')
        fp.readline()  # dummy line
        for line in fp:
            try:
                sptlist = line.split('\t')
                reaction_id = sptlist[0].strip()
                reversibility = True

                stoich_dict = {}
                equation = sptlist[1].strip()

                if sptlist[3].strip() == 'true':
                    mass_balance = 'balanced'
                elif sptlist[3].strip() == 'false':
                    mass_balance = 'unbalanced'

                ec_number = sptlist[4].strip()

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
                            'balance': mass_balance,
                            'ec': ec_number_list,
                            'reversibility': reversibility}
            except:
                logging.debug('Cannot parse MNXR: %s' %line)

        self.reaction_info = reaction_info
        fp.close()
        return

    def get_cobra_reactions(self):
        logging.debug('Creating MNXR reactions: time-consuming')

        cobra_reactions = []

        for each_reaction in self.reaction_info:
            reaction_name = each_reaction
            metabolites = self.reaction_info[each_reaction]['stoichiometry']
            mass_balance = self.reaction_info[each_reaction]['balance']
            ec_number_list = self.reaction_info[each_reaction]['ec']
            reaction_reversibility = self.reaction_info[each_reaction]['reversibility']

            compartment_list = ['c']
            compartment_list = list(set(compartment_list))
            flag = True

            for each_compartment in compartment_list:
                if each_compartment == 'e':
                    continue

                new_reaction_metabolite_obj = {}
                for each_metabolite in metabolites:
                    metabolite_id = '%s_%s' % (each_metabolite.strip(), each_compartment)
                    coeff = metabolites[each_metabolite]
                    new_metabolite_obj = Metabolite(
                            str(metabolite_id),
                            name=str(each_metabolite),
                            compartment=str(each_compartment))
                    new_reaction_metabolite_obj[new_metabolite_obj] = float(coeff)

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
                    new_reaction_id = new_reaction_id.replace(':', '_')
                    new_reaction_id = new_reaction_id + '_' + each_compartment
                    reaction_obj = Reaction(new_reaction_id)
                    reaction_obj.name = str(reaction_name)
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

        return cobra_reactions

    def make_cobra_model(self, mnx_file):
        cobra_model = Model('mnxref_model')
        self.read_reaction_info(mnx_file)
        cobra_reactions = self.get_cobra_reactions()

        reaction_list = []
        for reaction_id in cobra_model.reactions:
            reaction_list.append(reaction_id.id)

        for each_cobra_reaction in cobra_reactions:
            if each_cobra_reaction.id not in reaction_list:
                print each_cobra_reaction, each_cobra_reaction.notes
                cobra_model.add_reaction(each_cobra_reaction)

        logging.debug('%i reactions in model' % len(cobra_model.reactions))
        logging.debug('%i metabolites in model' % len(cobra_model.metabolites))
        logging.debug('%i genes in model' % len(cobra_model.genes))

        return cobra_model


def unzip_tsv_files():
    tsv_files = glob.glob(join(mnxref_dir, '*.tsv'))
    if len(tsv_files) ==  0:
        zip = zipfile.ZipFile(join(mnxref_dir, 'mnxref.zip'))
        zip.extractall(mnxref_dir)


def pickle_metabolite_id_conversions():
    mnxm_parser = ParseMNXM()
    mnxm_parser.read_metabolite_info(join(mnxref_dir, 'reac_prop.tsv'))

def convert_mnxref_to_sbml():
    mnxr_parser = ParseMNXR()
    cobra_model = mnxr_parser.make_cobra_model(join(mnxref_dir, 'reac_prop.tsv'))
    write_sbml_model(cobra_model,
        join(mnxref_dir, 'universal_model.xml'), use_fbc_package=False)
    return cobra_model


def pickle_sbml(cobra_model):
    with open(join(mnxref_dir, 'universal_model.p'), 'wb') as f:
        pickle.dump(cobra_model, f, protocol=pickle.HIGHEST_PROTOCOL)


def remove_tsv_files():
    tsv_files = glob.glob(join(mnxref_dir, '*.tsv'))
    for tsv_file in tsv_files:
        os.remove(tsv_file)


if __name__ == '__main__':
    import logging

    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')

    unzip_tsv_files()
    cobra_model = convert_mnxr_to_sbml()
    pickle_sbml(cobra_model)
    remove_tsv_files()
