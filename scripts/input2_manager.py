#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cobra
import copy
import glob
import logging
import os
import pandas as pd
import pickle
import shutil
import sys
import time
import zipfile
from cobra import Model, Reaction, Metabolite
from os.path import join, abspath, dirname

sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import gmsm


input2_dir = join(os.pardir, 'gmsm', 'io', 'data', 'input2')
input2_tmp_dir = join(dirname(abspath(__file__)), 'input2_data')

class ParseMNXref(object):

    # Based on King et al. (2016) in NAR
    def fix_legacy_id_using_BiGGModels(self):
        pickle_dir = join(input2_tmp_dir, 'bigg_old_new_dict.p')

        if os.path.isfile(pickle_dir):
            with open(pickle_dir, 'rb') as f:
                bigg_old_new_dict = pickle.load(f)

        elif not os.path.isfile(pickle_dir):
            bigg_old_new_dict = {}

            zip = zipfile.ZipFile(join(input2_tmp_dir, 'King-etal-2016_fix_legacy_id.zip'))
            zip.extractall(input2_tmp_dir)

            # File 1: 'King-etal-2016_fix_legacy_id_metabolites.txt',
            # File 2: 'King-etal-2016_fix_legacy_id_reactions.txt']
            file_list = glob.glob(join(input2_tmp_dir, '*.txt'))

            for filename in file_list:
                if 'legacy' in filename:
                    f = open(join(input2_tmp_dir, filename),'r')
                    f.readline()

                    for line in f:
                        try:
                            id_list = line.split('\t')
                            old_id = id_list[0].strip()
                            new_id = id_list[1].strip()
                            bigg_old_new_dict[old_id] = new_id
                            logging.debug('%s; %s' %(old_id, new_id))
                        except:
                            logging.debug('Cannot parse BiGG Models IDs: %s' %line)

                    f.close()
                    os.remove(join(input2_tmp_dir, filename))

            zip.close()
        return bigg_old_new_dict


    # This function cannot parse the initial version of NMXref data
    def read_chem_xref(self, bigg_old_new_dict, df):
        """ Process lines from cross reference file chem_xref.

        Output
        ------------------------------------------
        mnxm_bigg_compound_dict
        mnxm_kegg_compound_dict

        Parameters
        ------------------------------------------
        bigg_old_new_dict : dict
            dict to reference old and new BiGG IDs Based on King et al. (2016) in NAR
        df : pandas.DataFrame
            DataFrame object which contains chem_xref information

        # Note : Ignore the information from the databases other than bigg and kegg for metabolites
        """

        mnxm_bigg_compound_dict = {}
        mnxm_kegg_compound_dict = {}

        logging.debug('Parsing of chem_xref initiated...')

        for i in df.index:
            try:
                xref = df.loc[i, '#source']
                if xref.startswith(('bigg', 'kegg')):
                    xref_list = xref.split(':')
                    xref_db = xref_list[0]
                    xref_id = xref_list[1]
                    mnxm = df.loc[i, 'ID']

                    # xref_id with 'M_' is secondary/obsolete/fantasy identifier in MNXref
                    if xref_db == 'biggM' and not xref_id.startswith('M_'):
                        if xref_id in bigg_old_new_dict.keys():
                            mnxm_bigg_compound_dict[mnxm] = bigg_old_new_dict[xref_id]
                        else:
                            mnxm_bigg_compound_dict[mnxm] = xref_id
                        logging.debug('%s; %s; %s' %(xref_db, xref_id, mnxm))

                   # use only compound (keggC) identifiers, not other kegg identifiers (keggD, keggE, keggG)
                    elif xref_db == 'keggC' and not xref_id.startswith('M_'):
                        if mnxm not in mnxm_kegg_compound_dict.keys():
                            mnxm_kegg_compound_dict[mnxm] = xref_id
                        logging.debug('%s; %s; %s' %(xref_db, xref_id, mnxm))

            except Exception as e:
                logging.debug(e)
                logging.debug('Cannot parse MNXM: %s' %df.loc[i,:].values)

        logging.debug('Cross reference dictionary for MNXM compounds to bigg has %d compounds' % len(mnxm_bigg_compound_dict))
        logging.debug('Cross reference dictionary for MNXM compounds to kegg has %d compounds' % len(mnxm_kegg_compound_dict))

        logging.debug('Parsing of chem_xref completed')

        self.mnxm_bigg_compound_dict = mnxm_bigg_compound_dict
        self.mnxm_kegg_compound_dict = mnxm_kegg_compound_dict

        return mnxm_bigg_compound_dict


    def read_chem_prop(self, df):
        ''' Create metabolites information dictionary 
        chem_prop has following fields :
        #ID	name	reference	formula	charge	mass	InChI	InChIKey	SMILES
        At the moment we store information of #ID, name, formula and charge as required for GMSM

        Output
        ------------------------------------------
        mnxm_compoundInfo_dict
        {'MNXM128019': ['methyl 4-methoxycinnamate', 'C11H12O3', 0]}

        Parameters
        ------------------------------------------
        df: pandas.DataFrame
            DataFrame object which contains chem_prop information
        '''
        mnxm_compoundInfo_dict = {}

        for i in df.index:
            try:
                mnxm_id = df.loc[i, '#ID']
                mnxm_name = df.loc[i, 'name']
                mnxm_formula = df.loc[i, 'formula']
                mnxm_charge = df.loc[i, 'charge']
                if mnxm_charge:
                    mnxm_charge = int(mnxm_charge)
                else:
                    mnxm_charge = 0
                mnxm_compoundInfo_dict[mnxm_id] = [mnxm_name, mnxm_formula, mnxm_charge]

            except Exception as e:
                logging.debug(e)
                logging.debug('Cannot parse MNXM: %s' %df.loc[i,:].values)

        logging.debug('Metabolite info dictionary has %d compounds' % len(mnxm_compoundInfo_dict))

        logging.debug('Parsing of chem_prop completed')
        self.mnxm_compoundInfo_dict = mnxm_compoundInfo_dict

        return mnxm_compoundInfo_dict


    def read_reac_xref(self, df):
        """ Process lines from cross reference file reac_xref.

        Output
        ------------------------------------------
        mnxr_kegg_dict, bigg_mnxr_dict, mnxr_name_dict

        Parameters
        ------------------------------------------
        df: pandas.DataFrame
            DataFrame object which contains reac_xref information

        # Note : Ignore the information from the databases other than bigg and kegg for reactions
        """
        mnxr_name_dict = {} # 1:n for {key:value}
        mnxr_kegg_dict = {} # 1:n for {key:value}
        bigg_mnxr_dict = {} # 1:1 for {key:value}

        logging.debug('Parsing of reac_xref initiated...')
        for i in df.index:
            try:
                xref = df.loc[i, '#source']

                if xref.startswith(('bigg', 'kegg')):
                    xref_list = xref.split(':')
                    xref_db = xref_list[0]
                    xref_id = xref_list[1]
                    mnxr = df.loc[i, 'ID']
                    name = df.loc[i, 'description'].split('||')[0]

                    # For reaction.name in MNXref.xml	
                    if (xref_db == 'biggR' and not xref_id.startswith('R_')) or xref_db == 'keggR':
                        if mnxr not in mnxr_name_dict:
                            mnxr_name_dict[mnxr] = [name]
                        elif mnxr in mnxr_name_dict:
                            mnxr_name_dict[mnxr].append(name)

                    if xref_db == 'keggR':
                        if mnxr not in mnxr_kegg_dict:
                            mnxr_kegg_dict[mnxr] = [xref_id]
                        elif mnxr in mnxr_kegg_dict:
                            mnxr_kegg_dict[mnxr].append(xref_id)

                        logging.debug('%s; %s; %s' %(xref_db, xref_id, mnxr))

                    if xref_db == 'biggR' and not xref_id.startswith('R_'):
                        bigg_mnxr_dict[xref_id] = mnxr

                        logging.debug('%s; %s; %s' %(xref_db, xref_id, mnxr))

            except Exception as e:
                logging.debug(e)
                logging.debug('Cannot parse MNXR: %s' %df.loc[i,:].values)

        logging.debug('Cross reference dictionary for bigg reactions to MNXR has %d reactions' % len(bigg_mnxr_dict))
        logging.debug('Cross reference dictionary for MNXR reactions to kegg has %d reactions' % len(mnxr_kegg_dict))

        logging.debug('Parsing of reac_xref completed')

        self.mnxr_name_dict = mnxr_name_dict
        self.bigg_mnxr_dict = bigg_mnxr_dict

        return mnxr_kegg_dict, bigg_mnxr_dict


    def read_reac_prop(self, df):
        ''' Create reaction information dictionary 
        reac_prop has following fields :
        #ID	mnx_equation	reference	classifs	is_balanced	is_transport

        Output
        ------------------------------------------
        reaction_info : dict

        Parameters
        ------------------------------------------
        df: pandas.DataFrame
            DataFrame object which contains reac_xref information
        '''
        reaction_info = {}
        trans_count = 0

        logging.debug('Parsing of reac_prop initiated...')

        # skip first row since no information is in there
        for i in df.index[1:]:
            try:
                reaction_id = df.loc[i, '#ID'].strip()
                stoich_dict = {}
                equation = df.loc[i, 'mnx_equation']

                if df.loc[i, 'is_balanced'] == 'B':
                    mass_balance = 'balanced'
                else:
                    mass_balance = 'unbalanced'

                ec_number = df.loc[i, 'classifs'].strip()
                reactant_info, product_info = self.parse_equation(equation)

                for metabolite in reactant_info:
                    stoich_dict[metabolite] = reactant_info[metabolite]
                for metabolite in product_info:
                    stoich_dict[metabolite] = product_info[metabolite]

                # skip reactions in which same metabolite is in both reactants and products such as transport
                # Example: pass reaction like this: 'MNXM1@MNXMD1 <=> MNXM1@MNXMD2'
                trans_compartment = False
                for each_reactant in reactant_info.keys():
                    if each_reactant in product_info.keys():
                        trans_compartment = True
                        trans_count += 1
                        break

                if trans_compartment == True:
                    continue

                # select only reactions with ec_number information
                if ec_number:
                    ec_number_list = ec_number.split(';')
                    reaction_info[reaction_id] = {
                            'stoichiometry': stoich_dict,
                            'balance': mass_balance,
                            'ec': ec_number_list
                            }
            except Exception as e:
                logging.debug(e)
                logging.debug('Cannot parse MNXR: %s' %df.loc[i,:].values)

        logging.debug('Reaction info dictionary has %d reactions' % len(reaction_info))
        logging.debug('Translocation reactions: %d' %trans_count)

        logging.debug('Parsing of reac_prop completed')

        self.reaction_info = reaction_info
        return


    def parse_equation(self, equation):
        ''' Parse reaction equation to give reactant and products
        Equation form: chemID is replaced by chemID@compID

        e.g.
        1 MNXM12@MNXD1 + 1 MNXM146442@MNXD1 = 1 MNXM32694@MNXD1 + 1 MNXM686@MNXD1 

        NOTE : At the moment, we ignore the compartmentalization feature introduced in version 4

        Output
        ------------------------------------------
        reactant_info, product_info 
        '''
        reactant_info = {}
        product_info = {}

        spt_equation = equation.split(' = ')
        reactant_str, product_str = spt_equation[0], spt_equation[1]

        reactant_list = reactant_str.split(' + ')
        for metabolite in reactant_list:
            coeff_metab = metabolite.strip().split(' ')
            coeff = coeff_metab[0]
            metab = coeff_metab[1].split('@')[0]
            reactant_info[metab] = float(coeff) * -1

        product_list = product_str.split(' + ')
        for metabolite in product_list:
            coeff_metab = metabolite.strip().split(' ')
            coeff = coeff_metab[0]
            metab = coeff_metab[1].split('@')[0]
            product_info[metab] = float(coeff)

        return reactant_info, product_info


    def get_cobra_reactions(self):
        ''' Creates Universal list of cobra reaction objects 	

        Output
        ------------------------------------------
        cobra_reactions
            list of MNXR to be added to MNXref model
        '''
        logging.debug('Creating MNXR reactions: time-consuming')

        cobra_reactions = []

        cnt = 0
        for each_reaction in self.reaction_info:
            cnt += 1

            logging.debug('Total reaction number %s; Reaction number covered %s; %s' \
                    %(len(self.reaction_info), cnt, each_reaction))
            reaction_name = each_reaction
            coeff_dict = self.reaction_info[each_reaction]['stoichiometry']
            mass_balance = self.reaction_info[each_reaction]['balance']
            ec_number_list = self.reaction_info[each_reaction]['ec']

            new_reaction_stoichiometry = {}
            for each_metabolite in coeff_dict:
                # Convert MNXM to BiGG IDs
                if each_metabolite in self.mnxm_bigg_compound_dict.keys() \
                        and each_metabolite in self.mnxm_compoundInfo_dict.keys():
                    metabolite_id = '%s_c'    %self.mnxm_bigg_compound_dict[each_metabolite]
                    metabolite_obj = Metabolite(
                            str(metabolite_id),
                            name = self.mnxm_compoundInfo_dict[each_metabolite][0],
                            formula = self.mnxm_compoundInfo_dict[each_metabolite][1],
                            compartment = 'c')
                # use MNXM IDs if BiGG IDs are not found but compound information is saved in dictionary
                elif each_metabolite not in self.mnxm_bigg_compound_dict.keys() \
                        and each_metabolite in self.mnxm_compoundInfo_dict.keys():
                    metabolite_id = '%s_c'    %each_metabolite
                    metabolite_obj = Metabolite(
                            str(metabolite_id),
                            name = self.mnxm_compoundInfo_dict[each_metabolite][0],
                            formula = self.mnxm_compoundInfo_dict[each_metabolite][1],
                            compartment = 'c')
                # use MNXM IDs if none of information is available
                elif each_metabolite not in self.mnxm_bigg_compound_dict.keys() \
                        and each_metabolite not in self.mnxm_compoundInfo_dict.keys():
                        metabolite_id = '%s_c'    %each_metabolite
                        metabolite_obj = Metabolite(
                            str(metabolite_id),
                            name = '',
                            formula = '',
                            compartment = 'c')

                # Insert other db IDs
                metabolite_obj.notes = {}
                metabolite_obj.notes['MNXM'] = each_metabolite
                if each_metabolite in self.mnxm_bigg_compound_dict:
                    metabolite_obj.notes['BiGG'] = \
                            self.mnxm_bigg_compound_dict[each_metabolite]
                if each_metabolite in self.mnxm_kegg_compound_dict:
                    metabolite_obj.notes['KEGG'] = \
                            self.mnxm_kegg_compound_dict[each_metabolite]

                coeff = coeff_dict[each_metabolite]
                new_reaction_stoichiometry[metabolite_obj] = float(coeff)

            # create a new reaction object
            reaction_obj = Reaction(each_reaction)

            # Insert bigg or kegg reaction names
            if each_reaction in self.mnxr_name_dict:
                reaction_obj.name = ';'.join(self.mnxr_name_dict[each_reaction])
            else:
                reaction_obj.name = ''

            reaction_obj.subsystem = ''
            reaction_obj.lower_bound = -1000.0
            reaction_obj.upper_bound = 1000.0
            reaction_obj.reversibility = 1
            reaction_obj.gene_reaction_rule = ''
            reaction_obj.add_metabolites(new_reaction_stoichiometry)

            # Currently writing to sbml not supported
            reaction_obj.notes = {}
            reaction_obj.notes['EC_number'] = ';'.join(ec_number_list)
            reaction_obj.notes['Balance'] = mass_balance
            cobra_reactions.append(copy.deepcopy(reaction_obj))

        self.cobra_reactions = cobra_reactions


    def make_cobra_model(self):
        ''' Creates Universal cobra model object

        Output
        ------------------------------------------
        cobra_model
        '''
        cobra_model = Model('mnxref_model')

        for each_cobra_reaction in self.cobra_reactions:
            logging.debug("Adding reaction %s to 'mnxref_model'"
                          %each_cobra_reaction.id)
            cobra_model.add_reaction(each_cobra_reaction)

        for metabolite in cobra_model.metabolites:
            old_metab = metabolite.id
            metabolite.id = self.fix_legacy_id(old_metab)
            if old_metab != metabolite.id:
                logging.debug("Fixed metabolite ID in 'mnxref_model': %s -> %s",
                          old_metab, metabolite.id)

        cobra_model = gmsm.utils.stabilize_model(
                cobra_model, input2_tmp_dir, 'MNXref', diff_name=True)

        logging.debug('%i reactions in model' % len(cobra_model.reactions))
        logging.debug('%i metabolites in model' % len(cobra_model.metabolites))
        logging.debug('%i genes in model' % len(cobra_model.genes))

        return cobra_model


    def fix_legacy_id(self, id):
        id = id.replace('_DASH_', '__')
        id = id.replace('_FSLASH_', '/')
        id = id.replace('_BSLASH_', "\\")
        id = id.replace('_LPAREN_', '(')
        id = id.replace('_LSQBKT_', '[')
        id = id.replace('_RSQBKT_', ']')
        id = id.replace('_RPAREN_', ')')
        id = id.replace('_COMMA_', ',')
        id = id.replace('_PERIOD_', '.')
        id = id.replace('_APOS_', "'")
        id = id.replace('&amp;', '&')
        id = id.replace('&lt;', '<')
        id = id.replace('&gt;', '>')
        id = id.replace('&quot;', '"')
        id = id.replace("-", "__")
        return id


def unzip_tsv_files():
    tsv_files = glob.glob(join(input2_tmp_dir, '*.tsv'))
    if len(tsv_files) ==  0:
        zip = zipfile.ZipFile(join(input2_tmp_dir, 'mnxref.zip'))
        zip.extractall(input2_tmp_dir)
        zip.close()


def run_ParseMNXref():
    chem_xref = pd.read_csv(join(input2_tmp_dir, 'chem_xref_4.0.tsv'),
                            delimiter='\t', skiprows=347, keep_default_na=False)
    chem_prop = pd.read_csv(join(input2_tmp_dir, 'chem_prop_4.0.tsv'),
                            delimiter='\t', skiprows=347, keep_default_na=False)
    reac_xref = pd.read_csv(join(input2_tmp_dir, 'reac_xref_4.0.tsv'),
                            delimiter='\t', skiprows = 347, keep_default_na=False)
    reac_prop = pd.read_csv(join(input2_tmp_dir, 'reac_prop_4.0.tsv'),
                            delimiter='\t', skiprows= 347, keep_default_na=False)
    
    mnx_parser = ParseMNXref()
    bigg_old_new_dict = mnx_parser.fix_legacy_id_using_BiGGModels()
    mnxm_bigg_compound_dict = mnx_parser.read_chem_xref(bigg_old_new_dict, chem_xref)
    mnxm_compoundInfo_dict = mnx_parser.read_chem_prop(chem_prop)
    mnxr_kegg_dict, bigg_mnxr_dict = mnx_parser.read_reac_xref(reac_xref)
    mnx_parser.read_reac_prop(reac_prop)
    mnx_parser.get_cobra_reactions()
    cobra_model = mnx_parser.make_cobra_model()

    # Write SBML file
    cobra.io.write_sbml_model(cobra_model,
            join(input2_tmp_dir, 'MNXref.xml'))

    # Write txt files
    with open(join(input2_tmp_dir, 'mnxm_bigg_compound_dict.txt'), 'w') as f:
        for k, v in mnxm_bigg_compound_dict.items():
            print('%s\t%s' %(k, v), file=f)

    with open(join(input2_tmp_dir, 'bigg_old_new_dict.txt'), 'w') as f:
        for k, v in bigg_old_new_dict.items():
            print('%s\t%s' %(k, v), file=f)

    with open(join(input2_tmp_dir, 'mnxm_compoundInfo_dict.txt'), 'w') as f:
        for k, v in mnxm_compoundInfo_dict.items():
            print('%s\t%s' %(k, v), file=f)

    with open(join(input2_tmp_dir, 'mnxr_kegg_dict.txt'), 'w') as f:
        for k, v in mnxr_kegg_dict.items():
            print('%s\t%s' %(k, v), file=f)

    with open(join(input2_tmp_dir, 'bigg_mnxr_dict.txt'), 'w') as f:
        for k, v in bigg_mnxr_dict.items():
            print('%s\t%s' %(k, v), file=f)

    # Create pickles in 'scripts/input2_data' and 'input2'
    if not os.path.isfile(join(input2_tmp_dir, 'bigg_old_new_dict.p')):
        with open(join(input2_tmp_dir, 'bigg_old_new_dict.p'), 'wb') as f:
            pickle.dump(bigg_old_new_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input2_dir, 'mnxm_compoundInfo_dict.p'), 'wb') as f:
        pickle.dump(mnxm_compoundInfo_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input2_dir, 'mnxr_kegg_dict.p'), 'wb') as f:
        pickle.dump(mnxr_kegg_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input2_dir, 'bigg_mnxr_dict.p'), 'wb') as f:
        pickle.dump(bigg_mnxr_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(join(input2_dir, 'MNXref.p'), 'wb') as f:
        pickle.dump(cobra_model, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Copy pickles to the destination
    shutil.copyfile(join(input2_dir, 'MNXref.p'),
            join(os.pardir, 'gmsm', 'tests', 'data', 'MNXref.p'))


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

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level='DEBUG')
    logger = logging.getLogger('')
    fomatter = logging.Formatter(
                '[%(levelname)s|%(filename)s:%(lineno)s] > %(message)s')
    fh = logging.FileHandler('logfile.log', mode = 'w')
    fh.setFormatter(fomatter)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    unzip_tsv_files()
    run_ParseMNXref()
    input2_tmp_dir_list2 = create_zip_file()
    remove_tsv_files(input2_tmp_dir_list2)

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
