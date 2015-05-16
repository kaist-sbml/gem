'''
2015 Kyu-Sang Hwang
2014-2015 Hyun Uk Kim
'''

from Bio import SeqIO
from sets import Set
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file
from MNX_checker2 import COBRA_TemplateModel_checking_MNXref_metabolites
from MNX_checker2 import fix_legacy_id
import pickle
import copy

def determine_domain(domain_comb):
## domain information(nrps) :
## Condensation    Condensation domain
## Condensation_DCL    Condensation domain that links L-amino acid to peptide ending with D-amino acid
## Condensation_LCL    Condensation domain that links L-amino acid to peptide ending with L-amino acid 
## Condensation_Dual     Dual condensation / epimerization domain 
## Condensation_Starter    Starter condensation domain
## CXglyc    Putatively inactive glycopeptide condensation-like domain
## Cglyc    Glycopeptide condensation domain 
## Heterocyclization    Heterocyclization domain
## Epimerization    Epimerization domain 
## AMP-binding    Adenylation domain
## A-OX     Adenylation domain with integrated oxidase
## PCP     Peptidyl-carrier protein domain 
## ACPS    4'-phosphopantetheinyl transferase
## NRPS-COM_Nterm    NRPS COM domain Nterminal
## NRPS-COM_Cterm    NRPS COM domain Cterminal

## domain information(pks) :
## AT    acyltransferase
## KS    ketosynthase
## ACP    acyl carrier protein
## KR    ketoreductase
## DH    dehydratase
## ER    enolase
## cMT    methyltransferase
## TD    thiolesterase domain

    #Exceptionsal cases
    if 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'A'

    #Starter unit_AMP-binding+PCP
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'A_PCP' #

    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'A_PCP'

    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'MT_A_PCP' ##

    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'MT_A_PCP'

    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'MT_A_PCP'

    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'MT_A_PCP'

    #####
    #Condensation starter domain : C-domain
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_PCP' #

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_PCP' #

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_PCP_E' ####

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_PCP_E' ####

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP' ###

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP' ###

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP'

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cs_A_MT_E_PCP'

    #Condensation domain linking D-amino acid to peptides ending with L-amino acid : C-domain
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_PCP' #

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_PCP' #

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_PCP' ##

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_PCP_E' ####

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_PCP_E' ####

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'C_A_MT_E_PCP' ###

    #Condensation domain linking D-amino acid to peptides ending with L-amino acid: C-domain
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_PCP' #

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_PCP' #

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_PCP_E' ####

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_PCP_E' ####

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cdcl_A_MT_E_PCP' ###

    #Condensation domain that links L-amino acid to peptides ending with L-amino acid: C-domain
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_PCP' #

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_PCP' #

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_PCP' ##

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_PCP_E' ####

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_PCP_E' ####

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Clcl_A_MT_E_PCP' ###

    #Condensation_dual: C-domain
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cd_A_PCP' #### + epimerization

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cd_A_PCP' #### + epimerization

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_com:
        discriminant = 'Cd_A_MT_PCP' ###

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cd_A_MT_PCP' ###

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cd_A_MT_PCP' ###

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cd_A_MT_PCP' ###

    #Glycopeptide condensation domain (O): C-domain
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_PCP' #

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_PCP' #

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_PCP' ##

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_PCP_E' ####

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_PCP_E' ####

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###
 
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'Cglyc_A_MT_E_PCP' ###

    #Glycopeptide condensation domain (X): C-domain
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_PCP' #

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_PCP' #

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_PCP' ##

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_PCP_E' ####

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_PCP_E' ####

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'CXglyc_A_MT_E_PCP' ###

    #Heterocyclization: C-domain
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'PCP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and 'ACP' in domain_comb:
        discriminant = 'HC_Aox_MT_E_PCP'

    #t1pks modules: Should be domains, and NOT modules?
    #Starter units 
    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KR_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_DH_KR_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_DH_ER_KR_ACP'

    #Exceptional cases
    elif 'PKS_AT' in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_ER_KR_ACP'

    elif 'PKS_AT' not in domain_comb and 'PKS_KS' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'ACP'

    #Extension units
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_KR_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_KR'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_KR_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_KR'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_ER_KR_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_DH_ER_KR'

    #Exeptional cases
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_ER_KR_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' in domain_comb and 'cMT' not in domain_comb:
        discriminant = 'AT_KS_ER_KR'

    #Methytransferase
    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_cMT_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' not in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_cMT'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_KR_cMT_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_KR_cMT'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb) and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT_ER_ACP'

    elif 'PKS_AT' in domain_comb and 'PKS_KS' in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb and 'PKS_KR' in domain_comb and 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb and 'cMT' in domain_comb:
        discriminant = 'AT_KS_DH_KR_cMT_ER'

    #Terminal reductase domain: terminal domain       
    elif 'TD' in domain_comb and 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'TD'

    else:
        discriminant = 'None'

    return discriminant


def extract_substrate_information_nrps(each_sub_set, discriminator):
    sptline2 = each_sub_set.split(';')
    whole_substrate_info = sptline2[1]
    print each_sub_set
    print whole_substrate_info
    list_participated_sustrate = []

    participated_substrates = whole_substrate_info.split(':')
    sptSubstrates = participated_substrates[1]
    print sptSubstrates

    if ', ' not in sptSubstrates:
        print "insufficient substrate_info"
        discriminator = "false"
        return list_participated_sustrate

    substrates = sptSubstrates.split(', ')
    print substrates

    for each_substrate in substrates:
        sptSubstrate = each_substrate.split('(')
        participated_substrate = sptSubstrate[0].strip()
        print participated_substrate

        list_participated_sustrate.append(participated_substrate)

    return list_participated_sustrate, discriminator


def extract_substrate_information_pks(each_sub_set):
    sptline2 = each_sub_set.split(';')
    whole_substrate_info = sptline2[1]

    participated_substrates = whole_substrate_info.split(':')
    sptSubstrates = participated_substrates[1]
    substrates = sptSubstrates.split(', ')

    list_participated_sustrate = []

    for each_substrate in substrates:
        sptSubstrate = each_substrate.split('(')
        participated_substrate = sptSubstrate[0].strip()

        list_participated_sustrate.append(participated_substrate)

    return list_participated_sustrate
