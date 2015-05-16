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
## domain information :
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

    #Exceptionsal cases
    if ('Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb) and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'A'

    #Starter unit_AMP-binding+PCP
    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'A_PCP' #

    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'MT_A_PCP' ##

    elif 'Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'MT_A_PCP'

    #Condensation starter domain : C-domain
    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_PCP' #   

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_PCP' ##

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_PCP_E' ####

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_E_PCP' ###

    elif 'Condensation_Starter' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cs_A_MT_E_PCP'

    #Condensation domain that links D-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_PCP' #
    
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_PCP' ##
        
    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_PCP' ##

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_PCP_E' ####

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_E_PCP' ###

    elif 'Condensation' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'C_A_MT_E_PCP' ###        

    #Condensation domain that links D-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_PCP' #

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_PCP' ##
        
    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_PCP' ##

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_PCP_E' ####

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_E_PCP' ###

    elif 'Condensation_DCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cdcl_A_MT_E_PCP' ###

    #Condensation domain that links L-amino acid to peptide ending with L-amino acid  : C-domain
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_PCP' #

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_PCP' ##

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_PCP' ##
 
    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_PCP_E' ####

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_E_PCP' ###

    elif 'Condensation_LCL' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Clcl_A_MT_E_PCP' ###

    #Condensation_dual : C-domain
    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cd_A_PCP' #### + epimerization

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cd_A_MT_PCP' ###

    elif 'Condensation_Dual' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cd_A_MT_PCP' ###

     #Glycopeptide condensation domain (O) : C-domain
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_PCP' #
  
    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_PCP' ##

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_PCP' ##

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_PCP_E' ####

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_E_PCP' ###

    elif 'Cglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'Cglyc_A_MT_E_PCP' ###

    #Glycopeptide condensation domain (X) : C-domain
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_PCP' #
   
    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_PCP' ##

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_PCP' ##

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_PCP_E' ####

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_E_PCP' ###

    elif 'CXglyc' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'CXglyc_A_MT_E_PCP' ###

    #Heterocyclization : C-domain
    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' not in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_A_MT_E_PCP' ###

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' in domain_comb and 'nMT' not in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_E_PCP'  

    elif 'Heterocyclization' in domain_comb and 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb and 'cMT' not in domain_comb and 'nMT' in domain_comb and 'Epimerization' in domain_comb and ('PCP' in domain_comb or 'ACP' in domain_comb):
        discriminant = 'HC_Aox_MT_E_PCP'  

    #Terminal reductase domain : terminal domain
    elif 'TD' in domain_comb and ('Condensation_Starter' not in domain_comb and 'Condensation' not in domain_comb and 'Condensation_DCL' not in domain_comb and 'Condensation_LCL' not in domain_comb and 'Condensation_Dual' not in domain_comb and 'Cglyc' not in domain_comb and 'CXglyc' not in domain_comb and 'Heterocyclization' not in domain_comb) and 'AMP-binding' not in domain_comb and 'A-OX' not in domain_comb and 'cMT' not in domain_comb and 'nMT' not in domain_comb and 'Epimerization' not in domain_comb and 'PCP' not in domain_comb and 'ACP' not in domain_comb:
        discriminant = 'TD'

    else:
        discriminant = 'None'
        
    return discriminant

# def Identifier_KR_activity(discriminant, each_module_KR_activity):
#     
#     print discriminant, each_module_KR_activity
# 
#     if each_module_KR_activity == 1 or each_module_KR_activity == 0:
#         discriminant_with_KRact = discriminant
#         return discriminant_with_KRact
#     
#     elif each_module_KR_activity == 2:
#         discriminant_with_KRact = discriminant.replace('KR','KR(inactive)')      
#         return discriminant_with_KRact

