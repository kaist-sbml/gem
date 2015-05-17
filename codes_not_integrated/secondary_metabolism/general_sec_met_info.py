'''
2015 Kyu-Sang Hwang
2014-2015 Hyun Uk Kim
'''

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


def get_biggid_from_aSid(each_substrate):

    if each_substrate == 'ala':
        met_name = 'ala_DASH_L'

    elif each_substrate == 'arg':
        met_name = 'arg_DASH_L'

    elif each_substrate == 'asn':
        met_name ='asn_DASH_L'

    elif each_substrate == 'asp':
        met_name = 'asp_DASH_L'

    elif each_substrate == 'cys':
        met_name = 'cys_DASH_L'

    elif each_substrate == 'gln':
        met_name = 'gln_DASH_L'

    elif each_substrate == 'glu':
        met_name = 'glu_DASH_L'

    elif each_substrate == 'gly':
        met_name = 'gly'

    elif each_substrate == 'his':
        met_name = 'his_DASH_L'

    elif each_substrate == 'leu':
        met_name = 'leu_DASH_L'

    elif each_substrate == 'lys':
        met_name = 'lys_DASH_L'

    elif each_substrate == 'met':
        met_name = 'met_DASH_L'

    elif each_substrate == 'phe':
        met_name = 'phe_DASH_L'

    elif each_substrate == 'pro':
        met_name = 'pro_DASH_L'

    elif each_substrate == 'ser':
        met_name = 'ser_DASH_L'

    elif each_substrate == 'thr':
        met_name = 'thr_DASH_L'

    elif each_substrate == 'trp':
        met_name = 'trp_DASH_L'

    elif each_substrate == 'tyr':
        met_name ='tyr_DASH_L'

    elif each_substrate == 'val':
        met_name = 'val_DASH_L'

    elif each_substrate == 'ile':
        met_name = 'ile_DASH_L'

    elif each_substrate == 'phg':
        met_name = 'phg_DASH_L'

    elif each_substrate == 'bht':
        met_name = 'bht_DASH_L'

    elif each_substrate == 'orn':
        met_name = 'orn'

    elif each_substrate == 'abu':
        met_name = 'abu'

    elif each_substrate == 'iva':
        met_name = 'iva'

    elif each_substrate == 'aad':
        met_name = 'L2aadp'

    elif each_substrate == 'hpg':
        met_name = 'hpg'

    elif each_substrate == 'dhb':
        met_name = '23dhb'

    elif each_substrate == 'dhpg':
        met_name = 'dhpg'

    elif each_substrate == 'hty':
        met_name = 'hty'

    elif each_substrate == 'cit':
        met_name = 'citr_DASH_L'

    elif each_substrate == 'pip':
        met_name = 'Lpipecol'

    elif each_substrate == 'b-ala':
        met_name = 'ala_DASH_B'

    elif each_substrate == 'dab':
        met_name = '24dab'

    elif each_substrate == 'phenylacetate' or each_substrate == 'Pha':
        met_name = 'pac'

    elif each_substrate == 'tcl':
        met_name = 'tcl'

    elif each_substrate == 'qa':
        met_name = 'qa'

    #t1pks substreate
    elif each_substrate == 'mal':
        met_name = 'malcoa'

    elif each_substrate == 'mmal':
        met_name = 'mmcoa_DASH_S'

    elif each_substrate == '2metbut':
        met_name = '2mbcoa'

    elif each_substrate == 'Ethyl_mal' or each_substrate == 'emal':
        met_name = 'emcoa_DASH_S'

    elif each_substrate == 'isobut':
        met_name = 'ibcoa'

    elif each_substrate == 'ace':
        met_name = 'accoa'

    elif each_substrate == 'prop':
        met_name = 'ppcoa'

    elif each_substrate == '3metbut':
        met_name = 'ivcoa'

    elif each_substrate == 'mxmal':
        met_name = 'mxmalacp'

    elif each_substrate == 'CHC-CoA':
        met_name = 'chccoa'

    elif each_substrate == 'N/A':
        met_name = 'N/A'

    else:
        print each_substrate
        raw_input('substrate_not_defined')

    return met_name

