
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
from cobra import Metabolite

def get_module_struct(domain_comb):
#TODO: category these lines according to core and optional domains
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
## ACP    4'-phosphopantetheinyl transferase
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

    #TODO: Place essential domains first, followed by optional domains
    #Exceptionsal cases
    if 'ACP' not in domain_comb \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation' not in domain_comb \
            and 'Condensation_Starter' not in domain_comb \
            and 'Condensation_DCL' not in domain_comb \
            and 'Condensation_LCL' not in domain_comb \
            and 'Condensation_Dual' not in domain_comb \
            and 'Cglyc' not in domain_comb \
            and 'CXglyc' not in domain_comb \
            and 'Epimerization' not in domain_comb \
            and 'Heterocyclization' not in domain_comb \
            and 'cMT' not in domain_comb \
            and 'nMT' not in domain_comb \
            and 'PCP' not in domain_comb:
        module_struct = 'A'

    #Starter unit_AMP-binding+PCP
    # Following combination applies to all cases
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation' not in domain_comb \
            and 'Condensation_Starter' not in domain_comb \
            and 'Condensation_DCL' not in domain_comb \
            and 'Condensation_LCL' not in domain_comb \
            and 'Condensation_Dual' not in domain_comb \
            and 'Cglyc' not in domain_comb \
            and 'CXglyc' not in domain_comb \
            and 'Epimerization' not in domain_comb \
            and 'Heterocyclization' not in domain_comb:
        if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
            module_struct = 'A_PCP'
        elif 'cMT' in domain_comb or 'nMT' in domain_comb:
            module_struct = 'MT_A_PCP'

    #Condensation starter domain : C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation_Starter' in domain_comb:
        if 'Epimerization' not in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Cs_A_PCP'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Cs_A_MT_PCP'
        elif 'Epimerization' in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Cs_A_PCP_E'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Cs_A_MT_E_PCP'

    #Condensation domain linking D-amino acid to peptides ending with L-amino acid: C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation' in domain_comb:
        if 'Epimerization' not in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'C_A_PCP'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'C_A_MT_PCP'
        elif 'Epimerization' in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'C_A_PCP_E'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'C_A_MT_E_PCP'

    #Condensation domain linking D-amino acid to peptides ending with L-amino acid: C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation_DCL' in domain_comb:
        if 'Epimerization' not in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Cdcl_A_PCP'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Cdcl_A_MT_PCP'
        elif 'Epimerization' in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Cdcl_A_PCP_E'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Cdcl_A_MT_E_PCP'

    #Condensation domain that links L-amino acid to peptides ending with L-amino acid: C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation_LCL' in domain_comb:
        if 'Epimerization' not in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Clcl_A_PCP'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Clcl_A_MT_PCP'
        elif 'Epimerization' in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Clcl_A_PCP_E'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Clcl_A_MT_E_PCP'

    #Condensation_dual: C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation_Dual' in domain_comb \
            and 'Epimerization' not in domain_comb:
        if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
            module_struct = 'Cd_A_PCP' #### + epimerization
        elif 'cMT' in domain_comb or 'nMT' in domain_comb:
            module_struct = 'Cd_A_MT_PCP'

    #Glycopeptide condensation domain (O): C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Cglyc' in domain_comb:
        if 'Epimerization' not in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Cglyc_A_PCP'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Cglyc_A_MT_PCP'
        elif 'Epimerization' in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'Cglyc_A_PCP_E'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'Cglyc_A_MT_E_PCP'

    #Glycopeptide condensation domain (X): C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'AMP-binding' in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'CXglyc' in domain_comb:
        if 'Epimerization' not in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'CXglyc_A_PCP'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'CXglyc_A_MT_PCP'
        elif 'Epimerization' in domain_comb:
            if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                module_struct = 'CXglyc_A_PCP_E'
            elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                module_struct = 'CXglyc_A_MT_E_PCP'

    #Heterocyclization: C-domain
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'Heterocyclization' in domain_comb:
        if 'AMP-binding' in domain_comb and 'A-OX' not in domain_comb:
            if 'Epimerization' not in domain_comb:
                if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                    module_struct = 'HC_A_PCP'
                elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                    module_struct = 'HC_A_MT_PCP'
            elif 'Epimerization' in domain_comb:
                module_struct = 'HC_A_MT_E_PCP'
        elif 'AMP-binding' not in domain_comb and 'A-OX' in domain_comb:
            if 'Epimerization' not in domain_comb:
                if 'cMT' not in domain_comb and 'nMT' not in domain_comb:
                    module_struct = 'HC_Aox_PCP'
                elif 'cMT' in domain_comb or 'nMT' in domain_comb:
                    module_struct = 'HC_Aox_MT_PCP'
            elif 'Epimerization' in domain_comb:
                module_struct = 'HC_Aox_MT_E_PCP'

    #Terminal reductase domain: terminal domain
    elif 'ACP' not in domain_comb and 'PCP' not in domain_comb \
            and 'AMP-binding' not in domain_comb \
            and 'A-OX' not in domain_comb \
            and 'Condensation' not in domain_comb \
            and 'Condensation_Starter' not in domain_comb \
            and 'Condensation_DCL' not in domain_comb \
            and 'Condensation_LCL' not in domain_comb \
            and 'Condensation_Dual' not in domain_comb \
            and 'Cglyc' not in domain_comb \
            and 'CXglyc' not in domain_comb \
            and 'Epimerization' not in domain_comb \
            and 'Heterocyclization' not in domain_comb \
            and 'cMT' not in domain_comb \
            and 'nMT' not in domain_comb \
            and 'TD' in domain_comb:
        module_struct = 'TD'

    #Starter units
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' not in domain_comb \
            and 'cMT' not in domain_comb:
        if 'PKS_KR' not in domain_comb:
            if 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb:
                module_struct = 'AT_ACP'
        elif 'PKS_KR' in domain_comb:
            if 'PKS_DH' not in domain_comb:
                if 'PKS_ER' not in domain_comb:
                    module_struct = 'AT_KR_ACP'
            elif 'PKS_DH' in domain_comb:
                if 'PKS_ER' not in domain_comb:
                    module_struct = 'AT_DH_KR_ACP'
                if 'PKS_ER' in domain_comb:
                    module_struct = 'AT_DH_ER_KR_ACP'

    #Exceptional cases
    # Following cases cannot be further reduced with if statements
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'PKS_AT' not in domain_comb \
            and 'PKS_KS' not in domain_comb \
            and 'PKS_KR' not in domain_comb \
            and 'PKS_DH' not in domain_comb \
            and 'PKS_ER' not in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'ACP'

    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' not in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' not in domain_comb \
            and 'PKS_ER' in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_ER_KR_ACP'

    #Extension units
    # Following cases cannot be further reduced with if statements
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and  'PKS_KR' not in domain_comb \
            and 'PKS_DH' not in domain_comb \
            and 'PKS_ER' not in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS_ACP'
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' not in domain_comb \
            and 'PKS_ER' not in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS_KR_ACP'
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' in domain_comb \
            and 'PKS_ER' not in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS_DH_KR_ACP'
    elif ('ACP' in domain_comb or 'PCP' in domain_comb) \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' in domain_comb \
            and 'PKS_ER' in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS_DH_ER_KR_ACP'

    elif 'ACP' not in domain_comb \
            and 'PCP' not in domain_comb \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' not in domain_comb \
            and 'PKS_DH' not in domain_comb \
            and 'PKS_ER' not in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS'
    elif 'ACP' not in domain_comb \
            and 'PCP' not in domain_comb \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' not in domain_comb \
            and 'PKS_ER' not in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS_KR'
    elif 'ACP' not in domain_comb \
            and 'PCP' not in domain_comb \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' in domain_comb \
            and 'PKS_ER' not in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS_DH_KR'
    elif 'ACP' not in domain_comb \
            and 'PCP' not in domain_comb \
            and 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' in domain_comb \
            and 'PKS_ER' in domain_comb \
            and 'cMT' not in domain_comb:
        module_struct = 'AT_KS_DH_ER_KR'

    #Exeptional cases
    elif 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'PKS_KR' in domain_comb \
            and 'PKS_DH' not in domain_comb \
            and'PKS_ER' in domain_comb \
            and 'cMT' not in domain_comb:
        if 'PCP' not in domain_comb and 'ACP' not in domain_comb:
            module_struct = 'AT_KS_ER_KR'
        elif ('PCP' in domain_comb or 'ACP' in domain_comb):
            module_struct = 'AT_KS_ER_KR_ACP'

    #Methytransferase
    elif 'PKS_AT' in domain_comb \
            and 'PKS_KS' in domain_comb \
            and 'cMT' in domain_comb:
        if 'PCP' not in domain_comb and 'ACP' not in domain_comb:
            if 'PKS_KR' not in domain_comb:
                if 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb:
                    module_struct = 'AT_KS_cMT'
            elif 'PKS_KR' in domain_comb:
                if 'PKS_DH' not in domain_comb:
                    if 'PKS_ER' not in domain_comb:
                        module_struct = 'AT_KS_KR_cMT'
                elif 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb:
                    module_struct = 'AT_KS_DH_KR_cMT'
                elif 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb:
                    module_struct = 'AT_KS_DH_KR_cMT_ER'
        elif ('PCP' in domain_comb or 'ACP' in domain_comb):
            if 'PKS_KR' not in domain_comb:
                if 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb:
                    module_struct = 'AT_KS_cMT_ACP'
            if 'PKS_KR' in domain_comb:
                if 'PKS_DH' not in domain_comb and 'PKS_ER' not in domain_comb:
                    module_struct = 'AT_KS_KR_cMT_ACP'
                elif 'PKS_DH' in domain_comb and 'PKS_ER' not in domain_comb:
                    module_struct = 'AT_KS_DH_KR_cMT_ACP'
                elif 'PKS_DH' in domain_comb and 'PKS_ER' in domain_comb:
                    module_struct = 'AT_KS_DH_KR_cMT_ER_ACP'


    else:
        module_struct = 'None'
    logging.debug('domain_comb: %s' %domain_comb)
    logging.debug('module_struct: %s' %module_struct)
    return module_struct


def get_kr_activity(each_locustag, domain_comb, locustag_kr_dict, module_struct):

    domain_kr_activity_dict = locustag_kr_dict[each_locustag]

    for each_domain in domain_comb:
        domain_name = each_domain[:-5]

        if domain_name == "PKS_KR" and domain_kr_activity_dict[each_domain] == 'inactive':
            module_struct_kr = module_struct.replace('KR','KR(inactive)')
            return module_struct_kr
        else:
            module_struct_kr = module_struct
            return module_struct_kr

    return module_struct_kr


def get_module_currency_metab_dict(module_struct, each_module, each_module_substrates, module_currency_metab_dict):

    if module_struct == 'A' or module_struct == 'Aox':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['h2o'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'A or Aox:', each_module_substrates

    elif module_struct == 'A_PCP' or module_struct == 'Cs_A_PCP' or module_struct == 'C_A_PCP' or module_struct == 'Cdcl_A_PCP' or module_struct == 'Clcl_A_PCP' or module_struct == 'Clcl_A_PCP' or module_struct == 'Cglyc_A_PCP' or module_struct == 'CXglyc_A_PCP':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['h2o'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'A_PCP or C_A_PCP:', each_module_substrates

    elif module_struct == 'MT_A_PCP' or module_struct == 'Cs_A_MT_PCP' or module_struct == 'C_A_MT_PCP' or module_struct == 'Cdcl_A_MT_PCP' or module_struct == 'Clcl_A_MT_PCP' or module_struct == 'Cglyc_A_MT_PCP' or module_struct == 'CXglyc_A_MT_PCP':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        each_module_substrates['h2o'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'MT-A-PCP or C-A-MT-PCP:', each_module_substrates

    elif module_struct == 'Cs_A_PCP_E' or module_struct == 'C_A_PCP_E' or module_struct == 'Cdcl_A_PCP_E' or module_struct == 'Clcl_A_PCP_E' or module_struct == 'Cd_A_PCP' or module_struct == 'Cglyc_A_PCP_E' or module_struct == 'CXglyc_A_PCP_E':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['h2o'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'C_A_PCP_E:', each_module_substrates

    elif module_struct == 'Cs_A_MT_E_PCP' or module_struct == 'C_A_MT_E_PCP' or module_struct == 'Cdcl_A_MT_E_PCP' or module_struct == 'Clcl_A_MT_E_PCP' or module_struct == 'Cglyc_A_MT_E_PCP' or module_struct == 'CXglyc_A_MT_E_PCP' or module_struct == 'Cd_A_MT_PCP':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        each_module_substrates['h2o'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'C_A_MT_E_PCP:', each_module_substrates

    elif module_struct == 'HC_A_PCP':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['h2o'] = 2
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'HC_A_PCP:', each_module_substrates

    elif module_struct == 'HC_Aox_PCP':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        each_module_substrates['fmn'] = 1
        each_module_substrates['fmnh2'] = -1
        each_module_substrates['h2o'] = 2
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'HC_Aox_PCP:', each_module_substrates

    elif module_struct == 'HC_A_MT_PCP' or module_struct == 'HC_A_MT_E_PCP':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['ppi'] = 1
        each_module_substrates['h2o'] = 2
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'HC_A_MT_PCP:', each_module_substrates

    elif module_struct == 'HC_Aox_MT_PCP' or module_struct == 'HC_Aox_MT_E_PCP':
        each_module_substrates['atp'] = -1
        each_module_substrates['amp'] = 1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        each_module_substrates['fmn'] = 1
        each_module_substrates['fmnh2'] = -1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        each_module_substrates['h2o'] = 2
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'HC_Aox_MT_PCP:', each_module_substrates

    elif module_struct == 'AT_ACP' or module_struct == 'AT_KR(inactive)_ACP' or module_struct == 'AT_DH_KR(inactive)_ACP' or module_struct == 'AT_DH_ER_KR(inactive)_ACP':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'AT-ACP:', each_module_substrates

    elif module_struct == 'AT_KS_ACP' or module_struct == 'AT_KS' or module_struct == 'AT_KS_KR(inactive)_ACP' or module_struct == 'AT_KS_KR(inactive)' or module_struct == 'AT_KS_DH_KR(inactive)_ACP' or module_struct == 'AT_KS_DH_KR(inactive)' or module_struct == 'AT_KS_DH_KR(inactive)_ACP' or module_struct == 'AT_KS_DH_ER_KR(inactive)':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'AT-KS-ACP:', each_module_substrates

    elif module_struct == 'AT_KS_KR_ACP' or module_struct == 'AT_KR_ACP' or module_struct == 'AT_KS_KR' or module_struct == 'AT_ER_KR_ACP' or module_struct == 'AT_KS_ER_KR' or module_struct == 'AT_KS_ER_KR_ACP':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        each_module_substrates['nadp'] = 1
        each_module_substrates['nadph'] = -1
        each_module_substrates['h'] = -1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'KS-AT-KR-ACP:', each_module_substrates

    elif module_struct == 'AT_KS_DH_KR_ACP' or module_struct == 'AT_DH_KR_ACP' or module_struct == 'AT_KS_DH_KR':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        each_module_substrates['h2o'] = 1
        each_module_substrates['nadp'] = 1
        each_module_substrates['nadph'] = -1
        each_module_substrates['h'] = -1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'KS-AT-DH-KR-ACP:', each_module_substrates

    elif module_struct == 'AT_KS_DH_ER_KR_ACP' or module_struct == 'AT_DH_ER_KR_ACP' or module_struct == 'AT_KS_DH_ER_KR':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        each_module_substrates['h2o'] = 1
        each_module_substrates['nadp'] = 2
        each_module_substrates['nadph'] = -2
        each_module_substrates['h'] = -2
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'KS-AT-DH-ER-KR-ACP:', each_module_substrates

    elif module_struct == 'AT_KS_cMT_ACP' or module_struct == 'AT_KS_cMT':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'AT_KS_cMT_ACP:', each_module_substrates

    elif module_struct == 'AT_KS_KR_cMT_ACP' or module_struct == 'AT_KS_KR_cMT':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        each_module_substrates['nadp'] = 1
        each_module_substrates['nadph'] = -1
        each_module_substrates['h'] = -1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'AT_KS_KR_cMT_ACP:', each_module_substrates

    elif module_struct == 'AT_KS_DH_KR_cMT_ACP' or module_struct == 'AT_KS_DH_KR_cMT':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        each_module_substrates['h2o'] = 1
        each_module_substrates['nadp'] = 1
        each_module_substrates['nadph'] = -1
        each_module_substrates['h'] = -1
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'AT_KS_DH_KR_cMT_ACP:', each_module_substrates

    elif module_struct == 'AT_KS_DH_KR_cMT_ER_ACP' or module_struct == 'AT_KS_DH_KR_cMT_ER':
        each_module_substrates['coa'] = 1
        each_module_substrates['hco3'] = 1
        each_module_substrates['h2o'] = 1
        each_module_substrates['nadp'] = 2
        each_module_substrates['nadph'] = -2
        each_module_substrates['h'] = -2
        each_module_substrates['amet'] = -1
        each_module_substrates['ahcys'] = 1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'AT_KS_DH_KR_cMT_ER_ACP:', each_module_substrates

    elif module_struct == 'TD':
        each_module_substrates['nadp'] = 1
        each_module_substrates['nadph'] = -1
        each_module_substrates['h'] = -1
        module_currency_metab_dict[each_module] = each_module_substrates
        #print 'reaction 8: HC_Aox_MT_PCP', each_module_substrates

    return module_currency_metab_dict


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

    #No MNXM, KEGG ID and bigg ID for "bht"
    elif each_substrate == 'bht':
        met_name = 'bht_DASH_L'

    elif each_substrate == 'orn':
        met_name = 'orn'

    #No bigg ID
    elif each_substrate == 'abu':
        met_name = 'MNXM17054'

    #No bigg ID
    elif each_substrate == 'iva':
        met_name = 'MNXM34821'

    elif each_substrate == 'aad':
        met_name = 'L2aadp'

    #No bigg ID
    elif each_substrate == 'hpg':
        met_name = 'MNXM4544'

    elif each_substrate == 'dhb':
        met_name = '23dhb'

    #No bigg ID
    elif each_substrate == 'dhpg':
        met_name = 'MNXM9962'

    #No bigg ID
    elif each_substrate == 'hty':
        met_name = 'MNXM59438'

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

    elif each_substrate == '2-3-diaminoproprionate':
        met_name = '23dappa'

    elif each_substrate == 'thr-4-cl':
        met_name = 'MNXM37380'

    #No bigg ID
    elif each_substrate == 'tcl':
        met_name = 'MNXM37380'

    #No bigg ID
    elif each_substrate == 'qa':
        met_name = 'MNXM80505'

    #No bigg ID
    elif each_substrate == 'trans-1,2-CPDA':
        met_name = '23cpda'

    #No bigg ID
    elif each_substrate == 'bmt':
        met_name = 'MNXM31446'

    elif each_substrate == 'sal':
        met_name = 'salc'

    #No bigg ID
    elif each_substrate == 'alaninol':
        met_name = 'MNXM8817'

    #t1pks substreate
    elif each_substrate == 'mal':
        met_name = 'malcoa'

    elif each_substrate == 'mmal':
        met_name = 'mmcoa_DASH_S'

    elif each_substrate == '2metbut':
        met_name = '2mbcoa'

    #Not available in bigg, but available in iMK1208
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

    #No bigg ID
    elif each_substrate == 'mxmal':
        met_name = 'MNXM61686'

    #No bigg ID
    elif each_substrate == 'CHC-CoA':
        met_name = 'MNXM5111'

    elif each_substrate == 'cemal':
        met_name = 'MNXM10927'

    #Taken from Supplementary Table of Minowa et al.
    elif each_substrate == 'qna':
        met_name = 'MNXM4797'

    #No bigg ID
    elif each_substrate == 'benz':
        met_name = 'MNXM240'

    #No bigg ID
    elif each_substrate == 'capreomycidine':
        met_name = 'MNXM18891'

    return met_name


def get_metab_coeff_dict():

    #Currency metabolites
    metab_coeff_dict = {}
    metab_coeff_dict['atp'] = 0
    metab_coeff_dict['amp'] = 0
    metab_coeff_dict['ppi'] = 0
    metab_coeff_dict['amet'] = 0
    metab_coeff_dict['ahcys'] = 0
    metab_coeff_dict['fmn'] = 0
    metab_coeff_dict['fmnh2'] = 0
    metab_coeff_dict['nadp'] = 0
    metab_coeff_dict['nadph'] = 0
    metab_coeff_dict['h'] = 0
    metab_coeff_dict['h2o'] = 0
    metab_coeff_dict['hco3'] = 0
    metab_coeff_dict['coa'] = 0

    metab_coeff_dict['ala_DASH_L'] = 0 #'L-alanine', 'C00041', 'MNXM32'
    metab_coeff_dict['arg_DASH_L'] = 0 #'L-Arginine', 'C00062', 'MNXM70'
    metab_coeff_dict['asn_DASH_L'] = 0 #'L-asparagine', 'C00152', 'MNXM147'
    metab_coeff_dict['asp_DASH_L'] = 0 #'L-Aspartate', 'C00049', 'MNXM42'
    metab_coeff_dict['cys_DASH_L'] = 0  #'L-Cysteine', 'C00097', 'MNXM55'
    metab_coeff_dict['gln_DASH_L'] = 0 #'L-Glutamine', 'C00064', 'MNXM37'
    metab_coeff_dict['glu_DASH_L'] = 0 #'L-Glutamate', 'C00025', 'MNXM89557'
    metab_coeff_dict['gly'] = 0 #'glycine', 'C00037', 'MNXM29'
    metab_coeff_dict['his_DASH_L'] = 0 #'L-Histidine', 'C00135', 'MNXM134'
    metab_coeff_dict['leu_DASH_L'] = 0 #'L-Leucine', 'C00123', 'MNXM140'
    metab_coeff_dict['lys_DASH_L'] = 0 #'L-Lysine', 'C00047 ', 'MNXM78'
    metab_coeff_dict['met_DASH_L'] = 0 #'L-Methionine', 'C00073', 'MNXM61'
    metab_coeff_dict['phe_DASH_L'] = 0 #'L-Phenylalanine', 'C00079', 'MNXM97'
    metab_coeff_dict['pro_DASH_L'] = 0 #'L-Proline', 'C00148', 'MNXM114'
    metab_coeff_dict['ser_DASH_L'] = 0 #'L-Serine', 'C00065', 'MNXM53'
    metab_coeff_dict['thr_DASH_L'] = 0 #'L-Threonine', 'C00188', 'MNXM142'
    metab_coeff_dict['trp_DASH_L'] = 0 #'L-Tryptophan', 'C00078', 'MNXM94'
    metab_coeff_dict['tyr_DASH_L'] = 0 #'L-Tyrosine', 'C00082', 'MNXM76'
    metab_coeff_dict['val_DASH_L'] = 0 #'L-Valine', 'C00183', 'MNXM199'
    metab_coeff_dict['ile_DASH_L'] = 0 #'L-Isoleucine', 'C00407', 'MNXM231'
    metab_coeff_dict['phg_DASH_L'] = 0 #'phenylglycine', 'C18623', 'MNXM59292'
    metab_coeff_dict['bht_DASH_L'] = 0 #'beta-hydroxyl-tyrosine', 'N/A', 'N/A'
    metab_coeff_dict['orn'] = 0 #'Ornithine', 'C01602', 'MNXM89689'
    metab_coeff_dict['MNXM17054'] = 0 #'abu', 'D-2-Aminobutyric acid', 'C02261', 'MNXM17054'
    metab_coeff_dict['MNXM34821'] = 0 #'iva', '2-Amino-2-methylbutanoate', 'C03571', 'MNXM34821'
    metab_coeff_dict['L2aadp'] = 0 #'L-2-Aminoadipic acid', 'C00956', 'MNXM268'
    metab_coeff_dict['MNXM4544'] = 0 #'hpg', 'D-4-Hydroxyphenylglycine', 'C03493', 'MNXM4544'
    metab_coeff_dict['23dhb'] = 0 #'2,3-Dihydroxybenzoic acid', 'C00196', 'MNXM455'
    metab_coeff_dict['MNXM9962'] = 0 #'dhpg', '3,5-Dihydroxy-phenylglycine', 'C12026', 'MNXM9962'
    metab_coeff_dict['MNXM59438'] = 0 #'hty', 'L-Homotyrosine', 'C18622', 'MNXM59438'
    metab_coeff_dict['citr_DASH_L'] = 0 #'L-citruline', 'C00327', 'MNXM211'
    metab_coeff_dict['Lpipecol'] = 0 #'L-pipecolate', 'C00408', 'MNXM684'
    metab_coeff_dict['ala_DASH_B'] = 0 #'beta-alanine zwitterion', 'C00099', 'MNXM144'
    metab_coeff_dict['24dab'] = 0 #'L-2,4-diazaniumylbutyrate', 'C03283', 'MNXM840'
    metab_coeff_dict['pac'] = 0 #'phenylacetate', 'C00548', 'MNXM497'
    metab_coeff_dict['23dappa'] = 0 #'3-aminoalanine zwitterion', 'C06393', 'MNXM91374'
    metab_coeff_dict['MNXM37380'] = 0 #'tcl', '4-Chlorothreonine', 'N/A', 'MNXM37380'
    metab_coeff_dict['tcl'] = 0 #(4S)-5,5,5-trichloro-leucine', 'N/A','N/A'
    metab_coeff_dict['MNXM80505'] = 0 #'qa', 'quinoxaline', 'C18575','MNXM80505'
    metab_coeff_dict['23cpda'] = 0 #'Trans-cyclopentane-(1R, 2R)-dicarboxylic acid', 'N/A','N/A'
    metab_coeff_dict['MNXM31446'] = 0 #'2-Butenyl-4-methyl-threonine', 'C12029','MNXM31446'
    metab_coeff_dict['salc'] = 0 #'salicylate', 'C00805','MNXM378'
    metab_coeff_dict['MNXM8817'] = 0 #'L-alaninol', 'N/A','MNXM8817'
    metab_coeff_dict['malcoa'] = 0 #'malonyl-CoA', 'C00083', 'MNXM40'
    metab_coeff_dict['mmcoa_DASH_S'] = 0 #'(S)-methylmalonyl-CoA(5-)','C00683', 'MNXM190'
    metab_coeff_dict['2mbcoa'] = 0 #'2-methylbutanoyl-CoA', C01033,'MNXM569'
    metab_coeff_dict['emcoa_DASH_S'] = 0 #'ethylmalonyl-CoA','C18026', 'MNXM2043'
    metab_coeff_dict['ibcoa'] = 0 #'2-Methylpropanoyl-CoA', 'C00630', 'MNXM470'
    metab_coeff_dict['accoa'] = 0 #'Acetyl-CoA', 'C00024', 'MNXM21'
    metab_coeff_dict['ppcoa'] = 0 #'Propionyl-CoA', 'C00100', 'MNXM86'
    metab_coeff_dict['ivcoa'] = 0 #'3-Methylbutanoyl-CoA', 'C02939', 'MNXM471'
    metab_coeff_dict['MNXM61686'] = 0 #'mxmalacp', 'Methoxymalonyl-[acp]', 'C18616', 'MNXM61686'
    metab_coeff_dict['MNXM5111'] = 0 #'chccoa', 'cyclohexane-1-carboxyl-CoA', 'C09823', 'MNXM5111'
    metab_coeff_dict['MNXM10927'] = 0 #'chloroethylmalonyl-CoA','N/A', 'MNXM10927'
    metab_coeff_dict['MNXM4797'] = 0 #'Quinaldinic acid','C06325', 'MNXM4797'
    metab_coeff_dict['MNXM240'] = 0 #'Benzoyl-CoA','C00512', 'MNXM240'
    metab_coeff_dict['MNXM18891'] = 0 #'L-Capreomycidine', 'C18472', 'MNXM18891'
    return metab_coeff_dict


#Add metabolite MNXM having no bigg ID to the model
def add_sec_met_mnxm_having_no_biggid_to_model(
        metab, metab_compt, mnxm_compoundInfo_dict):

    if metab == 'phg_DASH_L':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM59292'][1],
                name = mnxm_compoundInfo_dict['MNXM59292'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    #No MNXM, KEGG ID and bigg ID for "bht"
    elif metab == 'bht_DASH_L':
        metab_compt = Metabolite(metab_compt, compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM17054': #'abu'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM17054'][1],
                name = mnxm_compoundInfo_dict['MNXM17054'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM34821': #'iva'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM34821'][1],
                name = mnxm_compoundInfo_dict['MNXM34821'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM4544': #'hpg'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM4544'][1],
                name = mnxm_compoundInfo_dict['MNXM4544'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM9962': #'dhpg'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM9962'][1],
                name = mnxm_compoundInfo_dict['MNXM9962'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM59438': #'hty'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM59438'][1],
                name = mnxm_compoundInfo_dict['MNXM59438'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM37380': #'tcl'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM37380'][1],
                name = mnxm_compoundInfo_dict['MNXM37380'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    #No MNXM, KEGG ID and bigg ID for "bht"
    elif metab == 'tcl':
        metab_compt = Metabolite(metab_compt, compartment='c')

    elif metab == 'MNXM80505': #'qa'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM80505'][1],
                name = mnxm_compoundInfo_dict['MNXM80505'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == '23cpda':
        metab_compt = Metabolite(metab_compt, compartment='c')

    elif metab == 'MNXM31446':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM31446'][1],
                name = mnxm_compoundInfo_dict['MNXM31446'][0],
                compartment='c')

    elif metab == 'MNXM8817':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM8817'][1],
                name = mnxm_compoundInfo_dict['MNXM8817'][0],
                compartment='c')

    elif metab == 'emcoa_DASH_S':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM2043'][1],
                name = mnxm_compoundInfo_dict['MNXM2043'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM61686': #'mxmalacp'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM61686'][1],
                name = mnxm_compoundInfo_dict['MNXM61686'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM5111': #'chccoa'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM5111'][1],
                name = mnxm_compoundInfo_dict['MNXM5111'][0],
                compartment='c')
        #rxn.add_metabolites({metab_compt:metab_coeff_dict[metab]})

    elif metab == 'MNXM10927':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM10927'][1],
                name = mnxm_compoundInfo_dict['MNXM10927'][0],
                compartment='c')

    elif metab == 'MNXM4797':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM4797'][1],
                name = mnxm_compoundInfo_dict['MNXM4797'][0],
                compartment='c')

    elif metab == 'MNXM240':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM240'][1],
                name = mnxm_compoundInfo_dict['MNXM240'][0],
                compartment='c')

    elif metab == 'MNXM18891':
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM18891'][1],
                name = mnxm_compoundInfo_dict['MNXM18891'][0],
                compartment='c')

    return metab_compt
