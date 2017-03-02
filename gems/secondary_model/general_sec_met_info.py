
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
from cobra import Metabolite

def get_biggid_from_aSid(each_substrate):

    #'L-alanine', 'C00041', 'MNXM32'
    if each_substrate == 'ala':
        met_name = 'ala_DASH_L'

    #'L-Arginine', 'C00062', 'MNXM70'
    elif each_substrate == 'arg':
        met_name = 'arg_DASH_L'

    #'L-asparagine', 'C00152', 'MNXM147'
    elif each_substrate == 'asn':
        met_name ='asn_DASH_L'

    #'L-Aspartate', 'C00049', 'MNXM42'
    elif each_substrate == 'asp':
        met_name = 'asp_DASH_L'

    #'L-Cysteine', 'C00097', 'MNXM55'
    elif each_substrate == 'cys':
        met_name = 'cys_DASH_L'

    #'L-Glutamine', 'C00064', 'MNXM37'
    elif each_substrate == 'gln':
        met_name = 'gln_DASH_L'

    #'L-Glutamate', 'C00025', 'MNXM89557'
    elif each_substrate == 'glu':
        met_name = 'glu_DASH_L'

    #'glycine', 'C00037', 'MNXM29'
    elif each_substrate == 'gly':
        met_name = 'gly'

    #'L-Histidine', 'C00135', 'MNXM134'
    elif each_substrate == 'his':
        met_name = 'his_DASH_L'

    #'L-Leucine', 'C00123', 'MNXM140'
    elif each_substrate == 'leu':
        met_name = 'leu_DASH_L'

    #'L-Lysine', 'C00047 ', 'MNXM78'
    elif each_substrate == 'lys':
        met_name = 'lys_DASH_L'

    #'L-Methionine', 'C00073', 'MNXM61'
    elif each_substrate == 'met':
        met_name = 'met_DASH_L'

    #'L-Phenylalanine', 'C00079', 'MNXM97'
    elif each_substrate == 'phe':
        met_name = 'phe_DASH_L'

    #'L-Proline', 'C00148', 'MNXM114'
    elif each_substrate == 'pro':
        met_name = 'pro_DASH_L'

    #'L-Serine', 'C00065', 'MNXM53'
    elif each_substrate == 'ser':
        met_name = 'ser_DASH_L'

    #'L-Threonine', 'C00188', 'MNXM142'
    elif each_substrate == 'thr':
        met_name = 'thr_DASH_L'

    #'L-Tryptophan', 'C00078', 'MNXM94'
    elif each_substrate == 'trp':
        met_name = 'trp_DASH_L'

    #'L-Tyrosine', 'C00082', 'MNXM76'
    elif each_substrate == 'tyr':
        met_name ='tyr_DASH_L'

    #'L-Valine', 'C00183', 'MNXM199'
    elif each_substrate == 'val':
        met_name = 'val_DASH_L'

    #'L-Isoleucine', 'C00407', 'MNXM231'
    elif each_substrate == 'ile':
        met_name = 'ile_DASH_L'

    #'phenylglycine', 'C18623', 'MNXM59292'
    elif each_substrate == 'phg':
        met_name = 'phg_DASH_L'

    #No MNXM, KEGG ID and bigg ID for "bht"
    #'beta-hydroxyl-tyrosine', 'N/A', 'N/A'
    elif each_substrate == 'bht':
        met_name = 'bht_DASH_L'

    #'Ornithine', 'C01602', 'MNXM89689'
    elif each_substrate == 'orn':
        met_name = 'orn'

    #No bigg ID
    #'abu', 'D-2-Aminobutyric acid', 'C02261', 'MNXM155019'
    elif each_substrate == 'abu':
        met_name = 'MNXM155019'

    #No bigg ID
    #'iva', '2-Amino-2-methylbutanoate', 'C03571', 'MNXM34821'
    elif each_substrate == 'iva':
        met_name = 'MNXM34821'

    #'L-2-Aminoadipic acid', 'C00956', 'MNXM268'
    elif each_substrate == 'aad':
        met_name = 'L2aadp'

    #No bigg ID
    #'hpg', 'D-4-Hydroxyphenylglycine', 'C03493', 'MNXM4544'
    elif each_substrate == 'hpg':
        met_name = 'MNXM4544'

    #'2,3-Dihydroxybenzoic acid', 'C00196', 'MNXM455'
    elif each_substrate == 'dhb':
        met_name = '23dhb'

    #No bigg ID
    #'dhpg', '3,5-Dihydroxy-phenylglycine', 'C12026', 'MNXM9962'
    elif each_substrate == 'dhpg':
        met_name = 'MNXM9962'

    #No bigg ID
    #'hty', 'L-Homotyrosine', 'C18622', 'MNXM59438'
    elif each_substrate == 'hty':
        met_name = 'MNXM59438'

    #'L-citruline', 'C00327', 'MNXM211'
    elif each_substrate == 'cit':
        met_name = 'citr_DASH_L'

    #'L-pipecolate', 'C00408', 'MNXM684'
    elif each_substrate == 'pip':
        met_name = 'Lpipecol'

    #'beta-alanine zwitterion', 'C00099', 'MNXM144'
    elif each_substrate == 'b-ala':
        met_name = 'ala_DASH_B'

    #'L-2,4-diazaniumylbutyrate', 'C03283', 'MNXM840'
    elif each_substrate == 'dab':
        met_name = '24dab'

    #'phenylacetate', 'C00548', 'MNXM497'
    elif each_substrate == 'phenylacetate' or each_substrate == 'Pha':
        met_name = 'pac'

    #'3-aminoalanine zwitterion', 'C06393', 'MNXM91374'
    elif each_substrate == '2-3-diaminoproprionate':
        met_name = '23dappa'

    #'tcl', '4-Chlorothreonine', 'N/A', 'MNXM37380'
    elif each_substrate == 'thr-4-cl':
        met_name = 'MNXM37380'

    #No bigg ID
    #(4S)-5,5,5-trichloro-leucine', 'N/A','N/A'
    elif each_substrate == 'tcl':
        met_name = 'MNXM37380'

    #No bigg ID
    #'qa', 'quinoxaline', 'C18575','MNXM80505'
    elif each_substrate == 'qa':
        met_name = 'MNXM80505'

    #No bigg ID
    #'Trans-cyclopentane-(1R, 2R)-dicarboxylic acid', 'N/A','N/A'
    elif each_substrate == 'trans-1,2-CPDA':
        met_name = '23cpda'

    #No bigg ID
    #'2-Butenyl-4-methyl-threonine', 'C12029','MNXM31446'
    elif each_substrate == 'bmt':
        met_name = 'MNXM31446'

    #'salicylate', 'C00805','MNXM378'
    elif each_substrate == 'sal':
        met_name = 'salc'

    #No bigg ID
    #'L-alaninol', 'N/A','MNXM8817'
    elif each_substrate == 'alaninol':
        met_name = 'MNXM8817'

    #t1pks substreate
    #'malonyl-CoA', 'C00083', 'MNXM40'
    elif each_substrate == 'mal':
        met_name = 'malcoa'

    #'(S)-methylmalonyl-CoA(5-)','C00683', 'MNXM190'
    elif each_substrate == 'mmal':
        met_name = 'mmcoa_DASH_S'

    #'2-methylbutanoyl-CoA', C01033,'MNXM569'
    elif each_substrate == '2metbut':
        met_name = '2mbcoa'

    #Not available in bigg, but available in iMK1208
    #'ethylmalonyl-CoA','C18026', 'MNXM2043'
    elif each_substrate == 'Ethyl_mal' or each_substrate == 'emal':
        met_name = 'emcoa_DASH_S'

    #'2-Methylpropanoyl-CoA', 'C00630', 'MNXM470'
    elif each_substrate == 'isobut':
        met_name = 'ibcoa'

    #'Acetyl-CoA', 'C00024', 'MNXM21'
    elif each_substrate == 'ace':
        met_name = 'accoa'

    #'Propionyl-CoA', 'C00100', 'MNXM86'
    elif each_substrate == 'prop':
        met_name = 'ppcoa'

    #'3-Methylbutanoyl-CoA', 'C02939', 'MNXM471'
    elif each_substrate == '3metbut':
        met_name = 'ivcoa'

    #No bigg ID
    #'mxmalacp', 'Methoxymalonyl-[acp]', 'C18616', 'MNXM61686'
    elif each_substrate == 'mxmal':
        met_name = 'MNXM61686'

    #No bigg ID
    #'chccoa', 'cyclohexane-1-carboxyl-CoA', 'C09823', 'MNXM5111'
    elif each_substrate == 'CHC-CoA':
        met_name = 'MNXM5111'

    #'chloroethylmalonyl-CoA','N/A', 'MNXM10927'
    elif each_substrate == 'cemal':
        met_name = 'MNXM10927'

    #Taken from Supplementary Table of Minowa et al.
    #'Quinaldinic acid','C06325', 'MNXM4797'
    elif each_substrate == 'qna':
        met_name = 'MNXM4797'

    #No bigg ID
    #'Benzoyl-CoA','C00512', 'MNXM240'
    elif each_substrate == 'benz':
        met_name = 'MNXM240'

    #No bigg ID
    #'L-Capreomycidine', 'C18472', 'MNXM18891'
    elif each_substrate == 'capreomycidine':
        met_name = 'MNXM18891'

    return met_name


#Add metabolite MNXM having no bigg ID to the model
def add_sec_met_mnxm_having_no_biggid_to_model(
        metab, metab_compt, mnxm_compoundInfo_dict):

    if metab == 'phg_DASH_L':
        metab2 = options.mnxref.metabolites.get_by_id('_'.join('MNXM59292', 'c'))
        metab2.id = metab_compt

    #No MNXM, KEGG ID and bigg ID for "bht"
    elif metab == 'bht_DASH_L':
        metab_compt = Metabolite(metab_compt, compartment='c')

    elif metab == 'MNXM155019': #'abu'
        metab_compt = Metabolite(metab_compt,
                formula = mnxm_compoundInfo_dict['MNXM155019'][1],
                name = mnxm_compoundInfo_dict['MNXM155019'][0],
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
