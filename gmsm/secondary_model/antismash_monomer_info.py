
import cobra

# Metabolite ID priority: bigg > MNXM
# KEGG IDs are not used in the modeling
# Abbreviations in antiSMASH gbk file are used if standard IDs are not available
# Refer to:
# https://github.com/antismash/antismash/blob/master/antismash/modules/nrps_pks/name_mappings.py
# https://github.com/antismash/antismash/blob/master/antismash/modules/nrps_pks/data/aaSMILES.txt
# Some metabolites are named by our own

def get_std_id_from_antismash_id(each_substrate):

    # ===================================================
    # NOTE: BiGG IDs available both in MNXref and iKS1317
    # ===================================================

    met_name = ''

    # Acetyl-CoA: MNXM732448; C00024
    if each_substrate == 'ace':
        met_name = 'accoa'

    # L-alanine: MNXM733982; C00041
    elif each_substrate in ['ala', 'ala-d', 'Ala']:
        met_name = 'ala__L'

    # beta-Alanine; MNXM144; C00099
    elif each_substrate in ['ala-b', 'b-ala', 'beta-ala', 'bAla']:
        met_name = 'ala__B'

    # L-Arginine: MNXM739527; C00062
    elif each_substrate in ['arg', 'Arg']:
        met_name = 'arg__L'

    # L-asparagine: MNXM728945; C00152
    elif each_substrate in ['asn', 'Asn']:
        met_name ='asn__L'

    # L-Aspartate: MNXM735836; C00049
    elif each_substrate in ['asp', 'Asp']:
        met_name = 'asp__L'

    # L-Citrulline; MNXM732471; C00327
    elif each_substrate in ['cit', 'Cit']:
        met_name = 'citr__L'

    # L-Cysteine: MNXM738068; C00097
    elif each_substrate in ['cys', 'Cys']:
        met_name = 'cys__L'

    # L-Glutamine: MNXM37; C00064
    elif each_substrate in ['gln', 'Gln']:
        met_name = 'gln__L'

    # L-Glutamate: MNXM741173; C00025
    elif each_substrate in ['glu', 'Glu']:
        met_name = 'glu__L'

    # glycine: MNXM29; C00037
    elif each_substrate in ['gly', 'Gly']:
        met_name = 'gly'

    # L-Histidine: MNXM728101; C00135
    elif each_substrate in ['his', 'His']:
        met_name = 'his__L'

    # L-Isoleucine: MNXM728337; C00407
    elif each_substrate in ['ile', 'Ile']:
        met_name = 'ile__L'

    # Isobutyryl-CoA: MNXM736714; C00630
    elif each_substrate == 'isobut':
        met_name = 'ibcoa'

    # L-Leucine: MNXM140; C00123
    elif each_substrate in ['leu', 'Leu']:
        met_name = 'leu__L'

    # L-Lysine: MNXM78; C00047
    elif each_substrate in ['lys', 'Lys']:
        met_name = 'lys__L'

    # Malonyl-CoA: MNXM735439; C00083
    elif each_substrate in ['mal', 'ohmal', 'ccmal', 'redmal']:
        met_name = 'malcoa'

    # L-Methionine: MNXM738804; C00073
    elif each_substrate in ['met', 'Met']:
        met_name = 'met__L'

    # 2-Methylbutanoyl-CoA: MNXM730317; C01033
    elif each_substrate == '2metbut':
        met_name = '2mbcoa'

    # 3-Methylbutanoyl-CoA or Isovaleryl-CoA: MNXM736744; C02939
    elif each_substrate == '3metbut':
        met_name = 'ivcoa'

    # (R)-Methylmalonyl-CoA: MNXM738783; C01213
    elif each_substrate in ['mmal', 'ohmmal', 'ccmmal', 'redmmal']:
        met_name = 'mmcoa__R'

    # Ornithine: MNXM741175; C00077
    elif each_substrate in ['orn', 'Orn']:
        met_name = 'orn'

    # Phenylacetic acid: MNXM737265; C07086
    elif each_substrate in ['phenylacetate', 'Pha', 'phe-ac']:
        met_name = 'pac'

    # L-Phenylalanine: MNXM741664; C00079
    elif each_substrate in ['phe', 'Phe']:
        met_name = 'phe__L'

    # L-Proline: MNXM114; C00148
    elif each_substrate in ['pro', 'Pro']:
        met_name = 'pro__L'

    # Propanoyl-CoA: MNXM740852; C00100
    elif each_substrate == 'prop':
        met_name = 'ppcoa'

    # L-Serine: MNXM737787; C00065
    elif each_substrate in ['ser', 'Ser']:
        met_name = 'ser__L'

    # L-Threonine: MNXM142; C00188
    elif each_substrate in ['thr', 'Thr']:
        met_name = 'thr__L'

    # L-Tryptophan: MNXM741554; C00078
    elif each_substrate in ['trp', 'Trp']:
        met_name = 'trp__L'

    # L-Tyrosine: MNXM76; C00082
    elif each_substrate in ['tyr', 'Tyr']:
        met_name ='tyr__L'

    # L-Valine: MNXM199; C00183
    elif each_substrate in ['val', 'Val']:
        met_name = 'val__L'
        
    # L-allo-threonine: MNXM2125; C05519
    elif each_substrate == 'aThr':
        met_name = 'athr__L'
        
    # D-Alanine: MNXM156; C00133
    elif each_substrate == 'D-Ala':
        met_name = 'ala__D'
        
    # LL-2,6-Diaminoheptanedioate: MNXM644; C00666
    elif each_substrate == 'LDAP':
        met_name = '26dap_LL'
        
    # 3-Methyl-2-oxopentanoate: MNXM1103447; C00671
    elif each_substrate == 'MeHOval':
        met_name = '3mop'

    # 3-methyl-2-oxobutanoate: MNXM732866; C00141
    elif each_substrate == '2-oxo-isovaleric-acid':
        met_name = '3mob'
        
    # salicylate: MNXM730078; C00805
    elif each_substrate in ['sal', 'Sal']:
        met_name = 'salc'
        
    # ====================================================
    # NOTE: BiGG IDs available only in MNXref, not iKS1317
    # ====================================================

    # L-2-Aminoadipic acid: MNXM268; C00956
    elif each_substrate in ['aad', 'Aad']:
        met_name = 'L2aadp'

    # Benzoyl-CoA: MNXM732896; C00512
    elif each_substrate == 'benz':
        met_name = 'benzcoa'

    # L-2,4-Diaminobutanoate: MNXM840; C03283
    elif each_substrate in ['dab', 'Dab']:
        met_name = '24dab'

    # 2,3-Dihydroxybenzoic acid: MNXM455; C00196
    elif each_substrate in ['dhb', 'diOH-Bz']:
        met_name = '23dhb'

    # Beta-lysine: MNXM1865, C01142
    elif each_substrate in ['lys-b', 'bLys']:
        met_name = '36dahx'

    # D-2-hydroxyisovalerate: MNXM1106847; no kegg id
    elif each_substrate in ['hiv', 'hiv-d', 'D-Hiv']:
        met_name = '2hiv'

    # FIXME: Is this wrong spelling used in antismash?
    # 2,3-diaminopropionate: MNXM741395; C06393
    elif each_substrate in ['2-3-diaminoproprionate', 'Dpr']:
        met_name = '23dappa'

    # L-Pipecolate; MNXM684; C00408
    elif each_substrate in ['pip', 'Hpr']:
        met_name = 'Lpipecol'
        
    # N5-hydroxy-L-ornithine: MNXM6354; C20850
    elif each_substrate == 'OH-Orn':
        met_name = 'n5horn'

    # =======================================================
    # NOTE: BiGG IDs not available both in MNXref and iKS1317
    # =======================================================

    # L-alaninol: MNXM8817; no kegg id
    elif each_substrate in ['alaninol', 'Aol']:
        met_name = 'MNXM8817'

    # L-Capreomycidine: MNXM30936; C18472
    elif each_substrate in ['cap', 'capreomycidine', 'Cap']:
        met_name = 'cap'

    # chloroethylmalonyl-CoA: MNXM10927; no kegg id
    elif each_substrate == 'cemal':
        met_name = 'MNXM10927'

    # Cyclohexane-1-carboxyl-CoA or cyclohexylcarbonyl-CoA: MNXM5111; C09823
    elif each_substrate == 'CHC-CoA':
        met_name = 'chccoa'

    # 3,5-Dihydroxy-phenylglycine; MNXM9962; C12026
    elif each_substrate in ['dhpg', 'Dhpg']:
        met_name = 'dhpg'

    # L-4-Hydroxyphenylglycine: MNXM12045; C03493
    elif each_substrate in ['hpg', 'Hpg']:
        met_name = 'hpg'

    # L-Homotyrosine: MNXM59438; C18622
    elif each_substrate in ['hty', 'Hty']:
        met_name = 'hty'
        
    # 4-propyl-L-proline: MNXM1104981; no kegg id
    elif each_substrate == 'pPro':
        met_name = 'ppro'
        
    # 2,3-dehydroaminobutyric acid/(Z)-2-aminobutenoic acid: MNXM1454; C17234  
    elif each_substrate in ['dht', 'dhAbu', 'Dht']:
        met_name = 'dht'
        
    # 2-hydroxyvalerate/2-hydroxypentanoate: MNXM14267; no kegg id
    elif each_substrate == 'hyv-d':
        met_name = '2hv'
        
    # Lysergic acid: MNXM738781; C07541
    elif each_substrate == 'd-lyserg':
        met_name = 'lyserg'
        
    # (2S,3S)-3-methylphenylalanine: MNXM146324; C20895
    elif each_substrate in ['mephe', 'mePhe']:
        met_name = 'mephe'
        
    # Hydroxyasparagine: MNXM1127841; no kegg id
    elif each_substrate in ['hasn', 'hAsn']:
        met_name = 'hasn'
    
    # 9-N-methoxy-tryptophan: MNXM147393; no kegg id
    elif each_substrate == 's-nmethoxy-trp':
        met_name = 'mxtrp'
        
    # L-2-amino-8-oxodecanoate: MNXM148306; no kegg id
    elif each_substrate in ['aoda', 'n-oxoDec']:
        met_name = '2a8odca'
        
    # 2-hydroxyisocaproate: MNXM733822; C03264
    elif each_substrate in ['alpha-hydroxy-isocaproic-acid', '2S-Hic']:
        met_name = 'ahisocap'
            
    # L-allo-isoleucine, MNXM59280; C21096
    elif each_substrate == 'aIle':
        met_name = 'aile__L'
        
    # (S)-2-Aminobutanoate: MNXM17054; C02356
    elif each_substrate in ['abu', 'Abu']:
        met_name = 'abu'
        
    # isovaline: MNXM34821; C03571
    elif each_substrate == 'IVal':
        met_name = 'iva'
        
    # (2R)-Methoxymalonyl-[acp]: MNXM61686; C18616
    elif each_substrate in ['mxmal', 'ohmxmal', 'ccmxmal', 'redmxmal']:
        met_name = 'mxmal'

    # Phenylglycine: MNXM729975; C18623
    elif each_substrate in ['phg', 'Pgl']:
        met_name = 'phg'

    # Quinoxaline: MNXM80501; C18575
    elif each_substrate == 'qa':
        met_name = 'MNXM80501'

    # NOTE: Taken from Supplementary Table of Minowa et al.
    # Quinaldinic acid: MNXM4797; C06325
    elif each_substrate == 'qna':
        met_name = 'MNXM4797'

    # 4-Chlorothreonine: MNXM37380; no kegg id
    elif each_substrate == 'thr-4-cl':
        met_name = 'MNXM37380'

    # N5-formyl-N5-hydroxyornithine: MNXM12420; no kegg id
    elif each_substrate == 'Fo-OH-Orn':
        met_name = 'hforn'
        
    # N5-acetyl-N5-hydroxy-L-ornithine: MNXM19786; no kegg id
    elif each_substrate in ['haorn', 'Ac-OH-Orn']:
        met_name = 'haorn'

    # Beta-hydroxyl-tyrosine: MNXM149582; no kegg id
    elif each_substrate in ['bht', 'bOH-Tyr']:
        met_name = 'bht'

    # ================================================================
    # NOTE: BiGG IDs not available in MNXref, but adopted from iKS1317
    # ================================================================

    # (2S)-Ethylmalonyl-CoA: MNXM732335; C18026
    elif each_substrate in ['Ethyl_mal', 'emal', 'ohemal', 'ccemal', 'redemal']:
        met_name = 'emcoa__S'

    # ==========================================================
    # NOTE: No standard IDs available both in MNXref and iKS1317
    # ==========================================================

    # Trans-cyclopentane-(1R, 2R)-dicarboxylic acid
    elif each_substrate == 'trans-1,2-CPDA':
        met_name = '12cpda'

    # (4S)-5,5,5-Trichloro-leucine
    elif each_substrate in ['tcl', '3clLeu']:
        met_name = 'tcl'

    # L-3-Methylglutamate
    elif each_substrate == '3Me-Glu':
        met_name = '3meglu'

    # 1-(1,1-Dimethyl-2-propenyl)tryptophan
    elif each_substrate in ['N-(1,1-dimethyl-1-allyl)Trp', 'dmeaTrp']:
        met_name = 'dmalltrp'
    
    # 3,5-dichloro-4-hydroxy-phenylglycine
    elif each_substrate == 'Cl2-Hpg':
        met_name = 'dclhpg'
    
    # 2-hydroxy-3-methyl-pentanoic acid
    elif each_substrate in ['D-Hmp', 'Hmp']:
        met_name = 'hmp__D'
        
    # 2-amino-9,10-epoxi-8-oxodecanoic acid
    elif each_substrate == 'C10:0-NH2(2)-Ep(9)-oxo(8)':
        met_name = '2a8odca'
        
    # valinol
    elif each_substrate == 'Valol':
        met_name = 'valol'
        
    # 4-Butenyl-4-methylthreonine
    elif each_substrate in ['bmt', 'Bmt']:
        met_name = 'bmt'
        
    # 4-hydroxy-L-valine/2-Amino-4-hydroxy-3-methylbutanoic acid
    elif each_substrate == 'hyv':
        met_name = 'hval__L'
        
    # 2-Amino-3,3-dimethylbutanoic acid/tert-leucine/3-methyl-L-valine
    elif each_substrate == 'meval':
        met_name = 'meval'

    if met_name:
        #met_name = cobra.io.sbml.fix_legacy_id(met_name)
        return met_name
    else:
        return

