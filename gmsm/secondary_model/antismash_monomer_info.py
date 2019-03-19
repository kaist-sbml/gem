
import cobra

# Metabolite ID priority: bigg > MNXM
# KEGG IDs are not used in the modeling
# Abbreviations in antiSMASH gbk file are used if standard IDs are not available
def get_std_id_from_antismash_id(each_substrate):

    # ===================================================
    # NOTE: BiGG IDs available both in MNXref and iMK1208
    # ===================================================

    met_name = ''

    # Acetyl-CoA: MNXM21; C00024
    if each_substrate == 'ace':
        met_name = 'accoa'

    # L-alanine: MNXM32; C00041
    elif each_substrate == 'ala':
        met_name = 'ala_DASH_L'

    # NOTE: 'ala_B' without '_DASH_' in MNXref
    # beta-Alanine; MNXM144; C00099
    elif each_substrate == 'b-ala':
        met_name = 'ala_DASH_B'

    # L-Arginine: MNXM70; C00062
    elif each_substrate == 'arg':
        met_name = 'arg_DASH_L'

    # L-asparagine: MNXM147; C00152
    elif each_substrate == 'asn':
        met_name ='asn_DASH_L'

    # L-Aspartate: MNXM42; C00049
    elif each_substrate == 'asp':
        met_name = 'asp_DASH_L'

    # L-Citrulline; MNXM211; C00327
    elif each_substrate == 'cit':
        met_name = 'citr_DASH_L'

    # L-Cysteine: MNXM55; C00097'
    elif each_substrate == 'cys':
        met_name = 'cys_DASH_L'

    # L-Glutamine: MNXM37; C00064
    elif each_substrate == 'gln':
        met_name = 'gln_DASH_L'

    # L-Glutamate: MNXM89557; C00025
    elif each_substrate == 'glu':
        met_name = 'glu_DASH_L'

    # glycine: MNXM29; C00037
    elif each_substrate == 'gly':
        met_name = 'gly'

    # L-Histidine: MNXM134; C00135
    elif each_substrate == 'his':
        met_name = 'his_DASH_L'

    # L-Isoleucine: MNXM231; C00407
    elif each_substrate == 'ile':
        met_name = 'ile_DASH_L'

    # Isobutyryl-CoA: MNXM470; C00630
    elif each_substrate == 'isobut':
        met_name = 'ibcoa'

    # L-Leucine: MNXM140; C00123
    elif each_substrate == 'leu':
        met_name = 'leu_DASH_L'

    # L-Lysine: MNXM78; C00047
    elif each_substrate == 'lys':
        met_name = 'lys_DASH_L'

    # Malonyl-CoA: MNXM40; C00083
    elif each_substrate in ['mal', 'ohmal', 'ccmal', 'redmal', 'Malonyl-CoA']:
        met_name = 'malcoa'

    # L-Methionine: MNXM61; C00073
    elif each_substrate == 'met':
        met_name = 'met_DASH_L'

    # 2-Methylbutanoyl-CoA: MNXM569; C01033
    elif each_substrate == '2metbut':
        met_name = '2mbcoa'

    # 3-Methylbutanoyl-CoA or Isovaleryl-CoA: MNXM471; C02939
    elif each_substrate == '3metbut':
        met_name = 'ivcoa'

    # (R)-Methylmalonyl-CoA: MNXM608; C01213
    elif each_substrate in ['mmal', 'ohmmal', 'ccmmal', 'redmmal', 'Methylmalonyl-CoA']:
        met_name = 'mmcoa_DASH_R'

    # Ornithine: MNXM89689; C01602
    elif each_substrate == 'orn':
        met_name = 'orn'

    # Phenylacetic acid: MNXM497; C00548
    elif each_substrate in ['phenylacetate', 'Pha', 'phe-ac']:
        met_name = 'pac'

    # L-Phenylalanine: MNXM97; C00079
    elif each_substrate == 'phe':
        met_name = 'phe_DASH_L'

    # L-Proline: MNXM114; C00148
    elif each_substrate == 'pro':
        met_name = 'pro_DASH_L'

    # Propanoyl-CoA: MNXM86; C00100
    elif each_substrate == 'prop':
        met_name = 'ppcoa'

    # L-Serine: MNXM53; C00065
    elif each_substrate == 'ser':
        met_name = 'ser_DASH_L'

    # L-Threonine: MNXM142; C00188
    elif each_substrate == 'thr':
        met_name = 'thr_DASH_L'

    # L-Tryptophan: MNXM94; C00078
    elif each_substrate == 'trp':
        met_name = 'trp_DASH_L'

    # L-Tyrosine: MNXM76; C00082
    elif each_substrate == 'tyr':
        met_name ='tyr_DASH_L'

    # L-Valine: MNXM199; C00183
    elif each_substrate == 'val':
        met_name = 'val_DASH_L'

    # ====================================================
    # NOTE: BiGG IDs available only in MNXref, not iMK1208
    # ====================================================

    # L-2-Aminoadipic acid: MNXM268; C00956
    elif each_substrate == 'aad':
        met_name = 'L2aadp'

    # Benzoyl-CoA: MNXM240; C00512
    elif each_substrate == 'benz':
        met_name = 'benzcoa'

    # L-2,4-Diaminobutanoate: MNXM840; C03283
    elif each_substrate == 'dab':
        met_name = '24dab'

    # 2,3-Dihydroxybenzoic acid: MNXM455; C00196
    elif each_substrate == 'dhb':
        met_name = '23dhb'

    # FIXME: Is this wrong spelling used in antismash?
    # 2,3-diaminopropionate: MNXM91374; C06393
    elif each_substrate == '2-3-diaminoproprionate':
        met_name = '23dappa'

    # L-Pipecolate; MNXM684; C00408
    elif each_substrate == 'pip':
        met_name = 'Lpipecol'

    # salicylate: MNXM378; C00805
    elif each_substrate == 'sal':
        met_name = 'salc'

    # ================================================================
    # NOTE: BiGG IDs not available in MNXref, but adopted from iMK1208
    # ================================================================

    # (2S)-Ethylmalonyl-CoA: MNXM2043; C18026
    elif each_substrate in ['Ethyl_mal', 'emal', 'ohemal', 'ccemal', 'redemal']:
        met_name = 'emcoa_DASH_S'

    # =======================================================
    # NOTE: BiGG IDs not available both in MNXref and iMK1208
    # =======================================================

    # D-2-Aminobutyric acid: MNXM155019; C02261; no bigg id
    elif each_substrate == 'abu':
        met_name = 'MNXM155019'

    # L-alaninol: MNXM8817; no bigg id; no kegg id
    elif each_substrate == 'alaninol':
        met_name = 'MNXM8817'

    # NOTE: bmt is originally defined as `4-Butenyl-4-methylthreonine'
    # 2-Butenyl-4-methylthreonine: MNXM31446; C12029; no bigg id
    elif each_substrate == 'bmt':
        met_name = 'MNXM31446'

    # L-Capreomycidine: MNXM18891, C18472; no bigg id
    elif each_substrate == 'capreomycidine':
        met_name = 'MNXM18891'

    # chloroethylmalonyl-CoA: MNXM10927; no bigg id; no kegg id
    elif each_substrate == 'cemal':
        met_name = 'MNXM10927'

    # Cyclohexane-1-carboxyl-CoA or cyclohexylcarbonyl-CoA: MNXM5111; C09823; no bigg id
    elif each_substrate == 'CHC-CoA':
        met_name = 'MNXM5111'

    # 3,5-Dihydroxy-phenylglycine; MNXM9962; C12026; no bigg id
    elif each_substrate == 'dhpg':
        met_name = 'MNXM9962'

    # D-4-Hydroxyphenylglycine: MNXM4544; C03493; no bigg id
    elif each_substrate == 'hpg':
        met_name = 'MNXM4544'

    # L-Homotyrosine; MNXM59438; C18622; no bigg id
    elif each_substrate == 'hty':
        met_name = 'MNXM59438'

    # NOTE: iva is originally defined as `isovaline'
    # NOTE: 2-Amino-2-methylbutanoate is structurally similar to valine,
    #but a reason for the selection of this compound is not clear
    # 2-Amino-2-methylbutanoate: MNXM34821; C03571; no bigg id
    elif each_substrate == 'iva':
        met_name = 'MNXM34821'

    # (2R)-Methoxymalonyl-[acp]: MNXM61686; C18616; no bigg id
    elif each_substrate in ['mxmal', 'ohmxmal', 'ccmxmal', 'redmxmal']:
        met_name = 'MNXM61686'

    # phenylglycine: MNXM59292; C18623
    elif each_substrate == 'phg':
        met_name = 'MNXM59292'

    # Quinoxaline: MNXM80501; C18575; no bigg id
    elif each_substrate == 'qa':
        met_name = 'MNXM80501'

    # NOTE: Taken from Supplementary Table of Minowa et al.
    # Quinaldinic acid: MNXM4797; C06325
    elif each_substrate == 'qna':
        met_name = 'MNXM4797'

    # 4-Chlorothreonine or tcl: MNXM37380; no bigg id; no kegg id
    elif each_substrate == 'thr-4-cl' or each_substrate == 'tcl':
        met_name = 'MNXM37380'

    # ==========================================================
    # NOTE: No standard IDs available both in MNXref and iMK1208
    # ==========================================================

    # beta-hydroxyl-tyrosine: no mnxm id; no bigg id; no kegg id
    elif each_substrate == 'bht':
        met_name = 'bht'

    # Trans-cyclopentane-(1R, 2R)-dicarboxylic acid; no mnxm id; no bigg id; no kegg id
    elif each_substrate == 'trans-1,2-CPDA':
        met_name = '12cpda'

    if met_name:
        met_name = cobra.io.sbml.fix_legacy_id(met_name)
        return met_name
    else:
        return

