
import cobra

# Metabolite ID priority: bigg > MNXM
# KEGG IDs are not used in the modeling
# Abbreviations in antiSMASH gbk file are used if standard IDs are not available
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
    elif each_substrate in ['ala', 'ala-d']:
        met_name = 'ala__L'

    # beta-Alanine; MNXM144; C00099
    elif each_substrate in ['ala-b', 'b-ala', 'beta-ala']:
        met_name = 'ala__B'

    # L-Arginine: MNXM739527; C00062
    elif each_substrate == 'arg':
        met_name = 'arg__L'

    # L-asparagine: MNXM728945; C00152
    elif each_substrate == 'asn':
        met_name ='asn__L'

    # L-Aspartate: MNXM735836; C00049
    elif each_substrate == 'asp':
        met_name = 'asp__L'

    # L-Citrulline; MNXM732471; C00327
    elif each_substrate == 'cit':
        met_name = 'citr__L'

    # L-Cysteine: MNXM738068; C00097
    elif each_substrate == 'cys':
        met_name = 'cys__L'

    # L-Glutamine: MNXM37; C00064
    elif each_substrate == 'gln':
        met_name = 'gln__L'

    # L-Glutamate: MNXM741173; C00025
    elif each_substrate == 'glu':
        met_name = 'glu__L'

    # glycine: MNXM29; C00037
    elif each_substrate == 'gly':
        met_name = 'gly'

    # L-Histidine: MNXM728101; C00135
    elif each_substrate == 'his':
        met_name = 'his__L'

    # L-Isoleucine: MNXM728337; C00407
    elif each_substrate == 'ile':
        met_name = 'ile__L'

    # Isobutyryl-CoA: MNXM736714; C00630
    elif each_substrate == 'isobut':
        met_name = 'ibcoa'

    # L-Leucine: MNXM140; C00123
    elif each_substrate == 'leu':
        met_name = 'leu__L'

    # L-Lysine: MNXM78; C00047
    elif each_substrate == 'lys':
        met_name = 'lys__L'

    # Malonyl-CoA: MNXM735439; C00083
    elif each_substrate in ['mal', 'ohmal', 'ccmal', 'redmal']:
        met_name = 'malcoa'

    # L-Methionine: MNXM738804; C00073
    elif each_substrate == 'met':
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
    elif each_substrate == 'orn':
        met_name = 'orn'

    # Phenylacetic acid: MNXM737265; C07086
    elif each_substrate in ['phenylacetate', 'Pha', 'phe-ac']:
        met_name = 'pac'

    # L-Phenylalanine: MNXM741664; C00079
    elif each_substrate == 'phe':
        met_name = 'phe__L'

    # L-Proline: MNXM114; C00148
    elif each_substrate == 'pro':
        met_name = 'pro__L'

    # Propanoyl-CoA: MNXM740852; C00100
    elif each_substrate == 'prop':
        met_name = 'ppcoa'

    # L-Serine: MNXM737787; C00065
    elif each_substrate == 'ser':
        met_name = 'ser__L'

    # L-Threonine: MNXM142; C00188
    elif each_substrate == 'thr':
        met_name = 'thr__L'

    # L-Tryptophan: MNXM741554; C00078
    elif each_substrate == 'trp':
        met_name = 'trp__L'

    # L-Tyrosine: MNXM76; C00082
    elif each_substrate == 'tyr':
        met_name ='tyr__L'

    # L-Valine: MNXM199; C00183
    elif each_substrate == 'val':
        met_name = 'val__L'

    # ====================================================
    # NOTE: BiGG IDs available only in MNXref, not iKS1317
    # ====================================================

    # L-2-Aminoadipic acid: MNXM268; C00956
    elif each_substrate == 'aad':
        met_name = 'L2aadp'

    # Benzoyl-CoA: MNXM732896; C00512
    elif each_substrate == 'benz':
        met_name = 'benzcoa'

    # L-2,4-Diaminobutanoate: MNXM840; C03283
    elif each_substrate == 'dab':
        met_name = '24dab'

    # 2,3-Dihydroxybenzoic acid: MNXM455; C00196
    elif each_substrate == 'dhb':
        met_name = '23dhb'

    # Beta-lysine: MNXM1865, C01142
    elif each_substrate == 'lys-b':
        met_name = '36dahx'

    # D-2-hydroxyisovalerate: no mnxm id; no kegg id
    elif each_substrate in ['hiv', 'hiv-d', 'hyv-d']:
        met_name = '2hiv'

    # FIXME: Is this wrong spelling used in antismash?
    # 2,3-diaminopropionate: MNXM741395; C06393
    elif each_substrate == '2-3-diaminoproprionate':
        met_name = '23dappa'

    # L-Pipecolate; MNXM684; C00408
    elif each_substrate == 'pip':
        met_name = 'Lpipecol'

    # salicylate: MNXM730078; C00805
    elif each_substrate == 'sal':
        met_name = 'salc'

    # 2-hydroxyisocaproate: MNXM733822; C03264
    elif each_substrate == 'alpha-hydroxy-isocaproic-acid':
        met_name = 'ahisocap'

    # (S)-2-Aminobutanoate: MNXM17054; C02356
    elif each_substrate == 'abu':
        met_name = 'abu'

    # ================================================================
    # NOTE: BiGG IDs not available in MNXref, but adopted from iKS1317
    # ================================================================

    # (2S)-Ethylmalonyl-CoA: MNXM732335; C18026
    elif each_substrate in ['Ethyl_mal', 'emal', 'ohemal', 'ccemal', 'redemal']:
        met_name = 'emcoa__S'

    # =======================================================
    # NOTE: BiGG IDs not available both in MNXref and iKS1317
    # =======================================================

    # L-alaninol: MNXM8817; no bigg id; no kegg id
    elif each_substrate == 'alaninol':
        met_name = 'MNXM8817'

    # NOTE: bmt is originally defined as `4-Butenyl-4-methylthreonine'
    # 2-Butenyl-4-methylthreonine: MNXM31446; C12029; no bigg id
    elif each_substrate == 'bmt':
        met_name = 'MNXM31446'

    # L-Capreomycidine: MNXM30936, C18472; no bigg id
    elif each_substrate in ['cap', 'capreomycidine']:
        met_name = 'cap'

    # chloroethylmalonyl-CoA: MNXM10927; no bigg id; no kegg id
    elif each_substrate == 'cemal':
        met_name = 'MNXM10927'

    # Cyclohexane-1-carboxyl-CoA or cyclohexylcarbonyl-CoA: MNXM5111; C09823; no bigg id
    elif each_substrate == 'CHC-CoA':
        met_name = 'chccoa'

    # 3,5-Dihydroxy-phenylglycine; MNXM9962; C12026; no bigg id
    elif each_substrate == 'dhpg':
        met_name = 'dhpg'

    # L-4-Hydroxyphenylglycine: MNXM12045; C03493; no bigg id
    elif each_substrate == 'hpg':
        met_name = 'hpg'

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
        met_name = 'mxmal'

    # Phenylglycine: MNXM729975; C18623; no bigg id
    elif each_substrate == 'phg':
        met_name = 'phg'

    # Quinoxaline: MNXM80501; C18575; no bigg id
    elif each_substrate == 'qa':
        met_name = 'MNXM80501'

    # NOTE: Taken from Supplementary Table of Minowa et al.
    # Quinaldinic acid: MNXM4797; C06325
    elif each_substrate == 'qna':
        met_name = 'MNXM4797'

    # 4-Chlorothreonine or tcl: MNXM37380; no bigg id; no kegg id
    elif each_substrate == 'thr-4-cl':
        met_name = 'MNXM37380'

    # N5-acetyl-N5-hydroxy-L-ornithine: MNXM19786; no bigg id; no kegg id
    elif each_substrate == 'haorn':
        met_name = 'haorn'

    # Beta-hydroxyl-tyrosine: MNXM149582; no bigg id; no kegg id
    elif each_substrate == 'bht':
        met_name = 'bht'

    # ==========================================================
    # NOTE: No standard IDs available both in MNXref and iKS1317
    # ==========================================================

    # Trans-cyclopentane-(1R, 2R)-dicarboxylic acid; no mnxm id; no bigg id; no kegg id
    elif each_substrate == 'trans-1,2-CPDA':
        met_name = '12cpda'

    # (4S)-5,5,5-Trichloro-leucine: no mnxm id; no bigg id; no kegg id
    elif each_substrate == 'tcl':
        met_name = 'tcl'

    # L-3-Methylglutamate: no mnxm id, no bigg id, no kegg id
    elif each_substrate == '3-me-glu':
        met_name = '3meglu'

    # 1-(1,1-Dimethyl-2-propenyl)tryptophan: no mnxm id, no bigg id, no kegg id
    elif each_substrate == 'N-(1,1-dimethyl-1-allyl)Trp':
        met_name = 'dmalltrp'

    if met_name:
        #met_name = cobra.io.sbml.fix_legacy_id(met_name)
        return met_name
    else:
        return

