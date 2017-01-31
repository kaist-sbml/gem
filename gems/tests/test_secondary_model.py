
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from argparse import Namespace
from os.path import join
from gems.secondary_model.sec_met_rxn_generation import (
        get_cluster_location,
        get_cluster_product,
        get_cluster_info_from_seq_record,
        get_cluster_domain,
        get_cluster_monomers,
        get_all_metab_coeff
        )

class TestSecondary_model:
    """Test functions in gems.secondary_model"""

    def test_get_cluster_info_from_seq_record(self, seq_record, options):

        options.seq_record = seq_record

        # Hybrid cluster: nrps-transatpks-t1pks
        # locations: 341017 - 503094
        # Kirromycin biosynthetic gene cluster (81% of genes show similarity)
        cluster_nr = 3

        get_cluster_location(cluster_nr, options)
        get_cluster_info_from_seq_record(options)

        # Number of genes for Cluster 3 in Streptomyces collinus Tu 365
        assert len(options.cluster_info_dict) == 22

        for locustag in options.cluster_info_dict.keys():
            assert 'NRPS/PKS Domain' in str(options.cluster_info_dict[locustag])


    def test_get_cluster_product(self, seq_record, options):

        options.seq_record = seq_record
        cluster_nr = 3

        get_cluster_location(cluster_nr, options)
        get_cluster_product(cluster_nr, options)

        assert options.product == 'Cluster03_nrps_t1pks_transatpks'


    def test_get_cluster_domain(self, seq_record, options):

        options.seq_record = seq_record
        # 'cluster_info_dict' is too long to write in file
        cluster_nr = 3

        get_cluster_location(cluster_nr, options)
        get_cluster_info_from_seq_record(options)

        get_cluster_domain(options)

        assert len(options.locustag_domain_dict['B446_01400']) == 1
        assert 'PKS_ER' in str(options.locustag_domain_dict['B446_01400'])

        assert len(options.locustag_domain_dict['B446_01405']) == 1
        assert 'Polyketide_cyc' in str(options.locustag_domain_dict['B446_01405'])

        assert len(options.locustag_domain_dict['B446_01480']) == 7
        assert 'Condensation_Starter_D00' in options.locustag_domain_dict['B446_01480'] \
            and 'AMP-binding_D01' in options.locustag_domain_dict['B446_01480'] \
            and 'PCP_D02' in options.locustag_domain_dict['B446_01480'] \
            and 'Epimerization_D03' in options.locustag_domain_dict['B446_01480'] \
            and 'Condensation_DCL_D04' in options.locustag_domain_dict['B446_01480'] \
            and 'AMP-binding_D05' in options.locustag_domain_dict['B446_01480'] \
            and 'PCP_D06' in options.locustag_domain_dict['B446_01480']

        assert len(options.locustag_domain_dict['B446_01485']) == 10
        assert 'Condensation_LCL_D00' in options.locustag_domain_dict['B446_01485'] \
            and 'AMP-binding_D01' in options.locustag_domain_dict['B446_01485'] \
            and 'PCP_D02' in options.locustag_domain_dict['B446_01485'] \
            and 'Condensation_LCL_D03' in options.locustag_domain_dict['B446_01485'] \
            and 'AMP-binding_D04' in options.locustag_domain_dict['B446_01485'] \
            and 'PCP_D05' in options.locustag_domain_dict['B446_01485'] \
            and 'Condensation_LCL_D06' in options.locustag_domain_dict['B446_01485'] \
            and 'AMP-binding_D07' in options.locustag_domain_dict['B446_01485'] \
            and 'PCP_D08' in options.locustag_domain_dict['B446_01485'] \
            and 'Thioesterase_D09' in options.locustag_domain_dict['B446_01485']

        assert len(options.locustag_domain_dict['B446_01525']) == 1
        assert len(options.locustag_domain_dict['B446_01530']) == 4
        assert len(options.locustag_domain_dict['B446_01535']) == 3
        assert len(options.locustag_domain_dict['B446_01540']) == 2
        assert len(options.locustag_domain_dict['B446_01565']) == 7
        assert len(options.locustag_domain_dict['B446_01590']) == 1
        assert len(options.locustag_domain_dict['B446_01635']) == 1
        assert len(options.locustag_domain_dict['B446_01640']) == 1
        assert len(options.locustag_domain_dict['B446_01655']) == 2
        assert len(options.locustag_domain_dict['B446_01660']) == 3
        assert len(options.locustag_domain_dict['B446_01670']) == 9
        assert len(options.locustag_domain_dict['B446_01675']) == 8
        assert len(options.locustag_domain_dict['B446_01680']) == 20
        assert 'PKS_DH2_D00' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KR_D01' in options.locustag_domain_dict['B446_01680'] \
            and 'ACP_D02' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KS_D03' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_DH2_D04' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KR_D05' in options.locustag_domain_dict['B446_01680'] \
            and 'cMT_D06' in options.locustag_domain_dict['B446_01680'] \
            and 'ACP_D07' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KS_D08' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_DHt_D09' in options.locustag_domain_dict['B446_01680'] \
            and 'ACP_D10' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KS_D11' in options.locustag_domain_dict['B446_01680'] \
            and 'Trans-AT_docking_D12' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KR_D13' in options.locustag_domain_dict['B446_01680'] \
            and 'ACP_D14' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KS_D15' in options.locustag_domain_dict['B446_01680'] \
            and 'Trans-AT_docking_D16' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KR_D17' in options.locustag_domain_dict['B446_01680'] \
            and 'ACP_D18' in options.locustag_domain_dict['B446_01680'] \
            and 'PKS_KS_D19' in options.locustag_domain_dict['B446_01680']

        assert len(options.locustag_domain_dict['B446_01685']) == 4
        assert len(options.locustag_domain_dict['B446_01690']) == 14
        assert len(options.locustag_domain_dict['B446_01695']) == 6
        assert len(options.locustag_domain_dict['B446_01700']) == 1
        assert len(options.locustag_domain_dict['B446_01730']) == 1


    def test_get_cluster_monomers(self, seq_record, options):

        options.seq_record = seq_record
        cluster_nr = 3

        get_cluster_location(cluster_nr, options)
        get_cluster_info_from_seq_record(options)
        get_cluster_monomers(options)

        assert len(options.locustag_monomer_dict.keys()) == 17

        #Consensus monomers
        assert options.locustag_monomer_dict['B446_01480_M0'][3] == 'nrp'
        assert options.locustag_monomer_dict['B446_01480_M1'][3] == 'ser'
        assert options.locustag_monomer_dict['B446_01670_M0'][2] == 'ccmmal'


    def test_get_all_metab_coeff(self, options):

        locustag_monomer_dict = {
                'B446_01480_M0': ['orn,lys,arg', 'lys', 'leu', 'nrp'],
                'B446_01480_M1': ['ser', 'ser', 'ser', 'ser'],
                'B446_01485_M0': ['val', 'val', 'val', 'val'],
                'B446_01485_M1': ['asp,asn,glu,gln,aad', 'N/A', 'ala', 'nrp'],
                'B446_01485_M2': ['gly', 'gly', 'gly', 'gly'],
                'B446_01525_M0': ['hydrophilic', 'arg', 'arg', 'arg'],
                'B446_01530_M0': ['N/A', 'ser', 'dab', 'nrp'],
                'B446_01535_M0': ['val', 'val', 'val', 'val'],
                'B446_01565_M0': ['thr', 'thr', 'thr', 'thr'],
                'B446_01565_M1': ['ser', 'ser', 'ser', 'ser'],
                'B446_01635_M0': ['mmal', 'mmal', 'mmal'],
                'B446_01655_M0': ['mal', 'mal', 'mal'],
                'B446_01655_M1': ['mal', 'mal', 'mal'],
                'B446_01660_M0': ['gly', 'gly', 'gly', 'gly'],
                'B446_01670_M0': ['mmal', 'mmal', 'ccmmal'],
                'B446_01670_M1': ['mal', 'mal', 'mal'],
                'B446_01685_M0': ['gly', 'gly', 'gly', 'gly']}

        options.locustag_monomer_dict = locustag_monomer_dict
        options.product = 'Cluster03_nrps_t1pks_transatpks'
        get_all_metab_coeff(options)

        assert len(options.metab_coeff_dict) == 71
        assert options.metab_coeff_dict['atp'] == 0
        assert options.metab_coeff_dict['mmcoa_DASH_S'] == -2
        assert options.metab_coeff_dict['ser_DASH_L'] == -2
        assert options.metab_coeff_dict['val_DASH_L'] == -2
        assert options.metab_coeff_dict['Cluster03_nrps_t1pks_transatpks'] == 1

#{'asn_DASH_L': 0, '2mbcoa': 0, 'ppcoa': 0, 'MNXM34821': 0, 'ppi': 0, 'MNXM10927': 0, 'hco3': 0, 'MNXM59438': 0, 'ile_DASH_L': 0, '24dab': -1, 'citr_DASH_L': 0, 'MNXM9962': 0, 'leu_DASH_L': -1, 'gln_DASH_L': 0, 'fmn': 0, 'mmcoa_DASH_S': -2, 'pro_DASH_L': 0, 'MNXM37380': 0, 'cys_DASH_L': 0, 'malcoa': -3, 'lys_DASH_L': 0, 'MNXM4797': 0, 'MNXM5111': 0, 'nadph': 0, 'MNXM17054': 0, 'pac': 0, 'MNXM80505': 0, 'MNXM61686': 0, 'arg_DASH_L': -1, 'MNXM18891': 0, 'Lpipecol': 0, 'bht_DASH_L': 0, 'accoa': 0, '23dappa': 0, 'ibcoa': 0, 'MNXM8817': 0, 'ser_DASH_L': -2, 'ahcys': 0, 'ivcoa': 0, 'orn': 0, 'thr_DASH_L': -1, 'MNXM31446': 0, 'Cluster03_nrps_t1pks_transatpks': 1, 'amp': 0, 'nadp': 0, 'L2aadp': 0, 'amet': 0, 'val_DASH_L': -2, 'glu_DASH_L': 0, 'fmnh2': 0, 'phg_DASH_L': 0, 'trp_DASH_L': 0, 'phe_DASH_L': 0, 'h': 0, 'emcoa_DASH_S': 0, 'salc': 0, 'tcl': 0, '23cpda': 0, 'coa': 0, 'gly': -3, 'h2o': 0, 'ala_DASH_L': -1, 'MNXM4544': 0, 'atp': 0, 'tyr_DASH_L': 0, 'ala_DASH_B': 0, 'asp_DASH_L': 0, '23dhb': 0, 'MNXM240': 0, 'met_DASH_L': 0, 'his_DASH_L': 0}
