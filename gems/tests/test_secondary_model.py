
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from argparse import Namespace
from os.path import join
from gems.secondary_model.sec_met_rxn_generation import (
        get_cluster_location,
        get_cluster_info_from_seq_record,
        get_cluster_domain,
        get_cluster_module
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

