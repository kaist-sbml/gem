
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from argparse import Namespace
from os.path import join
from gems.secondary_model.sec_met_rxn_generation import (
        get_cluster_location,
        get_cluster_product,
        get_cluster_info_from_seq_record,
        get_cluster_monomers,
        get_all_metab_coeff,
        add_sec_met_rxn
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

        assert len(options.metab_coeff_dict) == 11
        assert options.metab_coeff_dict['mmcoa_DASH_R'] == -2
        assert options.metab_coeff_dict['ser_DASH_L'] == -2
        assert options.metab_coeff_dict['val_DASH_L'] == -2
        assert options.metab_coeff_dict['Cluster03_nrps_t1pks_transatpks'] == 1


    def test_add_sec_met_rxn_cluster3(self,
            seq_record, sci_primary_model, mnxref, tmpdir, options):

        options.seq_record = seq_record
        options.product = 'Cluster03_nrps_t1pks_transatpks'
        options.metab_coeff_dict = {
                '24dab': -1, 'leu_DASH_L': -1, 'mmcoa_DASH_R': -2, 'malcoa': -3,
                'arg_DASH_L': -1, 'ser_DASH_L': -2, 'thr_DASH_L': -1,
                'Cluster03_nrps_t1pks_transatpks': 1, 'val_DASH_L': -2, 'gly': -3,
                'ala_DASH_L': -1}

        # All the monomer for Cluster 3 are already present in model
        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {}

        assert 'Cluster03_nrps_t1pks_transatpks' not in sci_primary_model.reactions
        assert 'Cluster03_nrps_t1pks_transatpks_c' not in sci_primary_model.metabolites

        cluster_nr = 3
        get_cluster_location(cluster_nr, options)
        get_cluster_info_from_seq_record(options)
        model = add_sec_met_rxn(sci_primary_model, options)

        assert 'Cluster03_nrps_t1pks_transatpks' in model.reactions
        assert 'Cluster03_nrps_t1pks_transatpks_c' in model.metabolites


    def test_add_sec_met_rxn_cluster7(self,
            seq_record, sci_primary_model, mnxref, tmpdir, options):

        options.seq_record = seq_record
        options.product = 'Cluster07_nrps_t1pks'
        options.metab_coeff_dict = {
                'mmcoa_DASH_R': -1, 'malcoa': -1, 'Cluster07_nrps_t1pks': 1, '23dhb': -1}

        # Following metabolite is absent in the model
        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {'MNXM455':['2,3-dihydroxybenzoate', 'C7H5O4']}

        assert 'Cluster07_nrps_t1pks' not in sci_primary_model.reactions
        assert '23dhb_c' not in sci_primary_model.metabolites
        assert 'Cluster07_nrps_t1pks_c' not in sci_primary_model.metabolites

        cluster_nr = 7
        get_cluster_location(cluster_nr, options)
        get_cluster_info_from_seq_record(options)
        model = add_sec_met_rxn(sci_primary_model, options)

        assert 'Cluster07_nrps_t1pks' in model.reactions
        assert '23dhb_c' in model.metabolites
        assert 'Cluster07_nrps_t1pks_c' in model.metabolites

