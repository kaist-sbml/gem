
import warnings
from argparse import Namespace
from os.path import join
from gmsm.secondary_model.sec_met_rxn_generation import (
        get_cluster_location,
        get_cluster_product,
        get_cluster_info_from_seq_record,
        get_cluster_monomers,
        get_all_metab_coeff,
        add_sec_met_rxn
        )

warnings.filterwarnings("ignore")

class TestSecondary_model:
    """Test functions in gmsm.secondary_model"""

    def test_get_cluster_info_from_seq_record(self, seq_record, options):

        # Hybrid cluster: nrps-transatpks-t1pks
        # locations: 341017 - 503094
        # Kirromycin biosynthetic gene cluster (81% of genes show similarity)
        cluster_nr = 3

        get_cluster_location(seq_record, cluster_nr, options)
        get_cluster_info_from_seq_record(seq_record, options)

        # Number of genes for Cluster 3 in Streptomyces collinus Tu 365
        assert len(options.cluster_info_dict) == 22

        for locustag in options.cluster_info_dict.keys():
            assert 'NRPS/PKS Domain' in str(options.cluster_info_dict[locustag])


    def test_get_cluster_product(self, seq_record, options):

        cluster_nr = 3

        get_cluster_location(seq_record, cluster_nr, options)
        get_cluster_product(seq_record, cluster_nr, options)

        assert options.product == 'Cluster03_nrps_t1pks_transatpks'


    def test_get_cluster_monomers(self, seq_record, options):

        cluster_nr = 3

        get_cluster_location(seq_record, cluster_nr, options)
        get_cluster_info_from_seq_record(seq_record, options)
        get_cluster_monomers(options)

        assert len(options.locustag_monomer_dict.keys()) == 17

        #Consensus monomers
        assert options.locustag_monomer_dict['B446_01480_M0'][3] == 'nrp'
        assert options.locustag_monomer_dict['B446_01480_M1'][3] == 'ser'
        assert options.locustag_monomer_dict['B446_01670_M0'][2] == 'ccmmal'


    def test_get_all_metab_coeff_antismash3(self, options):

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

        assert len(options.metab_coeff_dict) == 10
        assert options.metab_coeff_dict['mmcoa__R'] == -2
        assert options.metab_coeff_dict['ser__L'] == -3
        assert options.metab_coeff_dict['val__L'] == -2
        assert options.metab_coeff_dict['Cluster03_nrps_t1pks_transatpks'] == 1


    def test_get_all_metab_coeff_antismash4(self, options):

        locustag_monomer_dict = {
                'B446_RS01450_M0': ['no_call', 'N/A', 'trp', 'NA-n/a', 'no_call'],
                'B446_RS01450_M1': ['dpr|ser', 'ser', 'ser', 'NA-n/a', 'ser'],
                'B446_RS01455_M0': ['val', 'val', 'val', 'NA-n/a', 'no_call'],
                'B446_RS01455_M1': ['no_call', 'N/A', 'gly', 'NA-n/a', 'no_call'],
                'B446_RS01455_M2': ['aoh-gly|gly', 'gly', 'gly', 'NA-n/a', 'gly'],
                'B446_RS01495_M0': ['no_call', 'N/A', 'gly', 'NA-n/a', 'no_call'],
                'B446_RS01500_M0': ['dab', 'N/A', 'dab', 'NA-n/a', 'no_call'],
                'B446_RS01505_M0': ['val', 'val', 'val', 'NA-n/a', 'val'],
                'B446_RS01535_M0': ['athr|deoxy-thr|dhab|dht|ser|thr', 'thr', 'thr', 'NA-n/a', 'thr'],
                'B446_RS01535_M1': ['gln|ser', 'ser', 'ser', 'NA-n/a', 'ser'],
                'B446_RS01600_M0': ['mmal', 'mmal', 'mmal'],
                'B446_RS01620_M0': ['mal', 'mal', 'mal'],
                'B446_RS01620_M1': ['mal', 'mal', 'mal'],
                'B446_RS01625_M0': ['gly', 'gly', 'gly', 'NA-n/a', 'gly'],
                'B446_RS01635_M0': ['mmal', 'mmal', 'ccmmal'],
                'B446_RS01635_M1': ['mal', 'mal', 'mal'],
                'B446_RS01650_M0': ['gly', 'gly', 'gly', 'NA-n/a', 'gly']
                }

        options.locustag_monomer_dict = locustag_monomer_dict
        options.product = 'Cluster03_nrps_t1pks_transatpks'
        get_all_metab_coeff(options)

        assert len(options.metab_coeff_dict) == 9
        assert options.metab_coeff_dict['mmcoa__R'] == -2
        assert options.metab_coeff_dict['ser__L'] == -2
        assert options.metab_coeff_dict['val__L'] == -2
        assert options.metab_coeff_dict['Cluster03_nrps_t1pks_transatpks'] == 1


    def test_add_sec_met_rxn_cluster3(self,
            seq_record, sci_primary_model, mnxref, options):

        options.product = 'Cluster03_nrps_t1pks_transatpks'
        options.metab_coeff_dict = {
                '24dab': -1, 'leu__L': -1, 'mmcoa__R': -2, 'malcoa': -3,
                'arg__L': -1, 'ser__L': -2, 'thr__L': -1,
                'Cluster03_nrps_t1pks_transatpks': 1, 'val__L': -2, 'gly': -3,
                'ala__L': -1}

        # All the monomer for Cluster 3 are already present in model
        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {}

        assert 'Cluster03_nrps_t1pks_transatpks' not in sci_primary_model.reactions
        assert 'Cluster03_nrps_t1pks_transatpks_c' not in sci_primary_model.metabolites

        cluster_nr = 3
        get_cluster_location(seq_record, cluster_nr, options)
        get_cluster_info_from_seq_record(seq_record, options)
        model = add_sec_met_rxn(sci_primary_model, options)

        assert 'Cluster03_nrps_t1pks_transatpks' in model.reactions
        assert 'Cluster03_nrps_t1pks_transatpks_c' in model.metabolites


    def test_add_sec_met_rxn_cluster7(self,
            seq_record, sci_primary_model, mnxref, options):

        options.product = 'Cluster07_nrps_t1pks'
        options.metab_coeff_dict = {
                'mmcoa__R': -1, 'malcoa': -1, 'Cluster07_nrps_t1pks': 1, '23dhb': -1}

        # Following metabolite is absent in the model
        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {'MNXM455':['2,3-dihydroxybenzoate', 'C7H5O4']}

        assert 'Cluster07_nrps_t1pks' not in sci_primary_model.reactions
        assert '23dhb_c' not in sci_primary_model.metabolites
        assert 'Cluster07_nrps_t1pks_c' not in sci_primary_model.metabolites

        cluster_nr = 7
        get_cluster_location(seq_record, cluster_nr, options)
        get_cluster_info_from_seq_record(seq_record, options)
        model = add_sec_met_rxn(sci_primary_model, options)

        assert 'Cluster07_nrps_t1pks' in model.reactions
        assert '23dhb_c' in model.metabolites
        assert 'Cluster07_nrps_t1pks_c' in model.metabolites

