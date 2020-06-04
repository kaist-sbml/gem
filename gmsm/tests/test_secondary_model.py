
import warnings
from argparse import Namespace
from os.path import join

from gmsm.config import load_config
from gmsm.secondary_model.run_secondary_modeling import (
        run_secondary_modeling, get_target_nonprod_monomers_for_gapfilling
        )
from gmsm.secondary_model.sec_met_rxn_generation import (
        get_region_location, get_cluster_location, get_region_info_from_seq_record,
        get_cluster_info_from_seq_record, get_region_product, get_cluster_product,
        get_region_monomers, get_cluster_monomers, get_all_metab_coeff, add_sec_met_rxn
        )

warnings.filterwarnings("ignore")

class TestSecondary_model:
    """Test functions in gmsm.secondary_model"""

    class TestRun_secondary_modeling:
    """Test functions in gmsm.secondary_model.run_secondary_modeling"""
    
    def test_run_secondary_modeling_anti5(self, seq_record_antismash5, sci_primary_model, options):
        
        load_config(options)
        options.seq_record_BGC_num_lists = [[seq_record_antismash5, 32]]
        options.temp_loc1 = 0
        options.anti_version = 5
        options.outputfolder5 = './tmp'
        model = run_secondary_modeling(sci_primary_model, options, options, options)
        
        get_target_nonprod_monomers_for_gapfilling(model, options, options, options)
        
        assert 'Region7_nrps_and_t1pks_and_bacteriocin' in model.reactions
        assert options.nonprod_sec_met_dict['Region10_t1pks_and_hgle_ks'] == ['malcoa', 'MNXM61686']
        assert options.adj_unique_nonprod_monomers_list == ['MNXM61686']
        
        
    def test_run_secondary_modeling_anti4(self, seq_record_antismash4, sci_primary_model, options):
        
        load_config(options)
        options.seq_record_BGC_num_lists = [[seq_record_antismash4, 32]]
        options.temp_loc1 = 0
        options.anti_version = 4
        options.outputfolder5 = './tmp'
        model = run_secondary_modeling(sci_primary_model, options, options, options)
        
        get_target_nonprod_monomers_for_gapfilling(model, options, options, options)
        
        assert 'Cluster7_t1pks_nrps' in model.reactions
        assert options.nonprod_sec_met_dict['Cluster23_t1pks'] == ['ala__L', 'malcoa', '23dhb']
        assert options.adj_unique_nonprod_monomers_list == ['23dhb', 'MNXM61686']
    
    def test_get_region_location(self, seq_record_antismash5, options):

        # Region 3 of NC_021985.1_antismash5.gbk file
        # Hybrid region: nrps-t1pks-transatpks
        # locations: 341093 - 489118
        # Kirromycin biosynthetic gene cluster (62% of genes show similarity)

        options.temp_loc1 = 207493
        get_region_location(seq_record_antismash5, options)

        assert options.region_loc1 == 341092
        assert options.region_loc2 == 489118
        assert options.temp_loc1 == 341092


    def test_get_cluster_location(self, seq_record_antismash4, options):

        # Cluster 3 of NC_021985.1_antismash4.gbk file
        # Hybrid cluster: nrps-t1pks-transatpks
        # locations: 341018 - 489943
        # Kirromycin biosynthetic gene cluster (79% of genes show similarity)

        options.temp_loc1 = 195827
        get_cluster_location(seq_record_antismash4, options)

        assert options.cluster_loc1 == 341017
        assert options.cluster_loc2 == 489943
        assert options.temp_loc1 == 341017


    def test_get_region_info_from_seq_record(self, seq_record_antismash5, options):

        options.temp_loc1 = 207493

        get_region_location(seq_record_antismash5, options)
        get_region_info_from_seq_record(seq_record_antismash5, options)

        # Number of genes involved in secondary metabolism for Region 3 in Streptomyces collinus Tu 365
        assert len(options.region_info_dict) == 26
        assert 'B446_RS01450' in str(options.region_info_dict.keys())
        assert 'AMP-binding' in str(options.region_info_dict['B446_RS01450'])
        assert 'B446_RS01600' in str(options.region_info_dict.keys())
        assert 'PKS_AT' in str(options.region_info_dict['B446_RS01600'])


    def test_get_cluster_info_from_seq_record(self, seq_record_antismash4, options):

        options.temp_loc1 = 195827

        get_cluster_location(seq_record_antismash4, options)
        get_cluster_info_from_seq_record(seq_record_antismash4, options)

        # Number of genes for involved in secondary metabolism Cluster 3 in Streptomyces collinus Tu 365
        assert len(options.cluster_info_dict) == 33
        assert 'B446_RS01450' in str(options.cluster_info_dict.keys())
        assert 'AMP-binding' in str(options.cluster_info_dict['B446_RS01450'])
        assert 'B446_RS01600' in str(options.cluster_info_dict.keys())
        assert 'PKS_AT' in str(options.cluster_info_dict['B446_RS01600'])


    def test_get_region_product(self, seq_record_antismash5, options):

        options.seq_record_BGC_num_lists = [[seq_record_antismash5, 32]]
        options.temp_loc1 = 207493
        order = 1
        region_nr = 3

        get_region_location(seq_record_antismash5, options)
        get_region_product(seq_record_antismash5, order, region_nr, options, options)

        assert options.product == 'Region3_nrps_and_transat_pks_and_t1pks'


    def test_get_cluster_product(self, seq_record_antismash4, options):

        options.temp_loc1 = 195827
        cluster_nr = 3

        get_cluster_location(seq_record_antismash4, options)
        get_cluster_product(seq_record_antismash4, cluster_nr, options)

        assert options.product == 'Cluster3_transatpks_t1pks_nrps'


    def test_get_region_monomers(self, seq_record_antismash5, options):

        options.temp_loc1 = 207493

        get_region_location(seq_record_antismash5, options)
        get_region_monomers(seq_record_antismash5, options)

        assert len(options.locustag_monomer_dict.keys()) == 17

        assert options.locustag_monomer_dict['B446_RS01450_M1'][0] == 'ser'
        assert options.locustag_monomer_dict['B446_RS01455_M2'][0] == 'gly'
        assert options.locustag_monomer_dict['B446_RS01600_M0'][0] == 'mmal'


    def test_get_cluster_monomers(self, seq_record_antismash4, options):

        options.temp_loc1 = 195827

        get_cluster_location(seq_record_antismash4, options)
        get_cluster_info_from_seq_record(seq_record_antismash4, options)
        get_cluster_monomers(options)

        assert len(options.locustag_monomer_dict.keys()) == 17

        #Consensus monomers
        assert options.locustag_monomer_dict['B446_RS01450_M0'][2] == 'trp'
        assert options.locustag_monomer_dict['B446_RS01450_M1'][4] == 'ser'
        assert options.locustag_monomer_dict['B446_RS01635_M0'][2] == 'ccmmal'


    #def test_get_all_metab_coeff_antismash3(self, options):

        #locustag_monomer_dict = {
                #'B446_01480_M0': ['orn,lys,arg', 'lys', 'leu', 'nrp'],
                #'B446_01480_M1': ['ser', 'ser', 'ser', 'ser'],
                #'B446_01485_M0': ['val', 'val', 'val', 'val'],
                #'B446_01485_M1': ['asp,asn,glu,gln,aad', 'N/A', 'ala', 'nrp'],
                #'B446_01485_M2': ['gly', 'gly', 'gly', 'gly'],
                #'B446_01525_M0': ['hydrophilic', 'arg', 'arg', 'arg'],
                #'B446_01530_M0': ['N/A', 'ser', 'dab', 'nrp'],
                #'B446_01535_M0': ['val', 'val', 'val', 'val'],
                #'B446_01565_M0': ['thr', 'thr', 'thr', 'thr'],
                #'B446_01565_M1': ['ser', 'ser', 'ser', 'ser'],
                #'B446_01635_M0': ['mmal', 'mmal', 'mmal'],
                #'B446_01655_M0': ['mal', 'mal', 'mal'],
                #'B446_01655_M1': ['mal', 'mal', 'mal'],
                #'B446_01660_M0': ['gly', 'gly', 'gly', 'gly'],
                #'B446_01670_M0': ['mmal', 'mmal', 'ccmmal'],
                #'B446_01670_M1': ['mal', 'mal', 'mal'],
                #'B446_01685_M0': ['gly', 'gly', 'gly', 'gly']}

        #options.locustag_monomer_dict = locustag_monomer_dict
        #options.product = 'Cluster03_nrps_t1pks_transatpks'
        #options.anti_version = 3
        #get_all_metab_coeff(options, options)

        #assert len(options.metab_coeff_dict) == 10
        #assert options.metab_coeff_dict['mmcoa__R'] == -2
        #assert options.metab_coeff_dict['ser__L'] == -3
        #assert options.metab_coeff_dict['val__L'] == -2
        #assert options.metab_coeff_dict['Cluster03_nrps_t1pks_transatpks'] == 1


    def test_get_all_metab_coeff_antismash5(self, options):

        locustag_monomer_dict = {
                'B446_RS01450_M0': ['X'],
                'B446_RS01450_M1': ['ser'],
                'B446_RS01455_M0': ['val'],
                'B446_RS01455_M1': ['X'],
                'B446_RS01455_M2': ['gly'],
                'B446_RS01495_M0': ['hydrophilic'],
                'B446_RS01500_M0': ['X'],
                'B446_RS01505_M0': ['val'],
                'B446_RS01535_M0': ['ser'],
                'B446_RS01535_M1': ['thr'],
                'B446_RS01600_M0': ['mmal', 'Methylmalonyl-CoA', 'mmal'],
                'B446_RS01620_M0': ['mal', 'Malonyl-CoA', 'mal'],
                'B446_RS01620_M1': ['mal', 'Malonyl-CoA', 'mal'],
                'B446_RS01625_M0': ['gly'],
                'B446_RS01635_M0': ['mal', 'Malonyl-CoA', 'mal'],
                'B446_RS01635_M1': ['mmal', 'Methylmalonyl-CoA', 'ccmmal'],
                'B446_RS01650_M0': ['gly']
                }

        options.locustag_monomer_dict = locustag_monomer_dict
        options.product = 'Region3_nrps_and_transat_pks_and_t1pks'
        options.anti_version = 5
        get_all_metab_coeff(options, options)

        assert len(options.metab_coeff_dict) == 7
        assert options.metab_coeff_dict['mmcoa__R'] == -2
        assert options.metab_coeff_dict['ser__L'] == -2
        assert options.metab_coeff_dict['thr__L'] == -1
        assert options.metab_coeff_dict['Region3_nrps_and_transat_pks_and_t1pks'] == 1


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
        options.product = 'Cluster3_nrps_t1pks_transatpks'
        options.anti_version = 4
        get_all_metab_coeff(options, options)

        assert len(options.metab_coeff_dict) == 9
        assert options.metab_coeff_dict['mmcoa__R'] == -2
        assert options.metab_coeff_dict['ser__L'] == -2
        assert options.metab_coeff_dict['val__L'] == -2
        assert options.metab_coeff_dict['Cluster3_nrps_t1pks_transatpks'] == 1


    def test_add_sec_met_rxn_region3(self,
            seq_record_antismash5, sci_primary_model, mnxref, options):

        options.anti_version = 5
        options.product = 'Region3_nrps_and_transat_pks_and_t1pks'
        options.metab_coeff_dict = {
                'gly': -3, 'mmcoa__R': -2, 'malcoa': -3, 'thr__L': -1, 
                'val__L': -2, 'ser__L': -2, 
                'Region3_nrps_and_transat_pks_and_t1pks': 1}

        #All the monomer for Region 3 are already present in model
        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {}

        assert 'Region3_nrps_and_transat_pks_and_t1pks' not in sci_primary_model.reactions
        assert 'Region3_nrps_and_transat_pks_and_t1pks_c' 

        options.temp_loc1 = 207493
        get_region_location(seq_record_antismash5, options)
        get_region_info_from_seq_record(seq_record_antismash5, options)
        model = add_sec_met_rxn(sci_primary_model, options, options)
        
        assert 'Region3_nrps_and_transat_pks_and_t1pks' in sci_primary_model.reactions
        assert 'Region3_nrps_and_transat_pks_and_t1pks_c' in sci_primary_model.metabolites

    def test_add_sec_met_rxn_region10(self,
            seq_record_antismash5, sci_primary_model, mnxref, options):

        options.anti_version = 5
        options.product = 'Region10_t1pks_and_hgle_ks'
        options.metab_coeff_dict = {'malcoa': -1, 'MNXM61686': -1, 
                'Region10_t1pks_and_hgle_ks': 1}

        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {}

        assert 'Region10_t1pks_and_hgle_ks' not in sci_primary_model.reactions
        assert 'Region10_t1pks_and_hgle_ks_c' not in sci_primary_model.metabolites

        options.temp_loc1 = 1500866
        get_region_location(seq_record_antismash5, options)
        get_region_info_from_seq_record(seq_record_antismash5, options)
        model = add_sec_met_rxn(sci_primary_model, options, options)

        assert 'Region10_t1pks_and_hgle_ks' in sci_primary_model.reactions
        assert 'Region10_t1pks_and_hgle_ks_c' in sci_primary_model.metabolites


    def test_add_sec_met_rxn_cluster3(self,
            seq_record_antismash4, sci_primary_model, mnxref, options):

        options.anti_version = 4
        options.product = 'Cluster3_nrps_t1pks_transatpks'
        options.metab_coeff_dict = {
                '24dab': -1, 'leu__L': -1, 'mmcoa__R': -2, 'malcoa': -3,
                'arg__L': -1, 'ser__L': -2, 'thr__L': -1,
                'Cluster3_nrps_t1pks_transatpks': 1, 'val__L': -2, 'gly': -3,
                'ala__L': -1}

        # All the monomer for Cluster 3 are already present in model
        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {}

        assert 'Cluster3_nrps_t1pks_transatpks' not in sci_primary_model.reactions
        assert 'Cluster3_nrps_t1pks_transatpks_c' not in sci_primary_model.metabolites

        options.temp_loc1 = 195827
        get_cluster_location(seq_record_antismash4, options)
        get_cluster_info_from_seq_record(seq_record_antismash4, options)
        model = add_sec_met_rxn(sci_primary_model, options, options)

        assert 'Cluster3_nrps_t1pks_transatpks' in model.reactions
        assert 'Cluster3_nrps_t1pks_transatpks_c' in model.metabolites


    def test_add_sec_met_rxn_cluster7(self,
            seq_record_antismash4, sci_primary_model, mnxref, options):

        options.anti_version = 4
        options.product = 'Cluster7_nrps_t1pks'
        options.metab_coeff_dict = {
                'mmcoa__R': -1, 'malcoa': -1, 'Cluster7_nrps_t1pks': 1, '23dhb': -1}

        # Following metabolite is absent in the model
        options.mnxref = mnxref
        options.mnxm_compoundInfo_dict = {'MNXM455':['2,3-dihydroxybenzoate', 'C7H5O4']}

        assert 'Cluster7_nrps_t1pks' not in sci_primary_model.reactions
        assert '23dhb_c' not in sci_primary_model.metabolites
        assert 'Cluster7_nrps_t1pks_c' not in sci_primary_model.metabolites

        options.temp_loc1 = 1073125
        get_cluster_location(seq_record_antismash4, options)
        get_cluster_info_from_seq_record(seq_record_antismash4, options)
        model = add_sec_met_rxn(sci_primary_model, options, options)

        assert 'Cluster7_nrps_t1pks' in model.reactions
        assert '23dhb_c' in model.metabolites
        assert 'Cluster7_nrps_t1pks_c' in model.metabolites

