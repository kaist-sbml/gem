
# Copyright 2017 BioInformatics Research Center, KAIST
# Copyright 2017 Novo Nordisk Foundation Center for Biosustainability, DTU

from os.path import join
from gems.secondary_model.sec_met_rxn_generation import (
        get_cluster_info_from_seq_record,
        get_cluster_domain
        )

class TestSecondary_model:
    """Test functions in gems.secondary_model"""

    def test_get_cluster_info_from_seq_record(self, seq_record):
        options = seq_record
        get_cluster_info_from_seq_record(options)

        # Number of genes for Cluster 3 in Streptomyces collinus Tu 365
        assert len(options.cluster_info_dict) == 22

        for locustag in options.cluster_info_dict.keys():
            assert 'NRPS/PKS Domain' in str(options.cluster_info_dict[locustag])

    def test_get_cluster_domain(self, options):
        get_cluster_domain(options)
