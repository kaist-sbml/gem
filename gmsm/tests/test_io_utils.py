
from Bio import SeqIO
from gmsm.io import io_utils
from os import makedirs
from os.path import isdir, isfile, join

class TestIo_utils:
    """Test functions in gmsm.io.io_utils"""

    
    def test_get_antismash_version_from_gbk(self, seq_record, options):
        
        io_utils.get_antismash_version_from_gbk(seq_record, options)
        
        assert options.anti_version == 4
       
    
    def test_get_features_from_gbk(self, seq_record, options):
        
        options.eficaz = False
        options.eficaz_file = False
        options.targetGenome_locusTag_aaSeq_dict = {}
        options.targetGenome_locusTag_ec_dict = {}
        options.targetGenome_locusTag_prod_dict = {}
        options.total_region = 0
        options.total_cluster = 0
        
        io_utils.get_features_from_gbk(seq_record, options, options)
        
        assert 'B446_17290' in options.targetGenome_locusTag_aaSeq_dict
        assert options.targetGenome_locusTag_prod_dict['B446_17290'] == 'integral membrane ATPase'
        assert options.targetGenome_locusTag_ec_dict['B446_24090'] == ['3.5.1.103']

 
    def test_get_features_from_fasta(self, input_fasta, options):
        
        options.targetGenome_locusTag_aaSeq_dict = {}
        options.targetGenome_locusTag_prod_dict = {}
        seq_records = list(SeqIO.parse(input_fasta, 'fasta'))
        
        for seq_record in seq_records:
            io_utils.get_features_from_fasta(seq_record, options)
        
        assert 'NSK_00001-RA' in options.targetGenome_locusTag_aaSeq_dict
        assert options.targetGenome_locusTag_prod_dict['NSK_00001-RA'] == \
               'NSK_00001-RA "Protein of unknown function" AED:0.27 eAED:0.27 QI:0|0|0|1|1|1|8|0|641'
        
        

    def test_get_target_fasta(self, options):
        
        options.targetGenome_locusTag_aaSeq_dict = \
        {'NSK_00003-RB' : 'MRRSLDDLFVPHGTSDLEAGALLYLLRLNAKTKTEVEDWLHQTIPCDLDPSRTISLPIRRDFLSGVQHLHGDL'}
        options.outputfolder2 = './tmp/2_blastp_results'
        if not isdir(options.outputfolder2):
            makedirs(options.outputfolder2)
            
        io_utils.get_target_fasta(options)
            
        assert options.target_fasta == './tmp/2_blastp_results/targetGenome_locusTag_aaSeq.fa'
        
        
    def test_get_temp_fasta(self, options):
        
        options.orgName = 'sco'
        
        io_utils.get_temp_fasta(options, options)
        
        assert options.input1 == './gmsm/io/data/input1/sco/'
        assert options.temp_fasta == './gmsm/io/data/input1/sco/tempModel_locusTag_aaSeq.fa'
        
        
      