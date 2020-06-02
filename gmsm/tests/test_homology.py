
import pickle
from gmsm.config import load_config
from gmsm.homology import bidirect_blastp_analysis, blastp_utils
from os.path import abspath, dirname, isfile, join

class TestHomology:
    """Test functions in gmsm.homology"""

    # Streptomyces collinus Tu 365 : target, Streptomyces coelicolor A3(2) : template
    def test_make_blastDB(self, temp_fasta, target_fasta, options):

        options.outputfolder2 = 'gmsm/tests/data'
        options.target_fasta = target_fasta
        options.temp_fasta = temp_fasta
        
        blastp_utils.make_blastDB(options)
        assert isfile('%s/targetBlastDB.dmnd' %options.outputfolder2) == True
        assert isfile('%s/tempBlastDB.dmnd' %options.outputfolder2) == True
    
    
    # Streptomyces collinus Tu 365 : target, Streptomyces coelicolor A3(2) : template    
    def test_parseBlaspResults(self, inputFile_parseBlaspResults, outputFile_parseBlaspResults):
        
        parseBlaspResults = blastp_utils.parseBlaspResults(inputFile_parseBlaspResults, outputFile_parseBlaspResults)
        
        assert type(parseBlaspResults) == dict
        assert {'score': 2567.0, 'query_locusTag': 'B446_RS26080', 'evalue': '2.3e-291', 'identity': 93.2, 'length': 531, 'db_locusTag': 'SCO5535'} in parseBlaspResults.values()
        
    
    # Streptomyces collinus Tu 365 : target, Streptomyces coelicolor A3(2) : template
    def test_makeBestHits_dict(self, inputFile_makeBestHits_dict):
        
        bestHits_dict = blastp_utils.makeBestHits_dict(inputFile_makeBestHits_dict)
        
        assert type(bestHits_dict) == dict
        assert bestHits_dict['B446_RS28295'] == ['SCO6005']
        assert bestHits_dict['B446_RS18115'] == ['SCO4142']
    
    
    # Streptomyces coelicolor A3(2) : target, Streptomyces coelicolor A3(2) : template
    def test_getBBH(self, bestHits_dict1, bestHits_dict2, options):
        
        blastp_utils.getBBH(bestHits_dict1, bestHits_dict2, options)
        
        assert 'SCO5753' in options.targetBBH_list
        assert 'SCO3014' in options.targetBBH_list
        assert options.temp_target_BBH_dict['SCO1667'] == ['SCO1667', 'SCO5157'] or ['SCO5157', 'SCO1667']
        assert options.temp_target_BBH_dict['SCO7407'] == ['SCO7407', 'SCO5689'] or ['SCO5689', 'SCO7407']
        
    
    # Streptomyces collinus Tu 365 : target
    def test_get_nonBBH(self, options):
        
        options.targetGenome_locusTag_ec_dict = {'B446_23835':['4.1.1.45', '3.5.2.3']}
        options.targetBBH_list = ['B446_27575']
        
        blastp_utils.get_nonBBH(options, options)
        
        assert 'B446_23835' in options.nonBBH_list


    # get_homologs function works for executing functions of blastp_utils which already have test functions
    # Therefore, test_get_homologs function only checks the availability of get_homologs function
    def test_get_homologs(self, temp_fasta, target_fasta, options):
        
        options.outputfolder2 = 'gmsm/tests/data'
        options.target_fasta = target_fasta
        options.temp_fasta = temp_fasta
        options.targetGenome_locusTag_ec_dict = {'B446_23835':['4.1.1.45', '3.5.2.3']}
        options.targetBBH_list = ['B446_27575']
        
        bidirect_blastp_analysis.get_homologs(options, options)