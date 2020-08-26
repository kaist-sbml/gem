
import shutil
from gmsm.config import load_config
from gmsm.io import input_file_manager
from os.path import isdir, join

class TestInput_file_manager:
    """Test functions in gmsm.io.input_file_manager"""

    
    def test_make_folder(self):
        
        folder = './tmp/test_folder'
        
        input_file_manager.make_folder(folder)
        
        assert isdir(folder) == True
        
        shutil.rmtree(folder)
    
    
    def test_setup_outputfolders(self, options):
        
        options.outputfolder = './tmp'
        options.eficaz = True
        options.pmr_generation = True
        options.smr_generation= True
        
        input_file_manager.setup_outputfolders(options, options)
        
        assert isdir(options.outputfolder1) == True
        assert isdir(options.outputfolder2) == True
        assert isdir(options.outputfolder3) == True
        assert isdir(options.outputfolder4) == True
        assert isdir(options.outputfolder5) == True
        assert isdir(options.outputfolder6) == True
        
        shutil.rmtree(options.outputfolder1)
        shutil.rmtree(options.outputfolder2)
        shutil.rmtree(options.outputfolder3)
        shutil.rmtree(options.outputfolder4)
        shutil.rmtree(options.outputfolder5)
        shutil.rmtree(options.outputfolder6)

    
    # show_input_options function is used for logging, so there is nothing to do assert
    def test_show_input_options(self, options):
        
        options.input = ''
        options.outputfolder = ''
        options.orgName = ''
        options.eficaz = ''
        options.pmr_generation = ''
        options.smr_generation = ''
        options.eficaz_file = ''
        options.comp = ''
        
        input_file_manager.show_input_options(options)
        
    
    def test_check_input_filetype(self, input_fasta, input_genbank, options):
        
        options.input = input_fasta
        
        return1 = input_file_manager.check_input_filetype(options)
        
        options.input = input_genbank
        
        return2 = input_file_manager.check_input_filetype(options)
        
        assert return1 == 'fasta'
        assert return2 == 'genbank'
        
    
    def test_get_target_genome_from_input(self, input_fasta, input_genbank, options):

        options.eficaz = False
        options.eficaz_file = False
        options.input = input_fasta
        input_file_manager.get_target_genome_from_input('fasta', options, options)

        assert len(options.targetGenome_locusTag_aaSeq_dict.keys()) == 10
        assert "Protein of unknown function" in options.targetGenome_locusTag_prod_dict['NSK_00001-RA']

        options.input = input_genbank
        input_file_manager.get_target_genome_from_input('genbank', options, options)

        assert options.seq_record_BGC_num_lists[0][0].id == 'NC_021985.1'
        assert options.total_cluster == 32
        assert options.targetGenome_locusTag_ec_dict['B446_RS01045'] == ['3.5.1.54']


    def test_get_eficaz_file(self, eficaz_file, options):

        options.eficaz_file = eficaz_file
        options.targetGenome_locusTag_ec_dict = {}
        input_file_manager.get_eficaz_file(options, options)

        assert len(options.targetGenome_locusTag_ec_dict) == 2
        assert options.targetGenome_locusTag_ec_dict['NSK_00005-RA'] == ['2.7.1.83']

        
    # get_fasta_files is used to execute functions of io_utils that have test functions, so there is nothing to do assert
    def test_get_fasta_files(self, options):
        
        options.targetGenome_locusTag_aaSeq_dict = None
        options.orgName = 'sco'
        options.temp_fasta = None
        
        input_file_manager.get_fasta_files(options, options)
        
    
    def test_get_pickles_prunPhase(self, options):
        
        options.input1 = 'gmsm/tests/data/'
        
        model = input_file_manager.get_pickles_prunPhase(options)
        
        assert options.tempModel_biggRxnid_locusTag_dict['GLNS'] == ['SCO2198', 'or', 'SCO2210']
        assert 'ALAD_L' in model.reactions
        assert 'ALAR' in model.reactions
        
        
    def test_get_pickles_augPhase(self, options):
        
        options.input1 = 'gmsm/tests/data/'
        
        input_file_manager.get_pickles_augPhase(options)
        
        assert options.bigg_mnxr_dict['UGMDDS2'] == 'MNXR105095'
        assert options.mnxm_compoundInfo_dict['MNXM44346'] == ["biphenyl-4,4'-diol", 'C12H10O2', 0]
        assert options.mnxr_kegg_dict['MNXR150233'] == ['R08385']
        assert 'MNXR120358' in options.mnxref.reactions
        assert options.template_exrxnid_flux_dict['EX_glc__D_e'] == -0.8
    

    def test_get_locustag_comp_dict(self, comp_file, options):

        options.comp = comp_file
        input_file_manager.get_locustag_comp_dict(options, options)

        assert len(options.locustag_comp_dict) == 8
        assert options.locustag_comp_dict['NSK_00004-RA'] == ['h']

