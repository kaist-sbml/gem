
import os
import shutil
import time
import warnings
from os.path import join, isfile
from gmsm.io import output_file_manager

warnings.filterwarnings("ignore")

class TestOutput_file_manager:
    """Test functions in gmsm.io.output_file_manager"""
    
    # generate_outputs is used to execute functions of utils and output_file_manager that have test functions, so there is nothing to do assert
    def test_generate_outputs(self, sci_primary_model, options):

        options.rxnToRemove_dict = {}
        options.rxnToRemove_dict['PAPA160'] = '0'
        options.temp_target_BBH_dict = {'SCO4594' : ['SCO4594', 'SCO6269']}
        options.mnxr_to_add_list = ['MNXR107657']
        options.targetGenome_locusTag_ec_nonBBH_dict = {'SCO1202' : ['6.5.1.1']}
        options.rxnid_info_dict = {'R02672' : {'ENZYME': '4.1.1.7', 'DEFINITION': '4-Hydroxyphenylglyoxylate <=> 4-Hydroxybenzaldehyde + CO2', 'EQUATION': 'C03590 <=> C00633 + C00011', 'NAME': '4-Hydroxybenzoylformate carboxy-lyase', 'PATHWAY': 'rn00627 Aminobenzoate degradation'}}
        options.rxnid_locusTag_dict = {'R00375' : ['SCO1380', 'SCO3541', 'SCO4067', 'SCO2064', 'SCO3878', 'SCO1827']}
        options.rxn_newComp_fate_dict = {}
        options.verbose = False
        options.debug = True
        options.cpus = 8
        options.input = 'input/NC_021985.1.final_antismash4.gbk'
        options.outputfolder = 'output'
        options.orgName = 'sco'
        options.eficaz = False
        options.pmr_generation = True
        options.smr_generation = True
        options.eficaz_file = False
        options.comp = False
        options.adj_unique_nonprod_monomers_list = ['MNXM61686']
        options.outputfolder2 = './tmp/2_blastp_results'
        options.outputfolder3 = './tmp/3_primary_metabolic_model'
        options.outputfolder6 = './tmp/tmp_data_files'
        if not os.path.isdir(options.outputfolder2):
            os.makedirs(options.outputfolder2)
        if not os.path.isdir(options.outputfolder3):
            os.makedirs(options.outputfolder3)
        if not os.path.isdir(options.outputfolder6):
            os.makedirs(options.outputfolder6)
        start = time.time()
        runtime = time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start))
        
        output_file_manager.generate_outputs(options.outputfolder3, runtime, options, options, options, options, options, cobra_model = sci_primary_model)
        
    
    def test_get_model_reactions(self, sci_secondary_model, options):
    
        options.rxnToRemove_dict = {}
        options.rxnToRemove_dict['PAPA160'] = '0'
        folder = './tmp/4_complete_model'
        if not os.path.isdir(folder):
            os.makedirs(folder)     
        
        num_essen_rxn, num_kegg_rxn, num_bgc_rxn = output_file_manager.get_model_reactions(folder, options, cobra_model = sci_secondary_model)
        
        assert num_essen_rxn == 1
        assert num_kegg_rxn == 10
        assert num_bgc_rxn == 12
        
        
    def test_get_model_metabolites(self, sci_secondary_model, options):
        
        options.adj_unique_nonprod_monomers_list = ['MNXM61686']
        folder = './tmp/4_complete_model'
        if not os.path.isdir(folder):
            os.makedirs(folder)
        
        output_file_manager.get_model_metabolites(folder, sci_secondary_model, options)
        
        assert isfile(join(folder,'model_metabolites.txt')) == True
        assert isfile(join(folder,'rmc_metabolites_gapfilling_needed.txt')) == True
    
    
    def test_get_model_genes(self, tmpdir, sci_primary_model, options):

        options.orgName = 'sco'
        folder = './tmp'
        
        template_model_gene_list, duplicate_gene_list = \
                output_file_manager.get_model_genes(folder, sci_primary_model, options)

        assert 'B446_29745' in sci_primary_model.genes
        assert 'SCO0549' in sci_primary_model.genes
        assert 'SCO0549' in template_model_gene_list
        assert 'B446_25310' in duplicate_gene_list
        
        
    def test_get_summary_report(self, sci_secondary_model, options):
        
        options.verbose = False
        options.debug = True
        options.cpus = 8
        options.input = 'input/NC_021985.1.final_antismash4.gbk'
        options.outputfolder = 'output'
        options.orgName = 'sco'
        options.eficaz = False
        options.pmr_generation = True
        options.smr_generation = True
        options.eficaz_file = False
        options.comp = False
        options.adj_unique_nonprod_monomers_list = ['MNXM61686']
        num_essen_rxn = 1
        num_kegg_rxn = 10
        num_bgc_rxn = 12
        template_model_gene_list = range(1162)
        duplicate_gene_list = range(428)
        start = time.time()
        runtime = time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)) 
        folder = './tmp/4_complete_model'
        if not os.path.isdir(folder):
            os.makedirs(folder)
        
        output_file_manager.get_summary_report(folder, sci_secondary_model, runtime, num_essen_rxn, num_kegg_rxn, \
                                               num_bgc_rxn, template_model_gene_list, duplicate_gene_list, options, options)
        
        assert isfile(join(folder,'summary_report.txt')) == True
        

    def test_write_data_for_debug(self, options):
        
        options.temp_target_BBH_dict = {'SCO4594' : ['SCO4594', 'SCO6269']}
        options.mnxr_to_add_list = ['MNXR107657']
        options.targetGenome_locusTag_ec_nonBBH_dict = {'SCO1202' : ['6.5.1.1']}
        options.rxnid_info_dict = {'R02672' : {'ENZYME': '4.1.1.7', 'DEFINITION': '4-Hydroxyphenylglyoxylate <=> 4-Hydroxybenzaldehyde + CO2', 'EQUATION': 'C03590 <=> C00633 + C00011', 'NAME': '4-Hydroxybenzoylformate carboxy-lyase', 'PATHWAY': 'rn00627 Aminobenzoate degradation'}}
        options.rxnid_locusTag_dict = {'R00375' : ['SCO1380', 'SCO3541', 'SCO4067', 'SCO2064', 'SCO3878', 'SCO1827']}
        options.comp = False
        options.rxn_newComp_fate_dict = {}
        options.outputfolder2 = './tmp/2_blastp_results'
        options.outputfolder6 = './tmp/tmp_data_files'
        if not os.path.isdir(options.outputfolder2):
            os.makedirs(options.outputfolder2)
        if not os.path.isdir(options.outputfolder6):
            os.makedirs(options.outputfolder6)
        
        output_file_manager.write_data_for_debug(options, options, options, options)
        
        assert isfile(join(options.outputfolder2,'temp_target_BBH_dict.txt')) == True
        assert isfile(join(options.outputfolder6,'mnxr_to_add_list.txt')) == True
        assert isfile(join(options.outputfolder6,'targetGenome_locusTag_ec_nonBBH_dict.txt')) == True
        assert isfile(join(options.outputfolder6,'rxnid_info_dict.txt')) == True
        assert isfile(join(options.outputfolder6,'rxnid_locusTag_dict.txt')) == True

        
    def test_remove_tmp_model_files(self, options):
        
        options.outputfolder5 = './tmp/tmp_model_files'
        if not os.path.isdir(options.outputfolder5):
            os.makedirs(options.outputfolder5)
            
        output_file_manager.remove_tmp_model_files(options)
        
        assert isfile(options.outputfolder5) == False
        
        