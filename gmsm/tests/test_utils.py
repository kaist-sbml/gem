
import warnings
from cobra import Reaction, Metabolite
from gmsm import utils
from gmsm.config import load_config
from os import remove
from os.path import isfile, join

warnings.filterwarnings("ignore")

class TestUtils:
    """Test functions in gmsm.utils"""
    
    def test_setup_logging(self, options):
        
        options.verbose = False
        options.debug = True
        options.outputfolder = './tmp'
        
        if isfile(join(options.outputfolder, 'gmsm.log')):
            remove(join(options.outputfolder, 'gmsm.log'))
        
        utils.setup_logging(options)
        
        assert isfile(join(options.outputfolder, 'gmsm.log')) == True
    
    
    def test_get_version(self):
    
        version = utils.get_version()
        
        assert version == '0.7.3'
        
        
    def test_get_git_log(self):
        
        output = utils.get_git_log()
        
        assert not output == None
        
        
    # check_input_options function is used for logging, so there is nothing to do assert    
    def test_check_input_options(self, options):
        
        options.input = ' '
        options.eficaz = False
        options.eficaz_file = ' '
        options.pmr_generation = True
        options.smr_generation = True
        options.comp = ' '
        
        utils.check_input_options(options)
        
        
    def test_locate_executable(self):
        
        name = 'diamond'
        
        output = utils.locate_executable(name)
        
        assert output == '/usr/bin/diamond'
        
        
    def test_execute(self):
        
        args = ['git', 'rev-parse', '--short', 'HEAD']
        
        out, err, retcode = utils.execute(args)
        
        assert not out == None
        assert err == ''
        assert retcode == 0
        
        
    def test_get_all_features_of_type(self, seq_record):
        
        types = 'CDS'
        test_feature1 = None
        test_feature2 = None
        
        for feature in seq_record.features:
            if feature.type == types:
                test_feature1 = feature
                break
            else:
                test_feature2 = feature
                
        features = utils.get_all_features_of_type(seq_record, types)
        
        assert test_feature1 in features
        assert not test_feature2 in features
    
    
    def test_get_cds_features(self, seq_record):
        
        test_feature1 = None
        test_feature2 = None
        
        for feature in seq_record.features:
            if feature.type == 'CDS':
                test_feature1 = feature
                break
            else:
                test_feature2 = feature
                
        features = utils.get_cds_features(seq_record)
        
        assert test_feature1 in features
        assert not test_feature2 in features        
        
        
    def test_get_gene_id(self, seq_record):
        
        test_feature = None
        
        for feature in seq_record.features:
            if feature.type == 'CDS':
                test_feature = feature
                break
        
        output = utils.get_gene_id(test_feature)
        
        assert output == 'B446_00005'

        
    def test_time_bomb(self, options):
        
        cache_file = join('./tmp', 'time_bomb.txt')
        
        if isfile(cache_file) == False:
            with open(cache_file, 'w') as f:
                f.write('test_time_bomb')
                
        load_config(options)
                
        utils.time_bomb(cache_file, options)
        
        assert isfile(cache_file) == True
    
    
    def test_get_keggid_from_mnxr(self, mnxr_kegg_dict, options):
        
        mnxr1 = 'MNXR109113'
        mnxr2 = 'MNXR111930'
        options.mnxr_kegg_dict = mnxr_kegg_dict
        options.rxnid_info_dict = {'R04783' : {'ENZYME': '3.2.1.23', 'DEFINITION': '3-Ketolactose + H2O <=> 3-Keto-beta-D-galactose + beta-D-Glucose', 'EQUATION': 'C05403 + C00001 <=> C05394 + C00221', 'NAME': '3-Ketolactose galactohydrolase'}}     
        
        kegg_id1 = utils.get_keggid_from_mnxr(mnxr1, options, options)
        kegg_id2 = utils.get_keggid_from_mnxr(mnxr2, options, options)
        
        assert kegg_id1 == 'R04783'
        assert kegg_id2 == 'R08385'
    
    
    def test_check_duplicate_rxn(self, sco_tmp_model, options):

        ala__L = Metabolite('ala_DASH_L_c')
        h2o = Metabolite('h2o_c')
        nad = Metabolite('nad_c')
        h = Metabolite('h_c')
        nadh = Metabolite('nadh_c')
        nh4 = Metabolite('nh4_c')
        pyr = Metabolite('pyr_c')

        rxn1 = Reaction('test_rxn1')
        rxn1.add_metabolites({
            ala__L: -1.0,
            h2o: -1.0,
            nad: -1.0,
            h: 1.0,
            nadh: 1.0,
            nh4: 1.0,
            pyr: 1.0
            })
        rxn2 = Reaction('test_rxn2')
        rxn2.add_metabolites({
            ala__L: -1.0,
            h2o: -1.0,
            nad: -1.0,
            nadh: 1.0,
            nh4: 1.0,
            pyr: 1.0
            })

        res1 = utils.check_duplicate_rxn(sco_tmp_model, rxn1)
        res2 = utils.check_duplicate_rxn(sco_tmp_model, rxn2)
        assert res1 == 'duplicate'
        assert res2 == 'unique'


    def test_compare_rxns(self, sco_tmp_model, options):

        #'ALAD_L': 'ala_DASH_L_c + h2o_c + nad_c --> h_c + nadh_c + nh4_c + pyr_c'
        rxn1 = sco_tmp_model.reactions[0]

        #'ALAR' = 'ala_DASH_L_c <=> ala_DASH_D_c'
        rxn2 = sco_tmp_model.reactions[1]

        ala__L = Metabolite('ala_DASH_L_c')
        h2o = Metabolite('h2o_c')
        nad = Metabolite('nad_c')
        h = Metabolite('h_c')
        nadh = Metabolite('nadh_c')
        nh4 = Metabolite('nh4_c')
        pyr = Metabolite('pyr_c')

        rxn3 = Reaction('test_rxn3')
        rxn3.add_metabolites({
            ala__L: -1.0,
            h2o: -1.0,
            nad: -1.0,
            h: 1.0,
            nadh: 1.0,
            nh4: 1.0,
            pyr: 1.0
            })

        res1 = utils.compare_rxns(rxn1, rxn2)
        res2 = utils.compare_rxns(rxn1, rxn3)
        
        assert res1 == 'different'
        assert res2 == 'same'

        
    def test_stablize_model(self, sci_primary_model):
        
        folder = './tmp'
        label = 'test'
        
        model = utils.stabilize_model(sci_primary_model, folder, label)
        
        assert len(model.reactions) == 1805
    
    
    def test_get_exrnxid_flux(self, sci_primary_model, sco_tmp_model_flux):
        
        target_exrxnid_flux_dict = utils.get_exrxnid_flux(sci_primary_model, sco_tmp_model_flux)
        
        assert 'Biomass_SCO' in target_exrxnid_flux_dict
        assert 'EX_o2_LPAREN_e_RPAREN_' in target_exrxnid_flux_dict
        assert 'EX_glc_LPAREN_e_RPAREN_' in target_exrxnid_flux_dict

        
    def test_check_exrxn_flux_direction(self, sco_tmp_model_flux, options):
        
        load_config(options)
        
        target_erxnid_flux_dict = {'EX_pi_LPAREN_e_RPAREN_': -0.11953457529767661, 'EX_co2_LPAREN_e_RPAREN_': 0.15033424301248413, 'EX_o2_LPAREN_e_RPAREN_': 1.4386047432726423, 'EX_nh4_LPAREN_e_RPAREN_': -1.0453876463751293, 'EX_h2o_LPAREN_e_RPAREN_': 0.0, 'EX_h_LPAREN_e_RPAREN_': 0.10504080912455104, 'Biomass_SCO': 0.11761613846320686, 'EX_glc_LPAREN_e_RPAREN_': -0.8}

        exrxn_flux_change_list = utils.check_exrxn_flux_direction(sco_tmp_model_flux, target_erxnid_flux_dict, options)
        
        assert 'F' in exrxn_flux_change_list
        
        
        
        
        
        
                              