
import warnings
from os.path import join
from gmsm.io.output_file_manager import get_model_genes

warnings.filterwarnings("ignore")

class TestOutput_file_manager:
    """Test functions in gmsm.io.output_file_manager"""

    def test_get_model_genes(self, tmpdir, sci_primary_model, options):

        folder = './tmp'
        options.orgName = 'sco'
        template_model_gene_list, duplicate_gene_list = \
                get_model_genes(folder, sci_primary_model, options, options)

        assert 'B446_29745' in sci_primary_model.genes
        assert 'SCO0549' in sci_primary_model.genes
        assert 'SCO0549' in template_model_gene_list
        assert 'B446_25310' in duplicate_gene_list
