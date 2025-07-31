#**GMSM**

***G***enome-scale metabolic ***M***odeling with ***S***econdary ***M***etabolism (GMSM) automatically generates secondary metabolite biosynthetic reactions in a genome-scale metabolic model (GEM) using antiSMASH output GenBank file. GMSM overall enables high-throughput modeling of both primary and secondary metabolism.

#Development
This project was initiated as a research collaboration between [Metabolic & Biomolecular Eng. Nat’l Research Laboratory (MBEL) & BioInformatics Research Center](http://mbel.kaist.ac.kr/) at KAIST and [Novo Nordisk Foundation Center for Biosustainability](http://www.biosustain.dtu.dk/english) at DTU.

#Current features
- Metabolic modeling for primary metabolism

- Metabolic modeling for secondary metabolism

#Installation

###Major dependencies

[biopython](http://biopython.org/wiki/Biopython)

[cobrapy](https://opencobra.github.io/cobrapy/) ([GitHub](https://github.com/opencobra/cobrapy); [Document](https://cobrapy.readthedocs.io/en/latest/))

###[DIAMOND](https://github.com/bbuchfink/diamond)

DIAMOND is a sequence alignment program for proteins, and is known for higher speed than BLAST.

conda install -c bioconda diamond

###Gurobi (optional, internal)
1. install Gurobi Optimizer for Python (via pip or Anaconda)

        conda install -c gurobi

2. Get a *Free Academic* license.
 
###Source
1. Clone the repository

    (HTTPS)

        git clone https://github.com/kaist-sbml/gem.git
    (SSH)

        git clone git@github.com:kaist-sbml/gem.git

2. Create and activate virtual environment

        conda create -n gmsm python=3.7
        conda activate gmsm

3. Install packages

        pip install -r requirements.txt
        
4. Clone large files in repository
        
        conda install -c conda-forge git-lfs
        git lfs install
        git lfs pull
        
5. Test GMSM

        tox

    
#Implementation
###General
- GMSM builds a GEM based on a template high-quality GEM. A default template GEM is the [high-quality GEM of Streptomyces coelicolor A3(2)](https://onlinelibrary.wiley.com/doi/full/10.1002/biot.201800180). Other template GEMs can be selected from the menu using `-m`.

- Select one or combination of modeling options using: `-e` (EC number annotation), `-p` (primary metabolism modeling) and/or `-s` (secondary metabolism modeling).
- Input file:

    Create an input directory at root of the GMSM directory.

    Input files can be a standard full GenBank file with sequences (recommended) or FASTA file.

    [antiSMASH](https://antismash.secondarymetabolites.org)-annotated GenBank file **MUST** be provided for secondary metabolism modeling.
    
    An EC number prediction file can be used with the `-e` option. The file format must consist of locus tags of the sequences and the predicted EC numbers, separated by a tab.

    EFICAz implementation and subcellular localizations (compartments) can be provided as additional inputs, with options `-E` and `-C`, respectively.

- Sample input files (available in `/gmsm/input/` in the source):

    `NC_021985.1.final_antismash8.gbk`: an output GenBank file of antiSMASH 8.0 (for Streptomyces collinus Tu 365)

    `NC_021985.1_deepec.txt`: an output file of DeepEC (for Streptomyces collinus Tu 365)

    `sample_compartment_info.txt`: a sample file containing subcellular localizations (compartments) for each locus tag

    `sample_eficaz_output.txt`: a sample output file of EFICAz

    `sample_input_ten_CDS.fasta`: a sample FASTA file having ten locus tags

    `sample_input_two_CDS.gb`: a sample GenBank file having two locus tags

- Output directory:

    Defining output directory is *optional*.

    Create an output directory at root of the GMSM directory.

    If output directory is not given, an output directory `output` is automatically generated at root of the GMSM repository. **Note**: New result files will override existing files in the default `output` directory.

- User's computer should be connected to the internet while modeling primary metabolism as GMSM accesses [KEGG](http://www.kegg.jp/kegg/rest/) to retrieve new reactions.
  
- For more information:

        python run_gmsm.py -h

###Examples

- Run modeling of primary metabolism (~30 min). This run will create the primary metabolism model necessary for secondary metabolism modeling.

        python run_gmsm.py -i input/NC_021985.1_antismash8.gbk -p -d

- Run modeling of secondary metabolism (only with antiSMASH output GenBank file, ~3 min). A GMSM-derived primary metabolism model should be available in an automatically generated folder `/output/3_primary_metabolic_model/`.

        python run_gmsm.py -i input/NC_021985.1_antismash8.gbk -s -d

- Run modeling of primary and secondary metabolism (~40 min).

        python run_gmsm.py -i input/NC_021985.1_antismash8.gbk -p -s -d

- Run modeling of primary metabolism using EC number prediction file and secondary metabolism (~40 min). 

        python run_gmsm.py -i input/NC_021985.1.final_antismash4.gbk -p -e input/NC_021985.1_deepec.txt -s -d

Note: Option `-d` is for displaying debugging statements during program running.

#Model refinement
Model draft created by GMSM should be refined to ensure its quality. Output files with prefix `rmc_` provide starting points for manual curation. `rmc_` stands for 'resource for manual curation'.

#Publication

