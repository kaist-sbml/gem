#**GMSM**
#Project
***G***enome-scale metabolic ***M***odeling with ***S***econdary ***M***etabolism (GMSM) automatically generates secondary metabolite biosynthetic reactions in a genome-scale metabolic model (GEM) using antiSMASH output GenBank file. GMSM overall enables high-throughput modeling of both primary and secondary metabolism.

#Development
This project was initiated as a research collaboration between [Metabolic & Biomolecular Eng. Nat�l Research Laboratory (MBEL) & BioInformatics Research Center](http://mbel.kaist.ac.kr/) at KAIST and [Novo Nordisk Foundation Center for Biosustainability](http://www.biosustain.dtu.dk/english) at DTU.

#Current features
- Metabolic modeling for primary metabolism

- Metabolic modeling for secondary metabolism

#Installation

###Major dependencies

[biopython](http://biopython.org/wiki/Biopython)

[cobrapy](https://opencobra.github.io/cobrapy/) ([GitHub](https://github.com/opencobra/cobrapy); [Document](https://cobrapy.readthedocs.io/en/latest/))

###[DIAMOND] (https://github.com/bbuchfink/diamond)
DIAMOND is a sequence alignment program for proteins, and is known for higher speed than BLAST.

conda install -c bioconda diamond

###Gurobi (optional, internal)
1. install Gurobi Optimizer for Python (via pip or Anaconda)

        conda install -c gurobi

2. Get a *Free Academic* license.

###Docker (optional)
Docker image is available at https://hub.docker.com/r/mbelinsilico/gmsm. Docker image contains all the major dependencies above and minimizes manutal setup. Currently light and full versions are available. All the Docker images are also tagged with GMSM versions.

1. *Light version*

    Light version has all the dependencies, including [blastp](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/), but **not** [eficaz2.5](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html).

    This version uses ~1.43 GB for disk space.

    Download the Docker image (~3 min):

        docker pull mbelinsilico/gmsm:0.4.6light

2. *Full version*

    Full version has all the dependencies, including [blastp](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/), **and** [eficaz2.5](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html).

    This version uses **~31 GB for disk space**.

    Download the Docker image **(~40 min)**:

        docker pull mbelinsilico/gmsm:0.4.6full
 
###Source
1. Clone the repository

    (HTTPS)

        git clone https://bitbucket.org/kaistsystemsbiology/gmsm.git
    (SSH)

        git clone git@bitbucket.org:kaistsystemsbiology/gmsm.git

2. Create and activate virtual environment

        conda create -n gmsm python=3.7
        conda activate gmsm

3. Install packages

        pip install -r requirements.txt

4. Test [GMSM](https://bitbucket.org/kaistmbel/gmsm)

        tox

    
#Implementation
###General
- [GMSM](https://bitbucket.org/kaistmbel/gmsm) builds a GEM based on a template high-quality GEM. A default template GEM is the [high-quality GEM of Streptomyces coelicolor A3(2)](https://onlinelibrary.wiley.com/doi/full/10.1002/biot.201800180). Other template GEMs can be selected from the menu using `-m`.

- Select one or combination of modeling options using: `-e` (EC number annotation), `-p` (primary metabolism modeling) and/or `-s` (secondary metabolism modeling).
- Input file:

    Create an input directory at root of the [GMSM](https://bitbucket.org/kaistmbel/gmsm) directory.

    Input files can be a standard full GenBank file with sequences (recommended) or FASTA file.

    antiSMASH-annotated GenBank file **MUST** be provided for secondary metabolism modeling.

    EFICAz output file and subcellular localizations (compartments) can be provided as additional inputs, with options `-E` and `-C`, respectively.

- Sample input files (available in `/gmsm/input/` in the source):

    `NC_021985.1.final_antismash4.gbk`: an output GenBank file of antiSMASH 4.0 (for Streptomyces collinus Tu 365)

    `NC_021985.1.final_ec_antismash3.gbk`: an output GenBank file of antiSMASH 3.0 with full EC numbers via EFICAz (for Streptomyces collinus Tu 365)

    `sample_compartment_info.txt`: a sample file containing subcellular localizations (compartments) for each locus tag

    `sample_eficaz_output.txt`: a sample output file of EFICAz

    `sample_input_ten_CDS.fasta`: a sample FASTA file having ten locus tags

    `sample_input_two_CDS.gb`: a sample GenBank file having two locus tags

- Output directory:

    Defining output directory is *optional*.

    Create an output directory at root of the [GMSM](https://bitbucket.org/kaistmbel/gmsm) directory.

    If output directory is not given, an output directory `output` is automatically generated at root of the [GMSM](https://bitbucket.org/kaistmbel/gmsm) repository. **Note**: New result files will override existing files in the default `output` directory.

- User's computer should be connected to the internet while modeling primary metabolism as GMSM accesses [KEGG](http://www.kegg.jp/kegg/rest/) to retrieve new reactions.
  
- For more information:

        run_gmsm.py -h
 
###Docker image
Upon download, run the Docker image (full version):

        docker run --rm -it -v $HOME/users_input_dir:/gmsm/input -v $HOME/users_output_dir:/gmsm/output mbelinsilico/gmsm:0.4.6full

- `users_input_dir`: User's defined directory where input data are stored.
- `users_output_dir`: User's defined directory where output data are stored.
- To run the light version, replace `full` with `light` in the above command. 


###Examples
Following examples can be executed using both Docker image and source. However, `python` may need to be inserted at the beginning, depending on user's system environment. Running each example below takes a few minutes (~1-10 min) except for the last example.


- Run modeling of primary metabolism. This run will create the primary metabolism model necessary for secondary metabolism modeling.

        run_gmsm.py -i input/NC_021985.1.final_antismash4.gbk -p -d

- Run modeling of secondary metabolism (only with antiSMASH output GenBank file). A GMSM-derived primary metabolism model should be available in an automatically generated folder `/output/3_primary_metabolic_model/`.

        run_gmsm.py -i input/NC_021985.1.final_antismash4.gbk -s -d

- Run modeling of primary and secondary metabolism.

        run_gmsm.py -i input/NC_021985.1.final_antismash4.gbk -p -s -d

- Run modeling of primary metabolism using FASTA, EFICAz and compartment data.

        run_gmsm.py -i input/sample_input_ten_CDS.fasta -m nsal -p -E input/sample_eficaz_output.txt -C input/sample_compartment_info.txt -d

- Run modeling of primary and secondary metabolism (~30 min). Modeling using this input file takes much longer because this GenBank file has comprehensive EC number annotations from EFICAz, and thus has more reactions to be retrieved from [KEGG](http://www.kegg.jp/kegg/rest/) and added to the GEM draft.

        run_gmsm.py -i input/NC_021985.1.final_ec_antismash3.gbk -p -s -d

Note: Option `-d` is for displaying debugging statements during program running.

#Model refinement
Model draft created by GMSM should be refined to ensure its quality. Output files with prefix `rmc_` provide starting points for manual curation. `rmc_` stands for 'resource for manual curation'.

#Publication

