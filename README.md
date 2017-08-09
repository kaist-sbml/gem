#**GMSM**
#Project
***G***nome-scale metabolic ***M***odeling with ***S***econdary ***M***etabolism (GMSM) automatically generates secondary metabolite biosynthetic reactions in a genome-scale metabolic model (GEM) using antiSMASH output GenBank file. GMSM overall enables high-throughput modeling of both primary and secondary metabolism.

#Development
This project was initiated as a research collaboration between [Metabolic & Biomolecular Eng. Nat’l Research Laboratory (MBEL) & BioInformatics Research Center](http://mbel.kaist.ac.kr/) at KAIST and [Novo Nordisk Foundation Center for Biosustainability](http://www.biosustain.dtu.dk/english) at DTU.

#Current features
- EC number annotation using [EFICAz](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html)
- Metabolic modeling for primary metabolism
- Metabolic modeling for secondary metabolism

#Installation
###Major dependencies
1. [biopython](http://biopython.org/wiki/Biopython)
2. [blastp](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
3. [eficaz2.5](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html)
4. [cobra](https://opencobra.github.io/cobrapy/) (version 0.6.2 or greater; [GitHub](https://github.com/opencobra/cobrapy); [Document](https://cobrapy.readthedocs.io/en/latest/))

###Gurobi (optional, internal)
1. Create a symbolic link for the [gurobipy](http://www.gurobi.com/) installed in `root`. 

        ln -s /usr/local/lib/python2.7/dist-packages/gurobipy/ $HOME/gmsm/venv/lib/python2.7/site-packages/

2. Get a *Free Academic* license.

###Docker
Docker image is available at [https://hub.docker.com/kaistmbel/gmsm]. Docker image contains all the major dependencies above and minimizes manutal setup. Currently light and full versions are available. All the Docker images are also tagged with GMSM versions.

1. *Light version*

    Light version has all the dependencies, including [blastp](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/), but **not** [eficaz2.5](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html).

    This version uses ~1.43 GB for disk space.

    Download the Docker image (~3 min):

        docker pull mbel/gmsm:0.4.5light

2. *Full version*

    Full version has all the dependencies, including [blastp](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/), **and** [eficaz2.5](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html).

    This version uses **~31 GB for disk space**.

    Download the Docker image **(~40 min)**:

        docker pull mbel/gmsm:0.4.5full
 
###Source
1. Clone the repository

    (HTTPS)

        git clone https://bitbucket.org/kaistmbel/gmsm.git
    (SSH)

        git clone git@bitbucket.org:kaistmbel/gmsm.git

2. Create and activate virtual environment

        virtualenv venv
        source venv/bin/activate

3. Install packages

        pip install pip --upgrade
        pip install -r requirements.txt

4. Test [GMSM](https://bitbucket.org/kaistmbel/gmsm)

        tox

5. [blastp](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) for bidirectional blastp hits

    Get these executables from [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/). Preferably, place them in `venv/bin`.

    Make sure to get access to these executables using `chmod`.

6. [EFICAz](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html) for EC number annotation (internal)

    Place [eficaz2.5](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html) in a directory and set up `PATH` in `.bashrc`, e.g.:

    	export EFICAz25_PATH="$HOME/gmsm/venv/bin/EFICAz2.5.1/bin/"
    	export PATH="${PATH}:${EFICAz25_PATH}"

    **Note**: Following statement causes a system error: `export PATH="$HOME/gmsm/venv/bin/EFICAz2.5.1/bin/"`.
    
#Implementation
###General
- Select one or combination of major implementation options: `-e`, `-p` and/or `-s`
- Input file:

    Create an input directory at root of the [GMSM](https://bitbucket.org/kaistmbel/gmsm) directory.

    Input files can be a standard full GenBank file with sequences (recommended) or FASTA file.

    antiSMASH-annotated GenBank file **MUST** be provided for secondary metabolism modeling.

    EFICAz output file and subcellular localizations (compartments) can be provided as additional inputs, with options `-E` and `-C`, respectively.

- Sample input files (available in `/gmsm/input/`):

    `NC_021985.1.final_antismash4.gbk`: an output GenBank file of antiSMASH 4.0

    `NC_021985.1.final_ec_antismash3.gbk`: an output GenBank file of antiSMASH 3.0 with full EC numbers via EFICAz

    `sample_compartment_info.txt`: a sample file containing subcellular localizations (compartments) for each locus tag

    `sample_eficaz_output.txt`: a sample output file of EFICAz

    `sample_input_ten_CDS.fasta`: a sample FASTA file having ten locus tags

    `sample_input_two_CDS.gb`: a sample GenBank file having two locus tags

- Output directory:

    Defining output directory is *optional*.

    Create an output directory at root of the [GMSM](https://bitbucket.org/kaistmbel/gmsm) directory.

    If output directory is not given, result files are automatically stored in a directory `output` at root of the [GMSM](https://bitbucket.org/kaistmbel/gmsm) directory. **Note**: New result files will override existing files in the default `output` directory.

- [GMSM](https://bitbucket.org/kaistmbel/gmsm) builds a GEM based on a template high-quality GEM. A default template GEM is the [high-quality GEM of Streptomyces coelicolor A3(2)](http://onlinelibrary.wiley.com/doi/10.1002/biot.201300539/abstract). Other template GEMs can be selected from the menu.
  
- For more information:

        run_gmsm.py -h
 
###Docker image
Upon download, run the Docker image:

        docker run --rm -it -v $HOME/users_input_dir:/gmsm/input  -v $HOME/users_output_dir:/gmsm/output mbel/gmsm:0.4.5full

- `users_input_dir`: User's defined directory where input data are stored.
- `users_output_dir`: User's defined directory where output data are stored.


###Examples
Following examples can be executed using both Docker image and source. However, `python` may need to be inserted, depending on user's system environment.

- Run EC number annotation (takes long time for a full genome, ~6-16 h)

        run_gmsm.py -i input/sample_input_two_CDS.gb -e -d

- Run modeling of primary metabolism

        run_gmsm.py -i input/NC_021985.1.final_antismash4.gbk -p -d

- Run modeling of primary metabolism using FASTA, EFICAz and compartment data

        run_gmsm.py -i input/sample_input_ten_CDS.fasta -m nsal -p -E input/sample_eficaz_output.txt -C input/sample_compartment_info.txt -d

- Run modeling of secondary metabolism (only with antiSMASH output GenBank file)

    Following command will not generate secondary metabolite biosynthetic reactions if a GMSM-derived primary metabolism model is not available in the designated folder (i.e., `3_primary_metabolic_model`).

        run_gmsm.py -i input/NC_021985.1.final_antismash4.gbk -s -d

- Run modeling of primary and secondary metabolism

        run_gmsm.py -i input/NC_021985.1.final_antismash4.gbk -p -s -d 

Note: Option `-d` is for displaying debugging statements during program running.

#Model refinement
Model draft created by GMSM should be refined to ensure its quality. Output files with prefix `rmc_` provide starting points for manual curation. `rmc_` stands for 'resource for manual curation'.

#Publication
Hyun Uk Kim, Jae Yong Ryu, Kyu-Sang Hwang, Tilmann Weber and Sang Yup Lee. ***GMSM***: ***G***enome-scale metabolic ***M***odeling with ***S***econdary ***M***etabolism.
