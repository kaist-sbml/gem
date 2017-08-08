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
1. [`biopython`](http://biopython.org/wiki/Biopython) (version 1.68 tested)
2. [`blastp`](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and [`makeblastdb`](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
3. [`eficaz2.5`](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html) (versions 2.5 and 2.5.1 tested)
4. [`cobra`](https://opencobra.github.io/cobrapy/) (version 0.6.2 or greater; [GitHub](https://github.com/opencobra/cobrapy); [Document](https://cobrapy.readthedocs.io/en/latest/))

###Gurobi (optional, internal)
1. Create a symbolic link for the [gurobipy](http://www.gurobi.com/) installed in `root`. 

        ln -s /usr/local/lib/python2.7/dist-packages/gurobipy/ $HOME/gmsm/venv/lib/python2.7/site-packages/

2. Get a *Free Academic* license, and place it in a `GMSM` directory.

###Procedure
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

4. Test GMSM

        tox

5. blastp and makeblastdb for bidirectional blastp hits

    Get these executables from [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/). Preferably, place them in `venv/bin`.

    Make sure to get access to these executables using `chmod`.

6. EFICAz for EC number annotation (internal)

    Place `eficaz2.5` in a directory and set up `PATH` in `.bashrc`, e.g.:

    	export EFICAz25_PATH="$HOME/gmsm/venv/bin/EFICAz2.5.1/bin/"
    	export PATH="${PATH}:${EFICAz25_PATH}"

    **Note**: Following statement causes a system error: `export PATH="$HOME/gmsm/venv/bin/EFICAz2.5.1/bin/"`.
    
#Implementation
###Docker image
To appear.

###General
- Select one or combination of major implementation options: `-e`, `-p` and/or `-s`
- Input file:

    Create an input directory at root of the `GMSM` directory.

    Input files can be a standard full GenBank file with sequences (recommended) or FASTA file.

    antiSMASH-annotated GenBank file **MUST** be provided for secondary metabolism modeling.

    EFICAz output file and subcellular localizations (compartments) can be provided as inputs, with options `-E` and `-C`, respectively.

- Output directory:

    Defining output directory is *optional*.

    Create an output directory at root of the `GMSM` directory.

    If output directory is not given, result files are automatically stored in a directory `output` at root of the `GMSM` directory. **Note**: New result files will override existing files in the default `output` directory.
    
###Examples
- Run EC number annotation and modeling of primary and secondary metabolism

        run_gmsm.py -e -p -s -d -i input/NC_021985.1.final.gbk

- Run modeling of primary and secondary metabolism

        run_gmsm.py -p -s -d -i input/NC_021985.1.final.gbk

- Run modeling of primary metabolism

        run_gmsm.py -p -d -i input/NC_021985.1.final.gbk

- Run modeling of secondary metabolism

        run_gmsm.py -s -d -i input/NC_021985.1.final.gbk

- Run EC number annotation

        run_gmsm.py -e -d -i input/NC_021985.1.final.gbk

Note: Option `-d` is for displaying debugging statements during program running.

#Model refinement
Model draft created by GMSM should be refined to ensure its quality. Output files with prefix `rmc_` provide starting points for manual curation. `rmc_` stands for 'resource for manual curation'.

#Publication
Hyun Uk Kim, Jae Yong Ryu, Kyu-Sang Hwang, Tilmann Weber and Sang Yup Lee. ***GMSM***: ***G***enome-scale metabolic ***M***odeling with ***S***econdary ***M***etabolism.
