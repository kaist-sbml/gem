#Project
***GE***nome-scale metabolic ***M***odeling with ***S***econdary metabolism (GEMS) automatically generates secondary metabolite biosynthetic reactions in a genome-scale metabolic model (GEM) using antiSMASH output .gbk file. GEMS overall enables high-throughput modeling of primary and secondary metabolism.

#Development
This project was initiated as a research collaboration between [Metabolic & Biomolecular Eng. Nat’l Research Laboratory (MBEL) & BioInformatics Research Center](http://mbel.kaist.ac.kr/) at KAIST and [Novo Nordisk Foundation Center for Biosustainability](http://www.biosustain.dtu.dk/english) at DTU.

#Current features
- Homology analysis (bidirectional blastp hits)
- EC number annotation using [EFICAz](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html)
- Metabolic modeling for primary metabolism
- Metabolic modeling for secondary metabolism

#Installation
###Major dependencies
1. [`biopython`](http://biopython.org/wiki/Biopython) (version 1.68 tested)
2. [`blastp`](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/) and [`makeblastdb`](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/)
3. [`eficaz2.5`](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html) (versions 2.5 and 2.5.1 tested)
4. [`cobra`](https://opencobra.github.io/cobrapy/) (**MUST** be version 0.5.11 at the moment; [GitHub](https://github.com/opencobra/cobrapy); [Document](https://cobrapy.readthedocs.io/en/latest/))

###Gurobi (optional)
1. Create a symbolic link for the [gurobipy](http://www.gurobi.com/) installed in `root`. 

        ln -s /usr/local/lib/python2.7/dist-packages/gurobipy/ $HOME/gems/venv/lib/python2.7/site-packages/

2. Get a *Free Academic* license, and place it in a `GEMS` directory.

###Procedure
1. Clone the repository

    (HTTPS)

        git clone https://ehukim@bitbucket.org/ehukim/gems.git
    (SSH)

        git clone git@bitbucket.org:ehukim/gems.git

2. Create and activate virtual environment

        virtualenv venv
        source venv/bin/activate

3. Install packages

        pip install pip --upgrade
        pip install -r requirements.txt

    Installation of `zmq` and `numpy` using `requirements.txt` often causes an error. In this case, just do: `pip install zmq` and `pip install numpy`.

4. Test `GEMS`

    At the root of the repository,

        tox

5. `blastp` and `makeblastdb` for bidirectional blastp hits

    Get these executables from [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/). Preferably, place them in `venv/bin`.

    Make sure to get access to these executables using `chmod`.

6. EFICAz for EC number annotation

    Place `eficaz2.5` in a directory and set up `PATH` in `.bashrc`, e.g.:

    `export EFICAz25_PATH="/home/edhyunukkim/gems/venv/bin/EFICAz2.5.1/bin/`
    `export PATH="${PATH}:${EFICAz25_PATH}"`

    **Note**: Following statement causes a system error: `export PATH="/home/edhyunukkim/gems/venv/bin/EFICAz2.5.1/bin/"`.
    
#Implementation
###General
- Select one or combination of major implementation options: `-e`, `-p` and/or `-s`
- Input file:

    Create an input directory at root of the `GEMS` directory.

    Input file **MUST** be a standard full GenBank file with sequences.

    antiSMASH-annotated GenBank file **MUST** be provided for secondary metabolic modeling.

- Output directory:

    Defining output directory is *optional*.

    Create an output directory at root of the `GEMS` directory.

    If output directory is not given, result files are automatically stored in a directory `output` at root of the `GEMS` directory. **Note**: New result files will override existing files in the default `output` directory.
    
###Examples
- Run EC number annotation and modeling of primary and secondary metabolism

        run-gems -e -p -s -d -i input/NC_021985.1.final.gbk

- Run modeling of primary and secondary metabolism

        run-gems -p -s -d -i input/NC_021985.1.final.gbk

- Run modeling of primary metabolism

        run-gems -p -d -i input/NC_021985.1.final.gbk

- Run modeling of secondary metabolism

        run-gems -s -d -i input/NC_021985.1.final.gbk

- Run EC number annotation

        run-gems -e -d -i input/NC_021985.1.final.gbk

Note: Option `-d` is for displaying debugging statements during program running.

#Publication
TBD

#License
TBD
