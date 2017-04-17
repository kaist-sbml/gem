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
###Procedure
1. Clone the repository

    `git clone https://ehukim@bitbucket.org/ehukim/gems.git` (HTTPS) or 
    `git clone git@bitbucket.org:ehukim/gems.git` (SSH)

2. Create and activate virtual environment
    ```
    virtualenv venv
    source venv/bin/activate
    ```
3. Install packages
    ```
    pip install pip --upgrade
    pip install -r requirements.txt
    ```
    - Installation of `zmq` and `numpy` using `requirements.txt` often causes an error. In this case, just do: `pip install zmq` and `pip install numpy`
4. `tox` at the root of the repository to test `GEMS`
5. Place `blastp` and `makeblastdb` downloadable from [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/) preferably in `venv/bin`: for bidirectional blastp hits
6. Place `eficaz2.5` in a directory and set up `PATH` in `.bashrc`, e.g.:

    `export EFICAz25_PATH="/home/edhyunukkim/gems/venv/bin/EFICAz2.5.1/bin/`

    `export PATH="${PATH}:${EFICAz25_PATH}"`
    
###Notes
- Get access to these executables: `blastp` and `makeblastdb`.
- `export PATH="/home/edhyunukkim/gems/venv/bin/EFICAz2.5.1/bin/"` causes a system error.
- When using gurobi, make a symbolic link for the [gurobipy](http://www.gurobi.com/) installed in `root` and a get license for it: for optimization solvers.
   `ln -s /usr/local/lib/python2.7/dist-packages/gurobipy/ $HOME/gems/venv/lib/python2.7/site-packages/`

#Publication
TBD

#License
TBD
