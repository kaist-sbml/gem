#Project
This source automatically generates secondary metabolite biosynthetic reactions in a genome-scale metabolic model (GEM) using antiSMASH output .gbk file. The source also enables overall high-throughput metabolic modeling.

#Development
This project was initiated as a research collaboration between [Metabolic & Biomolecular Eng. Nat�l Research Laboratory (MBEL) & BioInformatics Research Center](http://mbel.kaist.ac.kr/) at KAIST and [Novo Nordisk Foundation Center for Biosustainability](http://www.biosustain.dtu.dk/english) at DTU.

#Current features
* Homology analysis (bidirectional blastp hits)
* EC number annotation using [EFICAz](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html)
* Metabolic modeling for primary metabolism
* Metabolic modeling for secondary metabolism

#Installation
###Run following in a terminal
1. Fork the repository `https://ehukim@bitbucket.org/ehukim/genome_analysis.git` and `git clone` it
2. `virtualenv venv` and activate `venv`
3. `pip install biopython` for [biopython](http://biopython.org/): for genome data handling
4. `pip install "cobra[all]"` for [cobrapy](https://github.com/opencobra/cobrapy): for loading, editing and writing GEMs
5. Place `blastp.exe` and `makeblastdb.exe` from [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/) preferably in `venv/bin`
6. Place `eficaz2.5` in a directory and set up `PATH` in `.bashrc`, e.g.,
```
export EFICAz25_PATH="/home/edhyunukkim/gems/venv/bin/EFICAz2.5.1/bin/"
export PATH="${PATH}:${EFICAz25_PATH}"
```
7. [gurobi](http://www.gurobi.com/) or [glpk](https://www.gnu.org/software/glpk/): optimization solvers

**Note:** Make sure that user has access to `run_gems.py` (for convenience), `blastp` and `makeblastdb` (mandatory) in Linux system.

#Publication
TBD

#License
TBD

