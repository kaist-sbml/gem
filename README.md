#Project
This source automatically generates secondary metabolite biosynthetic reactions in a genome-scale metabolic model (GEM) using antiSMASH output .gbk file. The source also enables overall high-throughput metabolic modeling.

#Development
This project was initiated as a research collaboration between [Metabolic & Biomolecular Eng. Nat’l Research Laboratory (MBEL) & BioInformatics Research Center](http://mbel.kaist.ac.kr/) at KAIST and [Novo Nordisk Foundation Center for Biosustainability](http://www.biosustain.dtu.dk/english), DTU.

#Current features
* Homology analysis (bidirectional blastp hits)
* EC number annotation using [EFICAz](http://cssb.biology.gatech.edu/skolnick/webservice/EFICAz2/index.html)
* Metabolic modeling for primary metabolism
* Metabolic modeling for secondary metabolism

#Publication
TBD

#License
TBD

#Installation
##Run following in a terminal
```
git clone https://ehukim@bitbucket.org/ehukim/genome_analysis.git
```

##External tools
All the dependencies will be available in `requirements.txt`.
* [biopython](http://biopython.org/): for genome data handling
* [cobrapy](https://github.com/opencobra/cobrapy): for loading, editing and writing GEMs
* `blastp.exe` and `makeblastdb.exe` from [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/): for homology analysis
* [gurobi](http://www.gurobi.com/) or [glpk](https://www.gnu.org/software/glpk/): optimization solvers
* [libsbml](http://sbml.org/Main_Page): needed for cobrapy
* [numpy](http://www.numpy.org/): needed for cobrapy
* [scipy](http://scipy.org/): needed for cobrapy

**Note:** Make sure that user has access to `run_metabolicmodeling.py` (for convenience), `blastp.exe` and `makeblastdb.exe` (mandatory) in Linux system.

