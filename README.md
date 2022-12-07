# Simulator for single-cell Nanopore RNA sequencing reads

## Dependency

* Badread (add_error.py)
* pandas
* numpy
* Biopython
* tqdm

#  Single-cell Long-read Simulator (SLSim)

1. SLSsim randomly sample cDNA from transcript reference and attach 10x cell barcodes, UMI, poly(dT) and Template switching oligo (TSO) to form perfect reads (i.e. reads with no error.): `bin/perfect_read_generator.py`
2. SLSsim using the error model and qscore modle from Badread (one of the depencencies) to introduce error into the reads (`bin/error_simulator.py`). 
