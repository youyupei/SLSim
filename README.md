# Single-cell Long-read Simulator (SLSim)


## Dependency
* pandas
* numpy
* Biopython
* tqdm
* Badread (error_simulator.py)

#  Introduction
<img src="sim_wf.png" width="1000"/>

1. SLSsim randomly sample cDNA from transcript reference and attach 10x cell barcodes, UMI, poly(dT) and Template switching oligo (TSO) to form perfect reads (i.e. reads with no error.): `bin/perfect_read_generator.py`
2. SLSsim using the error model and qscore modle from Badread (one of the depencencies) to introduce error into the reads (`bin/error_simulator.py`). 

## Install
```
git clone https://github.com/youyupei/SLSim
cd SLSim
```

## Example code
1. Generate perfect reads (`template.fa`)
```
python3 bin/perfect_read_generator.py -r reference_transcript.fa -i bc_umi.csv -o template.fa
```

2. Simulate nanopore errors and create fastq file
```
python3 bin/error_simulator.py -t template.fa
```

## To do
`perfect_read_generator.py` is still single threading.
