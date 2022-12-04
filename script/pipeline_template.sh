# Simulating ~90M  perfect reads
python3 path/to/bin/simulator.py --amp-rate 4

# File name output from simulator.py. Please modify it if a non-default filename was used
# when running simulator.py.
TEMPLATE=template.fa
# squigulator path
SQUIGULATOR=path/to/squigulator
BUTTERY_ENV=path/to/<buttery-eel environment>
GUPPY_BIN=path/to/guppy/bin


# Output file for simulated squiggle (BLOW5) format
BLOW5_FN=out_signal.blow5
OUT_FASTQ=sim.fastq
# Thread
THREAD=16

# simulate squiggle using squigulator
$SQUIGULATOR $TEMPLATE -o $BLOW5_FN  -t $THREAD --full-contigs

source $BUTTERY_ENV
buttery-eel -i $BLOW5_FN -o $OUT_FASTQ \
    -g  $GUPPY_BIN \
    --config dna_r9.4.1_450bps_sup.cfg   \
    --port 5557 --use_tcp   --device 'cuda:all' \
    --chunks_per_runner 1000 --slow5_threads 20 \
    --procs 10 -q 10

