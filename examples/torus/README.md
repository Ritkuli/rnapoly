#Generate primary structure
../../snacseq.py torus.snac -t 2

#Generate OxDNA files:
../../snac2ox.py -rg torus_seq.snac -o simulation/torus -xi 10
