#Generate primary structure
../../snacseq.py icosahedron.snac -t 2

#Generate OxDNA files:
../../snac2ox.py -rg icosahedron_seq.snac -o simulation/icosahedron
