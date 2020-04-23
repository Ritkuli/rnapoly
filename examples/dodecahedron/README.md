#Generate primary structure
../../snacseq.py dodecahedron.snac

#Generate OxDNA files:
../../snac2ox.py -rg dodecahedron_seq.snac -o simulation/dodecahedron -xi 10
