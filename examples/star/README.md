#Generate primary structure
../../snacseq.py star.snac -t 2

#Generate OxDNA files:
../../snac2ox.py -rg star_seq.snac -o simulation/star
