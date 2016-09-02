python ../NeighbourDist.py -v -l 10000 -g 0,100,0,0.7 -d 0,0.7,0,2500 heavy_n_junction_35000.fasta dist_10k_

python ../ClusterSeqs.py heavy_n_junction_35000.fasta heavy_n_junction_35000_clstr.clstr timeline.txt 0.15 "A1,A2,A3,TP,TS,TLN" -l 1000 -v

python ../ClusterGraph.py heavy_n_junction.fasta heavy_n_junction_clstr.clstr sample_colours.txt heavy S orange 0.15 -l 500

python ../PlotGermline.py kappa_n_A1_1000.txt,kappa_n_A1_1000.txt,kappa_n_A1_1000.txt,kappa_n_A1_1000.txt -t "A1,A2,A3" "V-GENE and allele" G -d "size=" -c 2 -y 15000 -b green,red,blue -g -gh -gv 3 -l 12 -w 0.8 -s kappa_germlines.pdf

python ../Spectratype.py kappa_aa_A1_1000.txt,kappa_aa_A2_1000.txt,kappa_aa_A3_1000.txt -t  "A1,A2,A3" -c 1 -u -xmin 5 -xmax 20 -g -w 0.8 -gv 3 -b red,yellow,green -s kappa_CDR3_length.pdf

python ../AbIdentity.py sample_heavy_n.txt H122 sample_intermediates_n.txt "M93173|IGHV1S40*01|Oryctolagus" IGHV1S40_n.txt intermediates.csv

python ../PlotIdentity.py H1-36_repertoire.csv -p "H1-36_ab.csv,rv/m10/lr,H1-36_germline.csv,b+/m20/lb,H1-36_family.csv,ro/m4" -my 85 -mx 20 -c autumn_r -g white -b -s identity.pdf